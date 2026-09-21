"""Utilities for interpreting CADD PHRED-scaled scores.

CADD raw scores are model outputs. Their mapping to PHRED-scaled scores depends on
ranking against a release-specific reference distribution, so this module does not
invent a universal raw-to-PHRED transform.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any, Dict, Optional

# CADD documentation describes ~8.6 billion possible GRCh37 substitutions.
# It is used only for an approximate rank/count display and is configurable.
DEFAULT_REFERENCE_SNVS = 8_600_000_000
TOTAL_GENOME_SNVS = DEFAULT_REFERENCE_SNVS  # backward-compatible alias


class MolecularConsequence(str, Enum):
    SYNONYMOUS = "synonymous_variant"
    MISSENSE = "missense_variant"
    STOP_GAINED = "stop_gained"
    STOP_LOST = "stop_lost"
    FRAMESHIFT = "frameshift_variant"
    SPLICE_DONOR = "splice_donor_variant"
    SPLICE_ACCEPTOR = "splice_acceptor_variant"
    SPLICE_REGION = "splice_region_variant"
    UTR_5 = "5_prime_UTR_variant"
    UTR_3 = "3_prime_UTR_variant"
    INTRONIC = "intron_variant"
    INTERGENIC = "intergenic_variant"
    NON_CODING = "non_coding_transcript_variant"


class AcmgEvidenceTier(str, Enum):
    """Retained for API compatibility; this converter does not assign ACMG evidence."""

    NOT_ASSESSED = "Not assessed from CADD alone"
    BENIGN_STRONG_BP4 = "Legacy: BP4 assignment disabled"
    INTERMEDIATE = "Legacy: intermediate category disabled"
    DELETERIOUS_MODERATE = "Legacy: PP3 moderate assignment disabled"
    PATHOGENIC_SUPPORTING_PP3 = "Legacy: PP3 assignment disabled"
    PATHOGENIC_STRONG_IN_SILICO = "Legacy: PP3 strong assignment disabled"


@dataclass
class VariantLocation:
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    gene_symbol: Optional[str] = None
    transcript_id: Optional[str] = None

    def __str__(self) -> str:
        chrom = self.chromosome[3:] if self.chromosome.lower().startswith("chr") else self.chromosome
        return f"chr{chrom}:{self.position}_{self.ref_allele}>{self.alt_allele}"


@dataclass
class VariantInput:
    variant_id: str
    location: Optional[VariantLocation] = None
    raw_score: Optional[float] = None
    phred_score: Optional[float] = None
    consequence: MolecularConsequence = MolecularConsequence.MISSENSE
    phylop_score: Optional[float] = None
    gerp_score: Optional[float] = None
    reference_size: int = DEFAULT_REFERENCE_SNVS

    def validate(self) -> None:
        if self.phred_score is None:
            if self.raw_score is not None:
                raise ValueError(
                    "A raw CADD score cannot be converted to PHRED without a "
                    "release-specific CADD ranking table. Provide the PHRED score "
                    "reported by CADD."
                )
            raise ValueError("phred_score is required.")
        if not math.isfinite(self.phred_score) or self.phred_score < 0.0:
            raise ValueError("phred_score must be a finite value >= 0.")
        if self.raw_score is not None and not math.isfinite(self.raw_score):
            raise ValueError("raw_score must be finite when supplied.")
        if self.reference_size <= 0:
            raise ValueError("reference_size must be a positive integer.")


@dataclass
class CaddConversionResult:
    variant_id: str
    raw_score: Optional[float]
    phred_score: float
    top_fraction_p: float
    top_percent: float
    percentile: float
    rank_from_top: int
    reference_size: int
    acmg_tier: AcmgEvidenceTier
    acmg_code: str
    is_pathogenic_predicted: Optional[bool]
    clinical_interpretation: str
    consequence: str
    conservation_summary: Dict[str, Any]

    @property
    def rank_among_genome_snvs(self) -> int:
        """Backward-compatible alias; rank 1 is the most deleterious end."""
        return self.rank_from_top

    @property
    def total_snvs_above(self) -> int:
        """Backward-compatible approximate count in the same top-tail set."""
        return self.rank_from_top

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data["acmg_tier"] = self.acmg_tier.value
        data["rank_among_genome_snvs"] = self.rank_from_top
        data["total_snvs_above"] = self.rank_from_top
        return data

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


class CaddPhredConverter:
    """Convert CADD PHRED rank units to tail fractions, percentiles and ranks."""

    @staticmethod
    def phred_from_tail_fraction(p: float) -> float:
        if not math.isfinite(p) or not (0.0 < p <= 1.0):
            raise ValueError("Tail fraction p must be finite and in (0, 1].")
        return -10.0 * math.log10(p)

    @staticmethod
    def tail_fraction_from_phred(phred: float) -> float:
        if not math.isfinite(phred) or phred < 0.0:
            raise ValueError("PHRED score must be finite and >= 0.")
        return 10.0 ** (-phred / 10.0)

    @classmethod
    def percentile_from_phred(cls, phred: float) -> float:
        return 100.0 * (1.0 - cls.tail_fraction_from_phred(phred))

    @classmethod
    def phred_from_percentile(cls, percentile: float) -> float:
        if not math.isfinite(percentile) or not (0.0 <= percentile < 100.0):
            raise ValueError("Percentile must be finite and in [0, 100).")
        return cls.phred_from_tail_fraction(1.0 - percentile / 100.0)

    @classmethod
    def rank_from_phred(cls, phred: float, reference_size: int = DEFAULT_REFERENCE_SNVS) -> int:
        if reference_size <= 0:
            raise ValueError("reference_size must be a positive integer.")
        return max(1, min(reference_size, math.ceil(reference_size * cls.tail_fraction_from_phred(phred))))

    @classmethod
    def raw_to_phred(cls, raw: float) -> float:
        raise NotImplementedError(
            "No universal raw-to-PHRED formula exists. Use a CADD release-specific "
            "ranking table or the PHRED value supplied by CADD."
        )

    @classmethod
    def phred_to_raw(cls, phred: float) -> float:
        raise NotImplementedError(
            "PHRED-to-raw inversion is release-specific and is not implemented."
        )

    @staticmethod
    def rank_context(phred: float) -> str:
        if phred >= 40:
            return "approximately within the top 0.01% of possible substitutions"
        if phred >= 30:
            return "approximately within the top 0.1% of possible substitutions"
        if phred >= 20:
            return "approximately within the top 1% of possible substitutions"
        if phred >= 10:
            return "approximately within the top 10% of possible substitutions"
        return "below the top 10% of possible substitutions"

    @classmethod
    def evaluate_variant(cls, variant: VariantInput) -> CaddConversionResult:
        variant.validate()
        assert variant.phred_score is not None
        phred = float(variant.phred_score)
        tail_fraction = cls.tail_fraction_from_phred(phred)
        percentile = cls.percentile_from_phred(phred)
        rank = cls.rank_from_phred(phred, variant.reference_size)

        conservation: Dict[str, Any] = {}
        if variant.phylop_score is not None:
            conservation["phylop"] = variant.phylop_score
        if variant.gerp_score is not None:
            conservation["gerp"] = variant.gerp_score

        interpretation = (
            f"CADD PHRED {phred:g} is {cls.rank_context(phred)}. "
            "This is a rank-based deleteriousness context, not a pathogenicity "
            "classification and not an automatic ACMG/AMP PP3 or BP4 assignment."
        )

        return CaddConversionResult(
            variant_id=variant.variant_id,
            raw_score=variant.raw_score,
            phred_score=phred,
            top_fraction_p=tail_fraction,
            top_percent=tail_fraction * 100.0,
            percentile=percentile,
            rank_from_top=rank,
            reference_size=variant.reference_size,
            acmg_tier=AcmgEvidenceTier.NOT_ASSESSED,
            acmg_code="NOT_ASSESSED",
            is_pathogenic_predicted=None,
            clinical_interpretation=interpretation,
            consequence=variant.consequence.value,
            conservation_summary=conservation,
        )


def format_cadd_report(res: CaddConversionResult) -> str:
    raw_line = "not supplied" if res.raw_score is None else f"{res.raw_score:g} (reported; not converted)"
    return "\n".join(
        [
            f"CADD PHRED RANK REPORT: {res.variant_id}",
            f"PHRED: {res.phred_score:g}",
            f"Raw CADD: {raw_line}",
            f"Top-tail fraction: {res.top_fraction_p:.8g}",
            f"Top-tail percent: {res.top_percent:.6g}%",
            f"Percentile: {res.percentile:.6f}%",
            f"Approx. rank from top: {res.rank_from_top:,} / {res.reference_size:,}",
            "ACMG/AMP evidence: not assessed from CADD alone",
            res.clinical_interpretation,
        ]
    )
