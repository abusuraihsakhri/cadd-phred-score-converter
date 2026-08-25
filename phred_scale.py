#!/usr/bin/env python3
"""
CADD raw <-> PHRED scale conversion.

The CADD PHRED score of a variant is defined as

    PHRED = -10 * log10( p )

where p is the fraction of all ~8.8 billion possible single-nucleotide
substitutions whose raw CADD score is greater than or equal to the variant's.
Consequences implemented here:

  percentile_from_phred(phred) = 100 * (1 - 10**(-phred/10))
      PHRED 20 -> top 1% ; PHRED 30 -> top 0.1% ; PHRED 40 -> top 0.01%

  genome_count_above(phred) = 8.8e9 * 10**(-phred/10)
  rank_of_variant(phred)    = N - count_above + 1

Raw-to-PHRED conversion uses the exponential-tail model documented by the
CADD authors: under P(raw >= r) ~= 10**(-r/10), the PHRED value equals the
raw score (each +10 on either scale = 10x rarer). The legacy heuristic found
in older pipelines, PHRED_legacy = 4**raw / 10, is provided for compatibility
but diverges rapidly from the true scale and is clamped to [1, 99].
"""

import math
from typing import Dict

GENOME_SNV_COUNT = 8_800_000_000   # approximate possible SNVs in GRCh37/38


def phred_from_tail_fraction(p: float) -> float:
    """PHRED score for a rank fraction p in (0, 1]."""
    if not 0.0 < p <= 1.0:
        raise ValueError("p must be within (0, 1]")
    return -10.0 * math.log10(p)


def tail_fraction_from_phred(phred: float) -> float:
    return 10.0 ** (-phred / 10.0)


def percentile_from_phred(phred: float) -> float:
    """Genome-wide percentile rank of a PHRED score."""
    return 100.0 * (1.0 - tail_fraction_from_phred(phred))


def phred_from_percentile(percentile: float) -> float:
    if not 0.0 <= percentile < 100.0:
        raise ValueError("percentile must be in [0, 100)")
    return phred_from_tail_fraction(1.0 - percentile / 100.0)


def genome_count_above(phred: float) -> int:
    return max(0, round(GENOME_SNV_COUNT * tail_fraction_from_phred(phred)))


def rank_of_variant(phred: float) -> int:
    """1-based rank among all genome substitutions (higher score = lower rank)."""
    return GENOME_SNV_COUNT - genome_count_above(phred) + 1


def raw_to_phred_exponential_tail(raw: float) -> float:
    """Raw -> PHRED under the exponential-tail model (identity mapping)."""
    return float(raw)


def legacy_four_power_conversion(cadd_score: float) -> float:
    """Legacy heuristic PHRED = 4**cadd / 10, clamped to [1, 99].

    Kept only for compatibility with older annotation pipelines; it grows
    exponentially and does NOT reproduce the published PHRED scale beyond
    very small raw scores. Prefer raw_to_phred_exponential_tail.
    """
    val = (4.0 ** cadd_score) / 10.0
    return min(99.0, max(1.0, val))


def annotate_variant_score(raw: float = None, phred: float = None) -> Dict[str, float]:
    """Full conversion report from whichever scale was supplied."""
    if phred is None and raw is None:
        raise ValueError("provide raw or phred")
    if phred is None:
        phred = raw_to_phred_exponential_tail(raw)
    return {
        "cadd_raw": raw,
        "cadd_phred": round(phred, 3),
        "percentile": round(percentile_from_phred(phred), 4),
        "genome_variants_above": genome_count_above(phred),
        "rank_of_variant": rank_of_variant(phred),
        "legacy_heuristic_value": round(legacy_four_power_conversion(
            raw if raw is not None else phred), 2),
    }


if __name__ == "__main__":
    print(f"{'PHRED':>6} {'percentile':>12} {'#variants above':>17}")
    print("-" * 38)
    for p in (10, 15, 20, 25, 30, 35, 40):
        print(f"{p:>6} {percentile_from_phred(p):>11.3f}% "
              f"{genome_count_above(p):>17,}")

    print("\nVariant reports:")
    for label, kw in (("raw=20", {"raw": 20.0}),
                      ("raw=33", {"raw": 33.0}),
                      ("phred=28", {"phred": 28.0})):
        rep = annotate_variant_score(**kw)
        print(f"  {label}: PHRED={rep['cadd_phred']} pct={rep['percentile']}% "
              f"above={rep['genome_variants_above']:,} "
              f"legacy4^x={rep['legacy_heuristic_value']}")
