#!/usr/bin/env python3
"""
Variant Annotation Enrichment for CADD PHRED Score Converter.
Enriches CADD scores with functional consequence, gene impact, and allele frequency data.
"""

from typing import Dict, Any, Optional, List
from dataclasses import dataclass


VARIANT_CONSEQUENCES = {
    "missense_variant": {"impact": "MODERATE", "category": "coding_sequence_variant"},
    "synonymous_variant": {"impact": "LOW", "category": "coding_sequence_variant"},
    "stop_gained": {"impact": "HIGH", "category": "coding_sequence_variant"},
    "stop_lost": {"impact": "HIGH", "category": "coding_sequence_variant"},
    "frameshift_variant": {"impact": "HIGH", "category": "coding_sequence_variant"},
    "splice_donor_variant": {"impact": "HIGH", "category": "splicing_variant"},
    "splice_acceptor_variant": {"impact": "HIGH", "category": "splicing_variant"},
    "intron_variant": {"impact": "MODIFIER", "category": "non_coding"},
    "upstream_gene_variant": {"impact": "MODIFIER", "category": "regulatory"},
    "downstream_gene_variant": {"impact": "MODIFIER", "category": "regulatory"},
    "intergenic_region": {"impact": "MODIFIER", "category": "non_coding"},
    "regulatory_region_variant": {"impact": "MODIFIER", "category": "regulatory"},
}


@dataclass
class VariantAnnotation:
    """Annotation for a single variant."""
    chromosome: str
    position: int
    ref: str
    alt: str
    consequence: str
    gene: str = ""
    protein_change: str = ""
    allele_frequency: float = 0.0


def annotate_variant(variant: VariantAnnotation, cadd_phred: float,
                     cadd_raw: float = 0.0) -> Dict[str, Any]:
    """Enrich variant with functional annotation."""
    consequence_info = VARIANT_CONSEQUENCES.get(variant.consequence,
                                                 {"impact": "UNKNOWN", "category": "unknown"})

    pathogenicity_score = cadd_phred
    if consequence_info["impact"] == "HIGH":
        pathogenicity_score *= 1.3
    elif consequence_info["impact"] == "MODERATE":
        pathogenicity_score *= 1.1

    if variant.allele_frequency > 0.05:
        population_status = "common_polymorphism"
        pathogenicity_score *= 0.3
    elif variant.allele_frequency > 0.01:
        population_status = "low_frequency"
        pathogenicity_score *= 0.7
    elif variant.allele_frequency > 0.001:
        population_status = "rare"
        pathogenicity_score *= 1.0
    else:
        population_status = "ultra_rare_or_novel"
        pathogenicity_score *= 1.2

    if pathogenicity_score >= 25:
        classification = "likely_pathogenic"
    elif pathogenicity_score >= 15:
        classification = "uncertain_significance"
    elif pathogenicity_score >= 10:
        classification = "likely_benign"
    else:
        classification = "benign"

    return {
        "variant": f"{variant.chromosome}:{variant.position}{variant.ref}>{variant.alt}",
        "gene": variant.gene,
        "consequence": variant.consequence,
        "impact": consequence_info["impact"],
        "protein_change": variant.protein_change,
        "allele_frequency": variant.allele_frequency,
        "population_status": population_status,
        "cadd_phred": cadd_phred,
        "cadd_raw": cadd_raw,
        "pathogenicity_score": round(pathogenicity_score, 2),
        "classification": classification,
    }


class VariantAnnotationAgent:
    """Sub-agent for variant annotation enrichment."""

    def __init__(self):
        self.agent_name = "VariantAnnotationAgent"

    def evaluate(self, variant: VariantAnnotation, cadd_phred: float,
                 cadd_raw: float = 0.0) -> Dict[str, Any]:
        """Evaluate variant annotation."""
        result = annotate_variant(variant, cadd_phred, cadd_raw)
        alerts = []

        if result["classification"] == "likely_pathogenic":
            alerts.append({
                "type": "LIKELY_PATHOGENIC_VARIANT", "severity": "CRITICAL",
                "message": f"{result['variant']} ({result['gene']}) classified as likely pathogenic "
                           f"(score: {result['pathogenicity_score']:.1f}).",
                "recommendation": "Report to clinician. Consider functional validation."
            })

        if result["population_status"] == "ultra_rare_or_novel" and cadd_phred > 20:
            alerts.append({
                "type": "RARE_HIGH_CADD", "severity": "WARNING",
                "message": "Ultra-rare variant with high CADD score.",
                "recommendation": "Strong candidate for pathogenicity. Family segregation analysis recommended."
            })

        return {"annotation": result, "alerts": alerts}
