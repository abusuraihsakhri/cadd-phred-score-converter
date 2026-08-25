#!/usr/bin/env python3
"""
Clinical Correlation Module for CADD PHRED Score Converter.
Correlates CADD scores with clinical phenotype and genetic test results.
"""

from typing import Dict, Any, Optional, List


def correlate_clinical(cadd_phred: float, gene: str, hpo_terms: Optional[List[str]] = None,
                       previous_tests: Optional[List[str]] = None,
                       family_history: Optional[List[str]] = None) -> Dict[str, Any]:
    """Correlate CADD score with clinical context."""
    hpo_terms = hpo_terms or []
    previous_tests = previous_tests or []
    family_history = family_history or []

    clinical_score = cadd_phred
    clinical_factors = []

    gene_phenotype_map = {
        "BRCA1": ["hereditary_breast_ovarian_cancer", "breast_cancer", "ovarian_cancer"],
        "BRCA2": ["hereditary_breast_ovarian_cancer", "breast_cancer", "prostate_cancer"],
        "TP53": ["li_fraumeni_syndrome", "multiple_cancer_types"],
        "CFTR": ["cystic_fibrosis", "bronchiectasis"],
        "HBB": ["sickle_cell_disease", "beta_thalassemia"],
        "PKD1": ["polycystic_kidney_disease"],
    }

    phenotype_matches = gene_phenotype_map.get(gene, [])
    hpo_overlap = len(set(hpo_terms) & set(phenotype_matches))
    if hpo_overlap > 0:
        clinical_score *= 1.5
        clinical_factors.append(f"HPO overlap with {gene} phenotypes ({hpo_overlap} matches)")

    if "known_pathogenic_variant" in previous_tests:
        clinical_score *= 0.8
        clinical_factors.append("Known pathogenic variant in gene - consider compound heterozygosity")

    if "negative_segregation" in previous_tests:
        clinical_score *= 0.5
        clinical_factors.append("Family segregation negative")

    if any("autosomal_recessive" in fh for fh in family_history):
        clinical_score *= 1.2
        clinical_factors.append("Family history consistent with autosomal recessive pattern")

    if clinical_score >= 30:
        clinical_significance = "pathogenic_likely"
    elif clinical_score >= 20:
        clinical_significance = "uncertain_significant"
    elif clinical_score >= 10:
        clinical_significance = "uncertain"
    else:
        clinical_significance = "likely_benign"

    recommendations = []
    if clinical_significance in ("pathogenic_likely", "uncertain_significant"):
        recommendations.append("Report as variant of clinical significance")
        recommendations.append("Consider functional studies")
        if family_history:
            recommendations.append("Offer cascade genetic testing")

    return {
        "cadd_phred": cadd_phred,
        "gene": gene,
        "clinical_score": round(clinical_score, 2),
        "clinical_significance": clinical_significance,
        "clinical_factors": clinical_factors,
        "hpo_matches": hpo_overlap,
        "recommendations": recommendations,
    }


class ClinicalCorrelationAgent:
    """Sub-agent for clinical correlation."""

    def __init__(self):
        self.agent_name = "ClinicalCorrelationAgent"

    def evaluate(self, cadd_phred: float, gene: str, hpo_terms: Optional[List[str]] = None,
                 previous_tests: Optional[List[str]] = None,
                 family_history: Optional[List[str]] = None) -> Dict[str, Any]:
        """Evaluate clinical correlation."""
        result = correlate_clinical(cadd_phred, gene, hpo_terms, previous_tests, family_history)
        alerts = []

        if result["clinical_significance"] == "pathogenic_likely":
            alerts.append({
                "type": "CLINICALLY_SIGNIFICANT", "severity": "CRITICAL",
                "message": f"Variant in {gene} correlates with clinical phenotype. "
                           f"Score: {result['clinical_score']:.1f}.",
                "recommendation": "Report to ordering clinician. Consider variant reclassification."
            })

        return {"clinical_result": result, "alerts": alerts}
