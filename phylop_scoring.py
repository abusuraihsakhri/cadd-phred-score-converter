#!/usr/bin/env python3
"""
Phylop scoring for CADD PHRED Score Converter.
Calculates phylop-based conservation scores for variant pathogenicity assessment.
"""

from typing import Dict, Any, Optional


PHYLOP_THRESHOLDS = {
    "highly_conserved": 2.0,
    "moderately_conserved": 1.0,
    "weakly_conserved": 0.5,
    "neutral": -0.5,
}


def calculate_phylop_score(phylop_score: float, gerp_score: float = 0.0,
                            conservation_percentile: float = 50.0) -> Dict[str, Any]:
    """Calculate phylop-based conservation assessment."""
    if phylop_score >= PHYLOP_THRESHOLDS["highly_conserved"]:
        conservation_level = "highly_conserved"
        pathogenicity_boost = 1.5
        interpretation = "Strong evolutionary conservation. Deleterious variants in this region are likely pathogenic."
    elif phylop_score >= PHYLOP_THRESHOLDS["moderately_conserved"]:
        conservation_level = "moderately_conserved"
        pathogenicity_boost = 1.2
        interpretation = "Moderate conservation. Variants may have functional impact."
    elif phylop_score >= PHYLOP_THRESHOLDS["weakly_conserved"]:
        conservation_level = "weakly_conserved"
        pathogenicity_boost = 1.0
        interpretation = "Weak conservation. Pathogenicity uncertain from conservation alone."
    elif phylop_score >= PHYLOP_THRESHOLDS["neutral"]:
        conservation_level = "near_neutral"
        pathogenicity_boost = 0.8
        interpretation = "Near-neutral evolution. Less constraint on this position."
    else:
        conservation_level = "accelerated"
        pathogenicity_boost = 0.5
        interpretation = "Accelerated evolution. May indicate positive selection."

    combined_score = phylop_score * pathogenicity_boost
    if gerp_score > 0:
        combined_score = (phylop_score * 0.6 + gerp_score * 0.4) * pathogenicity_boost

    return {
        "phylop_score": round(phylop_score, 3),
        "gerp_score": round(gerp_score, 3),
        "conservation_level": conservation_level,
        "pathogenicity_boost": pathogenicity_boost,
        "combined_conservation_score": round(combined_score, 3),
        "conservation_percentile": round(conservation_percentile, 1),
        "interpretation": interpretation,
    }


class PhylopAgent:
    """Sub-agent for phylop scoring."""

    def __init__(self):
        self.agent_name = "PhylopAgent"

    def evaluate(self, phylop_score: float, gerp_score: float = 0.0,
                 conservation_percentile: float = 50.0) -> Dict[str, Any]:
        """Evaluate phylop conservation."""
        result = calculate_phylop_score(phylop_score, gerp_score, conservation_percentile)
        alerts = []

        if result["conservation_level"] == "highly_conserved":
            alerts.append({
                "type": "HIGH_CONSERVATION", "severity": "WARNING",
                "message": f"Highly conserved position (phylop: {phylop_score:.3f}).",
                "recommendation": "Variants at this position are likely deleterious. Prioritize for review."
            })
        elif result["conservation_level"] == "accelerated":
            alerts.append({
                "type": "ACCELERATED_EVOLUTION", "severity": "ADVISORY",
                "message": f"Accelerated evolution detected (phylop: {phylop_score:.3f}).",
                "recommendation": "May indicate positive selection. Interpret with caution."
            })

        return {"phylop_result": result, "alerts": alerts}
