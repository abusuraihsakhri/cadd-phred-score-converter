"""
CADD Converter Bridge Interface
===============================
Exports core CADD and PHRED conversion models and functions.
"""

from cadd_phred_converter import (
    TOTAL_GENOME_SNVS,
    MolecularConsequence,
    AcmgEvidenceTier,
    VariantLocation,
    VariantInput,
    CaddConversionResult,
    CaddPhredConverter,
    format_cadd_report,
)

__all__ = [
    "TOTAL_GENOME_SNVS",
    "MolecularConsequence",
    "AcmgEvidenceTier",
    "VariantLocation",
    "VariantInput",
    "CaddConversionResult",
    "CaddPhredConverter",
    "format_cadd_report",
]
