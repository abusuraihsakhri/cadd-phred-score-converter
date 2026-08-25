"""
Unit Test Suite for CADD Raw to PHRED Score Converter
=====================================================
Comprehensive verification of logarithmic PHRED conversions, genomic ranks,
empirical calibration curves, ACMG in-silico evidence codes, batch processing, and CLI.
"""

import csv
import json
import os
import tempfile
import unittest

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
import cli


class TestMathematicalPhredConversions(unittest.TestCase):
    """Test core logarithmic and percentile mathematical formulas."""

    def test_phred_from_tail_fractions_exact(self):
        self.assertAlmostEqual(CaddPhredConverter.phred_from_tail_fraction(1.0), 0.0, places=3)
        self.assertAlmostEqual(CaddPhredConverter.phred_from_tail_fraction(0.1), 10.0, places=3)
        self.assertAlmostEqual(CaddPhredConverter.phred_from_tail_fraction(0.01), 20.0, places=3)
        self.assertAlmostEqual(CaddPhredConverter.phred_from_tail_fraction(0.001), 30.0, places=3)
        self.assertAlmostEqual(CaddPhredConverter.phred_from_tail_fraction(0.0001), 40.0, places=3)

    def test_tail_fraction_from_phred_exact(self):
        self.assertAlmostEqual(CaddPhredConverter.tail_fraction_from_phred(0.0), 1.0, places=6)
        self.assertAlmostEqual(CaddPhredConverter.tail_fraction_from_phred(10.0), 0.1, places=6)
        self.assertAlmostEqual(CaddPhredConverter.tail_fraction_from_phred(20.0), 0.01, places=6)
        self.assertAlmostEqual(CaddPhredConverter.tail_fraction_from_phred(30.0), 0.001, places=6)

    def test_percentile_from_phred(self):
        self.assertAlmostEqual(CaddPhredConverter.percentile_from_phred(0.0), 0.0, places=2)
        self.assertAlmostEqual(CaddPhredConverter.percentile_from_phred(10.0), 90.0, places=2)
        self.assertAlmostEqual(CaddPhredConverter.percentile_from_phred(20.0), 99.0, places=2)
        self.assertAlmostEqual(CaddPhredConverter.percentile_from_phred(30.0), 99.9, places=2)
        self.assertAlmostEqual(CaddPhredConverter.percentile_from_phred(40.0), 99.99, places=2)

    def test_phred_from_percentile(self):
        self.assertAlmostEqual(CaddPhredConverter.phred_from_percentile(90.0), 10.0, places=1)
        self.assertAlmostEqual(CaddPhredConverter.phred_from_percentile(99.0), 20.0, places=1)
        self.assertAlmostEqual(CaddPhredConverter.phred_from_percentile(99.9), 30.0, places=1)

    def test_invalid_math_inputs_raise_errors(self):
        with self.assertRaises(ValueError):
            CaddPhredConverter.phred_from_tail_fraction(0.0)
        with self.assertRaises(ValueError):
            CaddPhredConverter.phred_from_tail_fraction(1.5)
        with self.assertRaises(ValueError):
            CaddPhredConverter.tail_fraction_from_phred(-5.0)
        with self.assertRaises(ValueError):
            CaddPhredConverter.phred_from_percentile(-1.0)
        with self.assertRaises(ValueError):
            CaddPhredConverter.phred_from_percentile(100.0)


class TestEmpiricalCalibration(unittest.TestCase):
    """Test empirical raw to PHRED calibration mappings."""

    def test_raw_to_phred_monotonicity(self):
        raws = [-1.0, 0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]
        phreds = [CaddPhredConverter.raw_to_phred(r) for r in raws]
        for i in range(len(phreds) - 1):
            self.assertLess(phreds[i], phreds[i + 1])

    def test_raw_to_phred_key_milestones(self):
        self.assertAlmostEqual(CaddPhredConverter.raw_to_phred(1.0), 10.0, places=1)
        self.assertAlmostEqual(CaddPhredConverter.raw_to_phred(3.0), 20.0, places=1)
        self.assertAlmostEqual(CaddPhredConverter.raw_to_phred(6.0), 30.0, places=1)

    def test_phred_to_raw_reversibility(self):
        for original_phred in [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]:
            approx_raw = CaddPhredConverter.phred_to_raw(original_phred)
            recovered_phred = CaddPhredConverter.raw_to_phred(approx_raw)
            self.assertAlmostEqual(original_phred, recovered_phred, delta=0.5)


class TestAcmgClinicalClassification(unittest.TestCase):
    """Test ACMG / AMP pathogenicity tiering and in silico evidence codes."""

    def test_benign_tier_bp4(self):
        var = VariantInput(variant_id="V-BENIGN", phred_score=6.5)
        res = CaddPhredConverter.evaluate_variant(var)
        self.assertEqual(res.acmg_tier, AcmgEvidenceTier.BENIGN_STRONG_BP4)
        self.assertEqual(res.acmg_code, "BP4")
        self.assertFalse(res.is_pathogenic_predicted)

    def test_intermediate_tier(self):
        var = VariantInput(variant_id="V-INTER", phred_score=12.5)
        res = CaddPhredConverter.evaluate_variant(var)
        self.assertEqual(res.acmg_tier, AcmgEvidenceTier.INTERMEDIATE)
        self.assertEqual(res.acmg_code, "NONE")
        self.assertFalse(res.is_pathogenic_predicted)

    def test_deleterious_moderate_tier(self):
        var = VariantInput(variant_id="V-MOD", phred_score=18.2)
        res = CaddPhredConverter.evaluate_variant(var)
        self.assertEqual(res.acmg_tier, AcmgEvidenceTier.DELETERIOUS_MODERATE)
        self.assertEqual(res.acmg_code, "PP3_MODERATE")
        self.assertTrue(res.is_pathogenic_predicted)

    def test_pathogenic_supporting_pp3(self):
        var = VariantInput(variant_id="V-PP3", phred_score=24.5)
        res = CaddPhredConverter.evaluate_variant(var)
        self.assertEqual(res.acmg_tier, AcmgEvidenceTier.PATHOGENIC_SUPPORTING_PP3)
        self.assertEqual(res.acmg_code, "PP3")
        self.assertTrue(res.is_pathogenic_predicted)
        self.assertGreaterEqual(res.percentile, 99.0)

    def test_pathogenic_strong_in_silico(self):
        var = VariantInput(variant_id="V-STRONG", phred_score=35.0)
        res = CaddPhredConverter.evaluate_variant(var)
        self.assertEqual(res.acmg_tier, AcmgEvidenceTier.PATHOGENIC_STRONG_IN_SILICO)
        self.assertEqual(res.acmg_code, "PP3_STRONG")
        self.assertTrue(res.is_pathogenic_predicted)
        self.assertGreaterEqual(res.percentile, 99.9)

    def test_genomic_rank_calculation(self):
        var = VariantInput(variant_id="V-RANK", phred_score=20.0)
        res = CaddPhredConverter.evaluate_variant(var)
        # PHRED 20 = top 1% -> total SNVs above = 88,000,000
        self.assertEqual(res.total_snvs_above, 88_000_000)
        self.assertEqual(res.rank_among_genome_snvs, TOTAL_GENOME_SNVS - 88_000_000 + 1)


class TestConservationAndFormatting(unittest.TestCase):
    """Test multi-metric conservation integration and report formatting."""

    def test_conservation_scores_populated(self):
        var = VariantInput(
            variant_id="V-CONS",
            phred_score=25.0,
            phylop_score=3.5,
            gerp_score=4.2
        )
        res = CaddPhredConverter.evaluate_variant(var)
        self.assertTrue(res.conservation_summary["phylop_conserved"])
        self.assertTrue(res.conservation_summary["gerp_conserved"])

    def test_formatting_and_json_serialization(self):
        var = VariantInput(variant_id="V-FMT", phred_score=22.0)
        res = CaddPhredConverter.evaluate_variant(var)
        # JSON
        d = res.to_dict()
        self.assertEqual(d["variant_id"], "V-FMT")
        js = res.to_json()
        parsed = json.loads(js)
        self.assertEqual(parsed["variant_id"], "V-FMT")
        # Text
        txt = format_cadd_report(res)
        self.assertIn("CADD & PHRED GENOMIC VARIANT DELETERIOUSNESS REPORT", txt)
        self.assertIn("V-FMT", txt)

    def test_variant_validation_empty_scores(self):
        var = VariantInput(variant_id="V-EMPTY")
        with self.assertRaises(ValueError):
            CaddPhredConverter.evaluate_variant(var)


class TestCLIExecution(unittest.TestCase):
    """Test CLI commands, demos, and batch CSV processing."""

    def test_cli_demos(self):
        self.assertEqual(cli.main(["--demo", "benign_synonymous"]), 0)
        self.assertEqual(cli.main(["--demo", "intermediate_missense"]), 0)
        self.assertEqual(cli.main(["--demo", "pathogenic_missense"]), 0)
        self.assertEqual(cli.main(["--demo", "severe_stopgain"]), 0)

    def test_cli_direct_args_json(self):
        ret = cli.main([
            "--variant-id", "CLI-VAR-01",
            "--phred", "28.5",
            "--consequence", "stop_gained",
            "--json"
        ])
        self.assertEqual(ret, 0)

    def test_cli_batch_csv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_in = os.path.join(tmpdir, "variants_in.csv")
            csv_out = os.path.join(tmpdir, "variants_out.csv")
            with open(csv_in, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["variant_id", "phred_score", "consequence"])
                writer.writeheader()
                writer.writerow({"variant_id": "V1", "phred_score": "5.0", "consequence": "synonymous_variant"})
                writer.writerow({"variant_id": "V2", "phred_score": "26.0", "consequence": "missense_variant"})

            ret = cli.main(["--batch-csv", csv_in, "--output", csv_out])
            self.assertEqual(ret, 0)
            self.assertTrue(os.path.exists(csv_out))


if __name__ == "__main__":
    unittest.main()
