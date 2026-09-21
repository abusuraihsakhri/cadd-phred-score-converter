import csv
import json
import math
from pathlib import Path

import cli
from cadd_phred_converter import (
    DEFAULT_REFERENCE_SNVS,
    AcmgEvidenceTier,
    CaddPhredConverter,
    VariantInput,
)


def test_phred_tail_fraction_landmarks():
    for phred, fraction in [(0, 1), (10, .1), (20, .01), (30, .001), (40, .0001)]:
        assert CaddPhredConverter.tail_fraction_from_phred(phred) == pytest_approx(fraction)
        assert CaddPhredConverter.phred_from_tail_fraction(fraction) == pytest_approx(phred)


def pytest_approx(value, rel=1e-12, abs=1e-12):
    import pytest
    return pytest.approx(value, rel=rel, abs=abs)


def test_percentile_round_trip():
    for percentile in [0, 90, 99, 99.9, 99.99]:
        phred = CaddPhredConverter.phred_from_percentile(percentile)
        assert CaddPhredConverter.percentile_from_phred(phred) == pytest_approx(percentile, abs=1e-9)


def test_rank_is_from_most_deleterious_end():
    assert CaddPhredConverter.rank_from_phred(20, DEFAULT_REFERENCE_SNVS) == 86_000_000
    assert CaddPhredConverter.rank_from_phred(30, DEFAULT_REFERENCE_SNVS) == 8_600_000


def test_invalid_inputs():
    import pytest
    for value in [-1, math.inf, math.nan]:
        with pytest.raises(ValueError):
            CaddPhredConverter.tail_fraction_from_phred(value)
    with pytest.raises(ValueError):
        CaddPhredConverter.phred_from_tail_fraction(0)
    with pytest.raises(ValueError):
        CaddPhredConverter.phred_from_percentile(100)


def test_raw_conversion_is_explicitly_not_fabricated():
    import pytest
    with pytest.raises(NotImplementedError):
        CaddPhredConverter.raw_to_phred(3.0)
    with pytest.raises(NotImplementedError):
        CaddPhredConverter.phred_to_raw(20.0)


def test_variant_evaluation_does_not_assign_pathogenicity():
    result = CaddPhredConverter.evaluate_variant(VariantInput("v", raw_score=3.2, phred_score=25.3))
    assert result.raw_score == 3.2
    assert result.acmg_tier is AcmgEvidenceTier.NOT_ASSESSED
    assert result.acmg_code == "NOT_ASSESSED"
    assert result.is_pathogenic_predicted is None
    assert "not a pathogenicity classification" in result.clinical_interpretation


def test_raw_only_variant_is_rejected():
    import pytest
    with pytest.raises(ValueError, match="release-specific"):
        CaddPhredConverter.evaluate_variant(VariantInput("v", raw_score=3.0))


def test_cli_json(capsys):
    assert cli.main(["--phred", "20", "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["top_percent"] == pytest_approx(1.0)
    assert data["percentile"] == pytest_approx(99.0)
    assert data["rank_from_top"] == 86_000_000


def test_cli_requires_one_input_mode(capsys):
    assert cli.main([]) == 2
    assert "provide exactly one" in capsys.readouterr().err


def test_batch_requires_phred_and_preserves_raw(tmp_path: Path):
    source = tmp_path / "in.csv"
    target = tmp_path / "out.csv"
    source.write_text("variant_id,cadd_raw,cadd_phred\nv1,1.2,20\n", encoding="utf-8")
    assert cli.process_batch_csv(str(source), str(target)) == 0
    row = next(csv.DictReader(target.open(encoding="utf-8")))
    assert row["cadd_raw"] == "1.2"
    assert float(row["top_percent"]) == pytest_approx(1.0)
    assert row["acmg_code"] == "NOT_ASSESSED"


def test_batch_missing_phred_fails(tmp_path: Path, capsys):
    source = tmp_path / "in.csv"
    source.write_text("variant_id,cadd_raw\nv1,1.2\n", encoding="utf-8")
    assert cli.process_batch_csv(str(source)) == 1
    assert "missing CADD PHRED" in capsys.readouterr().err
