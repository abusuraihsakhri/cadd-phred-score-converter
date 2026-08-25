#!/usr/bin/env python3
"""
Command-Line Interface for CADD Raw to PHRED Score Converter
============================================================
Supports single variant conversion, interactive lookup, batch CSV/TSV processing,
VCF variant annotation, and ACMG in-silico classification reporting.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import List, Optional

from cadd_phred_converter import (
    TOTAL_GENOME_SNVS,
    MolecularConsequence,
    VariantInput,
    VariantLocation,
    CaddPhredConverter,
    format_cadd_report,
)


def run_demo(scenario: str = "all") -> int:
    """Runs pre-configured clinical benchmark variants."""
    scenarios = {
        "benign_synonymous": VariantInput(
            variant_id="DEMO-VAR-BENIGN",
            location=VariantLocation("7", 117199644, "A", "G", "CFTR"),
            raw_score=0.15,
            consequence=MolecularConsequence.SYNONYMOUS,
            phylop_score=-0.42
        ),
        "intermediate_missense": VariantInput(
            variant_id="DEMO-VAR-INTERMEDIATE",
            location=VariantLocation("13", 32914438, "C", "T", "BRCA2"),
            raw_score=1.45,
            consequence=MolecularConsequence.MISSENSE,
            phylop_score=1.20
        ),
        "pathogenic_missense": VariantInput(
            variant_id="DEMO-VAR-PATHOGENIC",
            location=VariantLocation("17", 43094861, "G", "A", "BRCA1"),
            phred_score=26.4,
            consequence=MolecularConsequence.MISSENSE,
            phylop_score=4.85,
            gerp_score=4.20
        ),
        "severe_stopgain": VariantInput(
            variant_id="DEMO-VAR-STOPGAIN",
            location=VariantLocation("17", 7577121, "C", "T", "TP53"),
            phred_score=36.0,
            consequence=MolecularConsequence.STOP_GAINED,
            phylop_score=7.60,
            gerp_score=5.50
        )
    }

    selected = scenarios.items() if scenario == "all" else [(scenario, scenarios[scenario])] if scenario in scenarios else []
    if not selected:
        print(f"Unknown scenario: {scenario}. Choose from: {list(scenarios.keys())} or 'all'")
        return 1

    for name, var in selected:
        res = CaddPhredConverter.evaluate_variant(var)
        print(format_cadd_report(res))
        print("\n")
    return 0


def interactive_mode() -> int:
    """Guides user through entering variant parameters."""
    print("=" * 60)
    print(" CADD Raw <-> PHRED Score Converter - Interactive Lookup")
    print("=" * 60)
    try:
        var_id = input("Enter Variant ID [VAR-2026-001]: ").strip() or "VAR-2026-001"
        mode = input("Input mode: (1) Raw CADD score, (2) PHRED score [2]: ").strip() or "2"

        raw_val = None
        phred_val = None
        if mode == "1":
            raw_str = input("Enter CADD Raw Score [3.2]: ").strip() or "3.2"
            raw_val = float(raw_str)
        else:
            phred_str = input("Enter CADD PHRED Score [22.5]: ").strip() or "22.5"
            phred_val = float(phred_str)

        csq_str = input("Molecular Consequence [missense_variant]: ").strip() or "missense_variant"
        csq = MolecularConsequence(csq_str) if csq_str in [c.value for c in MolecularConsequence] else MolecularConsequence.MISSENSE

        phylop_str = input("PhyloP score (optional, enter to skip): ").strip()
        phylop = float(phylop_str) if phylop_str else None

        gerp_str = input("GERP++ score (optional, enter to skip): ").strip()
        gerp = float(gerp_str) if gerp_str else None

        variant = VariantInput(
            variant_id=var_id,
            raw_score=raw_val,
            phred_score=phred_val,
            consequence=csq,
            phylop_score=phylop,
            gerp_score=gerp
        )

        res = CaddPhredConverter.evaluate_variant(variant)
        print("\n" + format_cadd_report(res))
        return 0
    except Exception as e:
        print(f"Error in interactive lookup: {e}", file=sys.stderr)
        return 1


def process_batch_csv(input_csv: str, output_csv: Optional[str] = None) -> int:
    """Processes batch CSV file with variants."""
    try:
        with open(input_csv, mode="r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        results = []
        for r in rows:
            vid = r.get("variant_id") or r.get("id") or "VAR-001"
            raw_str = r.get("raw_score") or r.get("cadd_raw")
            phred_str = r.get("phred_score") or r.get("cadd_phred") or r.get("phred")
            raw_val = float(raw_str) if raw_str else None
            phred_val = float(phred_str) if phred_str else None

            if raw_val is None and phred_val is None:
                phred_val = 15.0

            csq_val = r.get("consequence", "missense_variant")
            csq = MolecularConsequence(csq_val) if csq_val in [c.value for c in MolecularConsequence] else MolecularConsequence.MISSENSE

            var = VariantInput(
                variant_id=vid,
                raw_score=raw_val,
                phred_score=phred_val,
                consequence=csq
            )
            rep = CaddPhredConverter.evaluate_variant(var)
            row_dict = dict(r)
            row_dict["cadd_raw"] = rep.raw_score
            row_dict["cadd_phred"] = rep.phred_score
            row_dict["percentile"] = rep.percentile
            row_dict["acmg_code"] = rep.acmg_code
            row_dict["is_pathogenic"] = rep.is_pathogenic_predicted
            results.append(row_dict)

        if output_csv:
            with open(output_csv, mode="w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
                writer.writeheader()
                writer.writerows(results)
            print(f"Successfully converted {len(results)} variants -> {output_csv}")
        else:
            print(json.dumps(results, indent=2))
        return 0
    except Exception as e:
        print(f"Error in batch conversion: {e}", file=sys.stderr)
        return 1


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="CADD Raw <-> PHRED Score Converter & ACMG Variant Classifier"
    )
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive converter mode")
    parser.add_argument("--demo", choices=["benign_synonymous", "intermediate_missense", "pathogenic_missense", "severe_stopgain", "all"], help="Run benchmark demo scenario")
    parser.add_argument("--variant-id", default="VAR-001", help="Variant identifier (e.g. rs12345 or chr17:43094861_G>A)")
    parser.add_argument("--raw", type=float, help="CADD raw score")
    parser.add_argument("--phred", type=float, help="CADD PHRED-scaled score")
    parser.add_argument("--consequence", default="missense_variant", help="Sequence Ontology consequence term")
    parser.add_argument("--phylop", type=float, help="PhyloP conservation score")
    parser.add_argument("--gerp", type=float, help="GERP++ RS score")
    parser.add_argument("--batch-csv", help="Input CSV file for batch conversion")
    parser.add_argument("--output", "-o", help="Output file path (CSV or JSON)")
    parser.add_argument("--file", "-f", help="Load variant JSON file")
    parser.add_argument("--json", "-j", action="store_true", help="Output result as JSON")

    args = parser.parse_args(argv)

    if args.interactive:
        return interactive_mode()

    if args.demo:
        return run_demo(args.demo)

    if args.batch_csv:
        return process_batch_csv(args.batch_csv, args.output)

    if args.file:
        with open(args.file, "r") as fp:
            data = json.load(fp)
        csq_val = data.get("consequence", "missense_variant")
        csq = MolecularConsequence(csq_val) if csq_val in [c.value for c in MolecularConsequence] else MolecularConsequence.MISSENSE
        var = VariantInput(
            variant_id=data.get("variant_id", "FILE-VAR"),
            raw_score=data.get("raw_score"),
            phred_score=data.get("phred_score"),
            consequence=csq,
            phylop_score=data.get("phylop_score"),
            gerp_score=data.get("gerp_score")
        )
    else:
        # Direct arguments
        if args.raw is None and args.phred is None:
            # Default to PHRED 20.0
            phred_val = 20.0
            raw_val = None
        else:
            raw_val = args.raw
            phred_val = args.phred

        csq = MolecularConsequence(args.consequence) if args.consequence in [c.value for c in MolecularConsequence] else MolecularConsequence.MISSENSE
        var = VariantInput(
            variant_id=args.variant_id,
            raw_score=raw_val,
            phred_score=phred_val,
            consequence=csq,
            phylop_score=args.phylop,
            gerp_score=args.gerp
        )

    report = CaddPhredConverter.evaluate_variant(var)

    if args.json:
        out_str = report.to_json()
    else:
        out_str = format_cadd_report(report)

    if args.output:
        with open(args.output, "w") as fp:
            fp.write(out_str)
    else:
        print(out_str)

    return 0


if __name__ == "__main__":
    sys.exit(main())
