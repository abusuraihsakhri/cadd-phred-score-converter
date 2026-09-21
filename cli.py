#!/usr/bin/env python3
"""Command-line interface for CADD PHRED rank interpretation."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from typing import List, Optional

from cadd_phred_converter import (
    DEFAULT_REFERENCE_SNVS,
    CaddPhredConverter,
    MolecularConsequence,
    VariantInput,
    format_cadd_report,
)


def _consequence(value: str) -> MolecularConsequence:
    try:
        return MolecularConsequence(value)
    except ValueError:
        return MolecularConsequence.MISSENSE


def process_batch_csv(input_csv: str, output_csv: Optional[str] = None, reference_size: int = DEFAULT_REFERENCE_SNVS) -> int:
    try:
        with open(input_csv, "r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            raise ValueError("Input CSV contains no data rows.")

        output_rows = []
        for index, row in enumerate(rows, start=2):
            phred_text = row.get("phred_score") or row.get("cadd_phred") or row.get("phred")
            if phred_text in (None, ""):
                raise ValueError(f"Row {index}: missing CADD PHRED value.")
            try:
                phred = float(phred_text)
            except ValueError as exc:
                raise ValueError(f"Row {index}: invalid CADD PHRED value {phred_text!r}.") from exc

            raw_text = row.get("raw_score") or row.get("cadd_raw") or row.get("raw")
            raw = float(raw_text) if raw_text not in (None, "") else None
            variant_id = row.get("variant_id") or row.get("id") or f"row-{index - 1}"
            result = CaddPhredConverter.evaluate_variant(
                VariantInput(
                    variant_id=variant_id,
                    raw_score=raw,
                    phred_score=phred,
                    consequence=_consequence(row.get("consequence", "missense_variant")),
                    reference_size=reference_size,
                )
            )
            enriched = dict(row)
            enriched.update(
                tail_fraction=result.top_fraction_p,
                top_percent=result.top_percent,
                percentile=result.percentile,
                rank_from_top=result.rank_from_top,
                reference_size=result.reference_size,
                acmg_code=result.acmg_code,
                is_pathogenic="",
                interpretation=result.clinical_interpretation,
            )
            output_rows.append(enriched)

        if output_csv:
            with open(output_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(output_rows[0].keys()))
                writer.writeheader()
                writer.writerows(output_rows)
        else:
            print(json.dumps(output_rows, indent=2))
        return 0
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


def interactive_mode(reference_size: int = DEFAULT_REFERENCE_SNVS) -> int:
    try:
        variant_id = input("Variant ID [VAR-001]: ").strip() or "VAR-001"
        phred = float(input("CADD PHRED score: ").strip())
        print(format_cadd_report(CaddPhredConverter.evaluate_variant(
            VariantInput(variant_id=variant_id, phred_score=phred, reference_size=reference_size)
        )))
        return 0
    except (ValueError, EOFError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


def run_demo(reference_size: int = DEFAULT_REFERENCE_SNVS) -> int:
    for phred in (5.0, 10.0, 20.0, 30.0):
        result = CaddPhredConverter.evaluate_variant(
            VariantInput(variant_id=f"CADD-{phred:g}", phred_score=phred, reference_size=reference_size)
        )
        print(format_cadd_report(result), "\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Interpret CADD PHRED-scaled rank values")
    parser.add_argument("--reference-size", type=int, default=DEFAULT_REFERENCE_SNVS,
                        help=f"reference variant count for approximate rank output (default: {DEFAULT_REFERENCE_SNVS:,})")
    subparsers = parser.add_subparsers(dest="command")
    batch = subparsers.add_parser("batch", help="Annotate a CSV that already contains CADD PHRED scores")
    batch.add_argument("-i", "--input", required=True)
    batch.add_argument("-o", "--output")
    batch.add_argument("--reference-size", type=int, default=DEFAULT_REFERENCE_SNVS)

    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--variant-id", default="VAR-001")
    parser.add_argument("--phred", type=float)
    parser.add_argument("--percentile", type=float)
    parser.add_argument("--tail-fraction", type=float)
    parser.add_argument("--raw", type=float, help="optional reported raw score; not converted")
    parser.add_argument("--consequence", default="missense_variant")
    parser.add_argument("--json", "-j", action="store_true")
    parser.add_argument("--output", "-o")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "batch":
        return process_batch_csv(args.input, args.output, args.reference_size)
    if args.interactive:
        return interactive_mode(args.reference_size)
    if args.demo:
        return run_demo(args.reference_size)

    modes = [args.phred is not None, args.percentile is not None, args.tail_fraction is not None]
    if sum(modes) != 1:
        parser.print_usage(sys.stderr)
        print("error: provide exactly one of --phred, --percentile, or --tail-fraction", file=sys.stderr)
        return 2

    try:
        if args.phred is not None:
            phred = args.phred
        elif args.percentile is not None:
            phred = CaddPhredConverter.phred_from_percentile(args.percentile)
        else:
            phred = CaddPhredConverter.phred_from_tail_fraction(args.tail_fraction)
        if not math.isfinite(phred):
            raise ValueError("PHRED must be finite")
        result = CaddPhredConverter.evaluate_variant(
            VariantInput(
                variant_id=args.variant_id,
                raw_score=args.raw,
                phred_score=phred,
                consequence=_consequence(args.consequence),
                reference_size=args.reference_size,
            )
        )
        text = result.to_json() if args.json else format_cadd_report(result)
        if args.output:
            with open(args.output, "w", encoding="utf-8") as handle:
                handle.write(text + "\n")
        else:
            print(text)
        return 0
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
