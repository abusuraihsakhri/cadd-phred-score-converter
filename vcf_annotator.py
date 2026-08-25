#!/usr/bin/env python3
"""
VCF annotator: add CADD_PHRED INFO fields to a Variant Call Format stream.

Parses VCF line-by-line with stdlib only:
  - preserves header/meta lines verbatim
  - splits multi-allelic ALT records into one output row per allele
  - appends CADD_PHRED=<score> to each record's INFO column
  - optional --filter-threshold keeps only variants at or above a cutoff

Scoring accepts either a precomputed PHRED value or a raw CADD score via
the exponential-tail model in phred_scale.py.
"""

from typing import Callable, Iterable, List, Optional

from phred_scale import raw_to_phred_exponential_tail


def parse_info(info_field: str) -> dict:
    out = {}
    if info_field in (".", ""):
        return out
    for item in info_field.split(";"):
        key, _, value = item.partition("=")
        out[key] = value if _ else True
    return out


def format_info(info_map: dict) -> str:
    parts = []
    for k, v in info_map.items():
        parts.append(k if v is True else f"{k}={v}")
    return ";".join(parts) if parts else "."


def score_variant(chrom: str, pos: int, ref: str, alt: str,
                  scorer: Callable[[str, int, str, str], float]) -> float:
    """Hook for real CADD models; default returns a deterministic surrogate
    derived from the variant string so demos are reproducible."""
    return scorer(chrom, pos, ref, alt)


def surrogate_scorer(chrom: str, pos: int, ref: str, alt: str) -> float:
    """Deterministic pseudo-CADD raw score for offline demonstration."""
    h = (hash((chrom, pos, ref, alt)) & 0xFFFF) / 0xFFFF
    return round(5.0 + 40.0 * h * h, 2)


def annotate_vcf(lines: Iterable[str],
                 threshold: Optional[float] = None,
                 scorer: Callable[[str, int, str, str], float] = surrogate_scorer
                 ) -> List[str]:
    """Return annotated VCF lines; adds ##CADD=<...> meta header."""
    out: List[str] = []
    meta_added = False
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#"):
            if not meta_added and line.startswith("#CHROM"):
                out.append("##INFO=<ID=CADD_PHRED,Number=1,Type=Float,"
                           "Description=\"CADD PHRED-scaled deleteriousness\">")
                meta_added = True
            out.append(line)
            continue
        fields = line.split("\t")
        chrom, pos_s, ref, alts = fields[0], int(fields[1]), fields[3], fields[4]
        info_idx = 7 if len(fields) > 7 else None
        for alt in alts.split(","):
            phred = raw_to_phred_exponential_tail(scorer(chrom, pos_s, ref, alt))
            if threshold is not None and phred < threshold:
                continue
            row = list(fields)
            existing = parse_info(row[info_idx]) if info_idx else {}
            existing["CADD_PHRED"] = f"{phred:.2f}"
            target_idx = info_idx if info_idx is not None else len(row) - 1
            row[target_idx] = format_info(existing)
            out.append("\t".join(row))
    return out


def summarize_annotated(annotated_lines: List[str]) -> dict:
    scores = []
    for line in annotated_lines:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 7:
            info = parse_info(f[7])
            if "CADD_PHRED" in info:
                scores.append(float(info["CADD_PHRED"]))
    scores.sort(reverse=True)
    n = len(scores)
    median = scores[n // 2] if n else 0.0
    return {"variants": n,
            "max_phred": max(scores) if scores else None,
            "median_phred": median,
            "top_quartile_min": scores[n // 4] if n >= 4 else None}


if __name__ == "__main__":
    mini_vcf = "\n".join([
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t1000500\trs123\tA\tG\t.\tPASS\t.",
        "chr2\t20007500\trs456\tC\tT,G\t.\tPASS\t.",
        "chr7\t117120017\t.\tG\tA\t.\tPASS\tAC=1",
    ])
    annotated = annotate_vcf(mini_vcf.splitlines(), threshold=15.0)
    print("Annotated VCF:")
    print("\n".join(annotated))
    print("\nSummary:", summarize_annotated(annotated))
