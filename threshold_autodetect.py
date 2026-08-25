#!/usr/bin/env python3
"""
Automatic PHRED threshold detection for pathogenic-vs-tolerant separation.

Implements three stdlib-only methods on a cohort of CADD PHRED scores:

  1. Otsu's method: choose the bin cut maximizing between-class variance of
     the histogram (the standard document-thresholding algorithm applied to
     score distributions).
  2. Bimodality coefficient: (skew^2 + 1) / (excess kurtosis + 3*(n-1)^2 /
     ((n-2)(n-3))). Values > 0.555 (log-normal reference) indicate bimodality.
  3. Histogram valley detection: the minimum-density bin between the two
     largest peaks.

Returns an auto-cutoff usable as --filter-threshold in the VCF annotator.
"""

import math
from statistics import mean, median
from typing import Dict, List


def _moments(values: List[float]) -> Dict[str, float]:
    n = len(values)
    m = mean(values)
    var = sum((v - m) ** 2 for v in values) / n
    sd = math.sqrt(var) or 1e-12
    skew = sum(((v - m) / sd) ** 3 for v in values) / n
    kurt = sum(((v - m) / sd) ** 4 for v in values) / n - 3.0   # excess kurtosis
    return {"mean": m, "sd": sd, "skew": skew, "excess_kurtosis": kurt}


def bimodality_coefficient(values: List[float]) -> float:
    """SAS-style BC; >0.555 suggests a bimodal distribution."""
    if len(values) < 4:
        return 0.0
    mo = _moments(values)
    n = len(values)
    denom = mo["excess_kurtosis"] + 3.0 * ((n - 1) ** 2) / ((n - 2) * (n - 3))
    return (mo["skew"] ** 2 + 1.0) / denom if denom else 0.0


def _histogram(values: List[float], bins: int = 64):
    lo, hi = min(values), max(values)
    width = (hi - lo) / bins or 1e-9
    counts = [0] * bins
    for v in values:
        idx = min(bins - 1, int((v - lo) / width))
        counts[idx] += 1
    return lo, hi, width, counts


def otsu_threshold(values: List[float], bins: int = 64) -> float:
    """Otsu's between-class variance maximization over a score histogram."""
    lo, _, width, counts = _histogram(values, bins)
    total = len(values)
    total_sum = sum((lo + (i + 0.5) * width) * c for i, c in enumerate(counts))
    w_b = 0
    sum_b = 0.0
    best_t, best_var = None, -1.0
    for t in range(bins - 1):
        w_b += counts[t]
        if w_b == 0:
            continue
        w_f = total - w_b
        if w_f == 0:
            break
        sum_b += (lo + (t + 0.5) * width) * counts[t]
        mu_b, mu_f = sum_b / w_b, (total_sum - sum_b) / w_f
        between_var = w_b * w_f * (mu_b - mu_f) ** 2
        if between_var > best_var:
            best_var = between_var
            best_t = t
    return round(lo + (best_t + 1) * width, 3)


def valley_threshold(values: List[float], bins: int = 64,
                     smoothing: int = 3) -> Dict[str, float]:
    """Find the density valley between the two tallest histogram peaks."""
    lo, hi, width, counts = _histogram(values, bins)
    smooth = counts[:]
    for i in range(len(counts)):
        window = counts[max(0, i - smoothing):i + smoothing + 1]
        smooth[i] = sum(window) / len(window)
    peak1 = max(range(len(smooth)), key=lambda i: smooth[i])
    peak2 = max((i for i in range(len(smooth)) if abs(i - peak1) > smoothing),
                key=lambda i: smooth[i], default=None)
    if peak2 is None or smooth[peak2] == 0:
        return {"valley_bin_center": otsu_threshold(values),
                "method": "otsu_fallback", "bimodal": False}
    left, right = sorted((peak1, peak2))
    # Search strictly between the peaks so edge artifacts next to a peak
    # cannot masquerade as the density valley.
    lo_search = min(left + smoothing + 1, right)
    hi_search = max(right - smoothing - 1, left)
    if lo_search > hi_search:
        return {"valley_bin_center": otsu_threshold(values),
                "method": "otsu_fallback", "bimodal": False}
    valley = min(range(lo_search, hi_search + 1), key=lambda i: smooth[i])
    center = lo + (valley + 0.5) * width
    return {
        "valley_bin_center": round(center, 3),
        "peak_bins": sorted([left, right]),
        "method": "histogram_valley",
        "bimodal": True,
        "range": [round(lo, 2), round(hi, 2)],
    }


def detect_threshold(scores: List[float]) -> Dict[str, float]:
    bc = bimodality_coefficient(scores)
    result: Dict[str, float] = {
        "otsu_threshold": otsu_threshold(scores),
        "bimodality_coefficient": round(bc, 4),
        "is_bimodal": bc > 0.555,
        "median": median(scores),
    }
    if result["is_bimodal"]:
        v = valley_threshold(scores)
        result["valley_threshold"] = v["valley_bin_center"]
        result["recommended_cutoff"] = v["valley_bin_center"]
    else:
        result["recommended_cutoff"] = result["otsu_threshold"]
    return result


if __name__ == "__main__":
    import random
    rng = random.Random(11)
    tolerant = [rng.gauss(8, 4) for _ in range(600)]
    pathogenic = [rng.gauss(28, 6) for _ in range(250)]
    cohort = [max(1.0, s) for s in tolerant + pathogenic]

    report = detect_threshold(cohort)
    print("CADD PHRED cohort:", len(cohort), "variants")
    print(f"  bimodality coefficient : {report['bimodality_coefficient']} "
          f"(bimodal={report['is_bimodal']})")
    print(f"  Otsu threshold         : {report['otsu_threshold']}")
    print(f"  valley threshold       : {report.get('valley_threshold')}")
    print(f"  recommended cutoff     : {report['recommended_cutoff']}")

    truth = set(id(i) for i in pathogenic)
    above = [s for s in cohort if s >= report['recommended_cutoff']]
    tp = sum(1 for s in pathogenic if s >= report['recommended_cutoff'])
    fp = len(above) - tp
    fn = len(pathogenic) - tp
    precision = tp / (tp + fp) if tp + fp else 0
    recall = tp / (tp + fn)
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
    print(f"  separation quality     : P={precision:.3f} R={recall:.3f} F1={f1:.3f}")
