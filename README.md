# CADD Raw to PHRED Score Converter & Genomic Variant Deleteriousness Engine

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python: 3.10 | 3.11 | 3.12](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Testing: Pytest](https://img.shields.io/badge/Testing-Pytest-brightgreen.svg?logo=pytest&logoColor=white)](tests/)
[![Genomics: CADD v1.6/v1.7](https://img.shields.io/badge/Genomics-CADD%20v1.6%2Fv1.7-blueviolet.svg)](#formulation--mathematical-principles)

A clinical-grade computational genomics engine for converting Combined Annotation Dependent Depletion (**CADD**) raw scores to **PHRED-scaled** scores, calculating genome-wide rank percentiles across all 8.8 billion possible single nucleotide variants (SNVs) in the human reference genome (GRCh37/GRCh38), and classifying variants according to **ACMG/AMP sequence variant interpretation guidelines** (`PP3` / `BP4` in-silico evidence codes).

---

## 🔬 Formulation & Mathematical Principles

### CADD Modeling Paradigm

Combined Annotation Dependent Depletion (**Kircher et al., Nat Genet 2014**; **Rentzsch et al., Nucleic Acids Res 2019, 2021**) measures variant deleteriousness by contrasting evolutionary derived alleles fixed in the human lineage against simulated de novo mutations. A high-dimensional support vector machine (SVM) / logistic regression classifier integrates >60 functional genomic annotations (conservation scores, epigenetic profiles, sequence context, transcript structure, and protein impacts).

### The PHRED Transformation

Raw CADD scores represent the continuous decision function output from the classifier. Because raw scores have non-linear relative pathogenicity across different annotation regimes, CADD transforms raw scores into a **PHRED-like scaling** based on their genome-wide rank order among all $N = 8.8 \times 10^9$ possible human SNVs:

$$\text{PHRED} = -10 \cdot \log_{10}(p) = -10 \cdot \log_{10}\left(1 - \frac{\text{rank}}{N}\right)$$

Where:
* $N = 8{,}800{,}000{,}000$: Total possible SNVs across the reference genome.
* $\text{rank}$: Variant rank in ascending order of raw score ($1$ to $N$).
* $p = 1 - \frac{\text{rank}}{N}$: The top tail fraction of most deleterious SNVs.

Inverting this formula yields the top tail fraction $p$ and genome-wide percentile:

$$p = 10^{-\frac{\text{PHRED}}{10}}$$

$$\text{Percentile} = 100 \cdot (1 - p) = 100 \cdot \left(1 - 10^{-\frac{\text{PHRED}}{10}}\right)$$

### Interpretation Thresholds & Benchmarks

| PHRED Score | Tail Fraction ($p$) | Genome Percentile | Total SNVs More Deleterious | Clinical Interpretation |
|:-----------:|:-------------------:|:-----------------:|:---------------------------:|:------------------------|
| **10** | Top $10\%$ ($0.10$) | $90.0\%$ | $\sim 880{,}000{,}000$ | Top 10% most deleterious in genome; baseline functional impact |
| **15** | Top $3.16\%$ ($10^{-1.5}$) | $96.84\%$ | $\sim 278{,}278{,}000$ | Potentially deleterious; exceeds screening baseline |
| **20** | Top $1\%$ ($0.01$) | $99.0\%$ | $\sim 88{,}000{,}000$ | Top 1% most deleterious; standard clinical pathogenic threshold |
| **25** | Top $0.316\%$ ($10^{-2.5}$) | $99.68\%$ | $\sim 27{,}828{,}000$ | Top 0.3% most deleterious; strong candidate in rare disease panels |
| **30** | Top $0.1\%$ ($0.001$) | $99.9\%$ | $\sim 8{,}800{,}000$ | Top 0.1% most deleterious; high-confidence pathogenic variant |
| **40** | Top $0.01\%$ ($0.0001$) | $99.99\%$ | $\sim 880{,}000$ | Top 0.01% most deleterious; severe loss-of-function candidate |

---

## 🏛️ ACMG/AMP Sequence Variant Interpretation Guidelines

In accordance with **Richards et al. (Genet Med 2015)** and ClinGen Sequence Variant Interpretation (SVI) recommendations for computational (*in silico*) prediction algorithms:

* **BP4 (Benign Supporting Evidence):** Assigned when $\text{PHRED} < 10.0$ (variant lies within the bottom 90% of genome-wide scores, showing no significant deleterious impact).
* **Indeterminate / Intermediate:** Assigned when $10.0 \le \text{PHRED} < 15.0$; does not meet criteria for either benign or pathogenic supporting evidence.
* **PP3_MODERATE (Potentially Deleterious):** Assigned when $15.0 \le \text{PHRED} < 20.0$; functional impact exceeds benign threshold but remains below the standard clinical PP3 threshold.
* **PP3 (Pathogenic Supporting Evidence):** Assigned when $20.0 \le \text{PHRED} < 30.0$ (variant is in the top 1% to 0.1% most deleterious SNVs in the genome). Multiple in silico algorithms agree on deleterious effect.
* **PP3_STRONG (Pathogenic Strong In-Silico):** Assigned when $\text{PHRED} \ge 30.0$ (variant is in the top 0.1% most deleterious SNVs). Indicates highly disruptive alterations such as premature stop codons or essential splice-site disruptions.

---

## 🚀 Python Quickstart

```python
from cadd_phred_converter import (
    CaddPhredConverter,
    VariantInput,
    VariantLocation,
    MolecularConsequence,
    format_cadd_report,
)

# 1. Direct mathematical transformations
phred = CaddPhredConverter.phred_from_tail_fraction(0.01)   # -> 20.0 (top 1%)
percentile = CaddPhredConverter.percentile_from_phred(20.0) # -> 99.0%
p = CaddPhredConverter.tail_fraction_from_phred(30.0)        # -> 0.001 (top 0.1%)

# 2. Empirical raw <-> PHRED conversion
cadd_phred = CaddPhredConverter.raw_to_phred(4.35)           # -> 24.5
approx_raw = CaddPhredConverter.phred_to_raw(20.0)           # -> 3.0

# 3. Comprehensive variant clinical evaluation
var = VariantInput(
    variant_id="BRCA1:c.5266dupC",
    location=VariantLocation("17", 43094861, "G", "A", gene_symbol="BRCA1"),
    phred_score=26.4,
    consequence=MolecularConsequence.MISSENSE,
    phylop_score=4.85,
    gerp_score=4.20,
)
result = CaddPhredConverter.evaluate_variant(var)

print(f"Variant: {result.variant_id}")
print(f"PHRED Score: {result.phred_score}")
print(f"Percentile: {result.percentile}% (Top {100 - result.percentile:.3f}%)")
print(f"ACMG Code: {result.acmg_code}")
print(f"ACMG Tier: {result.acmg_tier.value}")
print(f"Is Pathogenic: {result.is_pathogenic_predicted}")

# 4. Formatted clinical report
print(format_cadd_report(result))
```

---

## 💻 CLI Usage & Batch Workflow

### 1. Batch CSV Processing

Process an entire batch of genomic variants in CSV format:

```bash
# Process sample.csv and generate annotated CSV
python cli.py batch -i sample.csv -o annotated_variants.csv

# View batch output as JSON in terminal
python cli.py batch -i sample.csv
```

#### Input CSV Format (`sample.csv`)

```csv
chromosome,position,ref_allele,alt_allele,gene_symbol,consequence,cadd_raw,cadd_phred,acmg_recommendation
chr7,117199644,A,G,CFTR,synonymous_variant,-0.42,4.8,Benign Supporting (BP4: PHRED < 10)
chr13,32914438,C,T,BRCA2,missense_variant,1.45,12.25,Intermediate / Indeterminate (PHRED 10-15)
chr17,41244435,G,A,BRCA1,missense_variant,2.85,19.25,Potentially Deleterious (PHRED 15-20)
chr17,43094861,G,A,BRCA1,missense_variant,4.35,24.5,Pathogenic Supporting (PP3: PHRED 20-30)
chr17,7577121,C,T,TP53,stop_gained,8.5,36.25,Pathogenic Strong In-Silico (PP3: PHRED >= 30)
chr2,48033742,G,A,MSH6,frameshift_variant,9.2,38.0,Pathogenic Strong In-Silico (PP3: PHRED >= 30)
```

The batch runner automatically populates:
* `cadd_raw`: CADD raw score (calculated or preserved).
* `cadd_phred`: PHRED-scaled score.
* `percentile`: Genome-wide percentile rank.
* `acmg_code`: `BP4`, `PP3`, `PP3_MODERATE`, or `PP3_STRONG`.
* `acmg_tier`: Full descriptive ACMG evidence tier.
* `is_pathogenic`: Boolean flag indicating predicted pathogenicity.
* `clinical_interpretation`: Human-readable summary for reporting.

### 2. Interactive Mode

```bash
python cli.py --interactive
```

### 3. Single Variant Evaluation

```bash
# Evaluate by PHRED score
python cli.py --variant-id "chr17:43094861_G>A" --phred 26.4 --consequence missense_variant

# Evaluate by Raw CADD score
python cli.py --variant-id "chr17:7577121_C>T" --raw 8.5 --consequence stop_gained --json

# Benchmark demo scenarios
python cli.py --demo all
```

---

## 🧪 Testing & Verification

Run the test suite with pytest:

```bash
python -m pytest -p no:zarr -v
```

Run the end-to-end CLI batch verification:

```bash
python cli.py batch -i sample.csv -o out_smoke.csv
python -c "import csv; assert len(list(csv.DictReader(open('out_smoke.csv')))) == 6; print('Pass!')"
```

---

## 📚 References

1. **Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J.** (2014). *A general framework to estimate the relative pathogenicity of human genetic variants.* Nature Genetics, 46(3), 310–315. [doi:10.1038/ng.2892](https://doi.org/10.1038/ng.2892)
2. **Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M.** (2019). *CADD: predicting the deleteriousness of variants throughout the human genome.* Nucleic Acids Research, 47(D1), D886–D894. [doi:10.1093/nar/gky1016](https://doi.org/10.1093/nar/gky1016)
3. **Richards S, Aziz N, Bale S, et al.** (2015). *Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology.* Genetics in Medicine, 17(5), 405–424. [doi:10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)
