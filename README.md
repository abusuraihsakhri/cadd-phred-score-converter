# CADD Raw <-> PHRED Score Converter & Variant Classifier

A deterministic, high-throughput computational genomics engine for converting **Combined Annotation Dependent Depletion (CADD)** raw scores to **PHRED-scaled ranks**, genome-wide percentiles, and **ACMG/AMP** in-silico pathogenicity evidence tiers.

---

## Genomic & Theoretical Foundation

CADD quantitatively scores the deleteriousness of single nucleotide variants (SNVs) and insertion/deletions across the human reference genome (GRCh37/GRCh38, containing $\approx 8.8 \times 10^9$ possible SNVs).

### Mathematical PHRED Scaling Formula

The PHRED-scaled CADD score quantifies relative deleteriousness on a logarithmic scale:
$$\text{PHRED} = -10 \log_{10}(p) = -10 \log_{10}\left(\frac{\text{Rank}}{N}\right)$$
where:
- $p \in (0, 1]$ is the fraction of all $N \approx 8.8 \times 10^9$ possible SNVs with a raw score greater than or equal to the variant.
- $\text{Percentile} = 100 \times (1 - p) = 100 \times \left(1 - 10^{-\frac{\text{PHRED}}{10}}\right)$
- $\text{Rank} = N - \text{round}(N \times p) + 1$

### PHRED Thresholds & Genomic Rarity

| PHRED Score | Tail Fraction ($p$) | Genomic Percentile | Human Genome Rank |
| :--- | :--- | :--- | :--- |
| **$\text{PHRED } 10$** | $10^{-1} = 0.1$ | $90.0\%$ (Top $10\%$) | Top $8.8 \times 10^8$ SNVs |
| **$\text{PHRED } 20$** | $10^{-2} = 0.01$ | $99.0\%$ (Top $1\%$) | Top $8.8 \times 10^7$ SNVs |
| **$\text{PHRED } 30$** | $10^{-3} = 0.001$ | $99.9\%$ (Top $0.1\%$) | Top $8.8 \times 10^6$ SNVs |
| **$\text{PHRED } 40$** | $10^{-4} = 0.0001$ | $99.99\%$ (Top $0.01\%$) | Top $8.8 \times 10^5$ SNVs |

### ACMG/AMP Clinical Pathogenicity Guidelines

- **$\text{PHRED } < 10.0$**: **BP4 (Benign Supporting)**. In the bottom $90\%$ of mutations across the genome; favored as benign/tolerated.
- **$10.0 \le \text{PHRED} < 15.0$**: **Indeterminate / Neutral**. Intermediate impact; insufficient computational signal for classification.
- **$15.0 \le \text{PHRED} < 20.0$**: **Moderate In-Silico Impact**. Top $3.16\%$ to $1\%$ of variants; candidate for secondary review.
- **$20.0 \le \text{PHRED} < 30.0$**: **PP3 (Pathogenic Supporting)**. In the top $1\%$ to $0.1\%$ most deleterious variants; supports pathogenic assertion.
- **$\text{PHRED } \ge 30.0$**: **PP3 Strong In-Silico**. Top $0.1\%$ in genome; strong computational indication of functional disruption.

---

## Installation & Setup

Requires **Python 3.9+** (standard library only).

```bash
git clone https://github.com/abusuraihsakhri/cadd-phred-score-converter.git
cd cadd-phred-score-converter
```

---

## CLI Usage Examples

### 1. Run Pre-Configured Benchmark Variants

```bash
python cli.py --demo benign_synonymous
python cli.py --demo intermediate_missense
python cli.py --demo pathogenic_missense
python cli.py --demo severe_stopgain
```

### 2. Single Variant Conversion with JSON Output

```bash
python cli.py --variant-id chr17:7577121_C>T --phred 34.0 \
  --consequence stop_gained --phylop 7.2 --gerp 5.1 --json
```

### 3. Batch CSV Conversion

```bash
python cli.py --batch-csv input_variants.csv --output annotated_variants.csv
```

### 4. Interactive Lookup

```bash
python cli.py --interactive
```

---

## Python API Integration

```python
from cadd_phred_converter import (
    VariantInput,
    MolecularConsequence,
    CaddPhredConverter,
    format_cadd_report,
)

variant = VariantInput(
    variant_id="rs121913343",
    phred_score=28.5,
    consequence=MolecularConsequence.MISSENSE,
    phylop_score=5.20,
    gerp_score=4.80
)

result = CaddPhredConverter.evaluate_variant(variant)
print(format_cadd_report(result))
```

---

## Unit Testing

Run the automated test suite with 20 unit test cases:

```bash
python -m unittest test_cadd_phred_converter.py -v
```

---

## License

MIT License. Authored and maintained by Dr. Abu Suraih Sakhri.
