# Cadd Phred Score Converter

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

CADD Converter Bridge Interface
===============================
Exports core CADD and PHRED conversion models and functions.

CADD Raw to PHRED Score Converter & Genomic Variant Deleteriousness Engine
==========================================================================
Comprehensive computational genomics module implementing:
- Exact mathematical conversion between raw CADD scores, PHRED-scaled scores,
  genomic percentiles, and 8.8-billion SNV reference ranks.
- Empirical calibration model matching CADD v1.6 / v1.7 genome-wide distribution.
- ACMG/AMP variant interpretation guidelines (PP3 / BP4 in silico evidence codes).
- Multi-evidence integration with PhyloP, GERP++, and molecular consequence types.
- VCF (Variant Call Format) and TSV/CSV batch annotation.

Standards:
- Rentzsch et al. (Nucleic Acids Res 2019 / 2021) CADD v1.6 / v1.7
- Kircher et al. (Nat Genet 2014) "A general framework to estimate the relative pathogenicity of human genetic variants"
- Richards et al. (Genet Med 2015) ACMG/AMP Standards and Guidelines for Sequence Variant Interpretation

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`MolecularConsequence`** — dedicated module for molecular consequence evaluation and state verification.
- **`AcmgEvidenceTier`** — dedicated module for acmg evidence tier evaluation and state verification.
- **`VariantLocation`**: Genomic coordinate for a variant.
- **`VariantInput`**: Input representation of a sequence variant for CADD scoring.
- **`CaddConversionResult`**: Full quantitative breakdown of CADD conversion and clinical pathogenicity tier.
- **`CaddPhredConverter`**: Mathematical and empirical conversion engine between CADD raw scores,
PHRED-scaled ranks, and ACMG interpretation thresholds.

---

## 📐 Mathematical Formulation & Logic

```text
  Calculates PHRED score from tail fraction p:
  Formula: PHRED = -10 * log10(p)
  Calculates tail fraction p from PHRED score:
  Formula: p = 10^(-PHRED / 10)
  Calculates genome-wide percentile from PHRED score:
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --interactive <value> --demo <value> --variant-id <value> --raw <value>
```

### Parameter Reference
- `--interactive`: Specifies input measurement or parameter value.
- `--demo`: Specifies input measurement or parameter value.
- `--variant-id`: Specifies input measurement or parameter value.
- `--raw`: Specifies input measurement or parameter value.
- `--phred`: Specifies input measurement or parameter value.
- `--consequence`: Specifies input measurement or parameter value.
- `--phylop`: Specifies input measurement or parameter value.
- `--gerp`: Specifies input measurement or parameter value.
- `--batch-csv`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t cadd-phred-score-converter .
docker run -p 8000:8000 cadd-phred-score-converter
```
