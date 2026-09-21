# CADD PHRED Rank Converter

### [Open the Live Application →](https://abusuraihsakhri.github.io/cadd-phred-score-converter/)

A small Python and browser utility for interpreting CADD PHRED-scaled scores as tail fractions, percentiles, and approximate genome-wide ranks.

## What it does

- Converts CADD PHRED values to top-tail fraction and percentile.
- Estimates a rank/count from the most-deleterious end for a configurable reference size.
- Processes CSV files that already contain CADD PHRED values.
- Provides a static browser interface suitable for GitHub Pages.
- Keeps an optional reported raw CADD score as metadata without fabricating a raw-to-PHRED conversion.

## Important scope

CADD raw scores are model outputs. Their mapping to PHRED-scaled scores is obtained by ranking scores against a release-specific reference distribution; there is no universal raw-score formula implemented here. CADD PHRED values describe relative predicted deleteriousness and are not, by themselves, a pathogenic/benign classification.

This project therefore does **not** automatically assign ACMG/AMP PP3/BP4 evidence. ClinGen guidance requires calibrated predictor use and may be gene- or disease-specific.

## Browser use

Open the live application above or `index.html` directly. Enter a CADD PHRED score and select **Analyze**. The page runs locally in the browser and sends no entered values to a server.

The UI is responsive, defaults to a light theme, and includes a dark-mode toggle.

## Python use

No runtime dependencies are required.

```bash
python cli.py --phred 20
python cli.py --phred 25.3 --json
python cli.py --percentile 99
python cli.py --tail-fraction 0.001
```

Batch CSV input must contain one of `cadd_phred`, `phred_score`, or `phred`:

```bash
python cli.py batch -i sample.csv -o annotated.csv
```

An optional `cadd_raw`/`raw_score` column is preserved as reported metadata; it is not converted.

## Mathematical definition

For a PHRED-scaled CADD value `q`:

```text
tail_fraction = 10^(-q / 10)
top_percent   = 100 * tail_fraction
percentile    = 100 * (1 - tail_fraction)
```

Using a reference size `N`, the approximate rank from the most-deleterious end is:

```text
rank_from_top = ceil(N * tail_fraction)
```

The default `N = 8.6 billion` follows the approximate GRCh37 substitution count described in CADD documentation and can be changed with `--reference-size`.

## Development

```bash
python -m pip install pytest
python -m compileall -q .
python -m pytest -q
node tests/web_smoke.js
```

GitHub Actions runs Python tests on Python 3.10–3.13 and a Node smoke test for the browser calculation code.

## Technology

- Python standard library
- HTML, CSS, and vanilla JavaScript
- GitHub Actions and GitHub Pages

No Pyodide runtime is used: the browser calculations are small deterministic equations, so native JavaScript avoids loading a large WebAssembly/Python runtime while the Python CLI remains available for local and batch workflows.

## Browser compatibility

Current Chrome/Chromium, Firefox, Safari, and Edge versions with standard ES2020-era JavaScript support are expected to work. The layout is responsive for desktop and mobile screens.

## References

- Kircher M, et al. *Nature Genetics*. 2014;46:310–315. doi:10.1038/ng.2892.
- Rentzsch P, et al. *Nucleic Acids Research*. 2019;47:D886–D894. doi:10.1093/nar/gky1016.
- Pejaver V, et al. *American Journal of Human Genetics*. 2022. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.
- CADD documentation: scaled scores are derived from genome-wide rank and raw scores should not be treated as having a universal absolute interpretation.

## License

MIT. See [LICENSE](LICENSE).
