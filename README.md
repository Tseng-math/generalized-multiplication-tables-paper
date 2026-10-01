# Generalized multiplication tables in every fixed dimension

**Tseng**

[Read the paper](generalized-multiplication-tables.pdf)

The complete manuscript, bibliography, figure source, and exact rational-interval
checks for the worked examples are included in this repository.

## Contents

| Path | Contents |
| --- | --- |
| `main.tex` | Manuscript entry point |
| `sections/` | Main text |
| `appendices/` | Technical appendices |
| `figures/` | TikZ figure source |
| `references/references.bib` | Bibliography |
| `generalized-multiplication-tables.pdf` | Final paper |
| `build.py` | Build helper |
| `reproducibility/` | Verification script and reference results |

## Build

Requires Python 3.7 or later and Tectonic on `PATH`. Run from the repository root:

```sh
python3 build.py
```

This rebuilds `generalized-multiplication-tables.pdf`, including the bibliography.
Intermediate files are created in a temporary directory and removed automatically.
All manuscript inputs are included; Tectonic downloads standard TeX resources on
the first build, which requires an internet connection.

Once those resources are cached, build offline with:

```sh
python3 build.py --offline
```

If Tectonic is not on `PATH`, specify its executable:

```sh
python3 build.py --compiler /path/to/tectonic
```

## Reproduce the worked examples

```sh
python3 reproducibility/verify_examples.py
```

The script uses only the Python standard library. It verifies rational bounds,
root isolations, KKT signs, all 15 ray comparisons and all 120 terminal comparisons,
then prints `PASS` and regenerates `reproducibility/examples_results.json`.
The included JSON contains the reference results. Run Python without `-O`, since
the checks use assertions. The corresponding analytic proofs are in the
worked-examples section and appendix of the paper.

## AI Disclosure

The author selected the problem and managed the project; the mathematical derivations were carried out by artificial intelligence.
