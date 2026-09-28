# Generalized multiplication tables in every fixed dimension

**Tseng**

[Read the paper](generalized-multiplication-tables.pdf)

This repository contains the complete manuscript, including its technical
appendices, bibliography, vector figure, and code for reproducing the finite
rational-interval checks in the worked examples.

For each fixed integer $k\geq 6$, the paper determines the order of rectangular
$(k+1)$-dimensional multiplication tables for arbitrary ordered side lengths,
and gives a uniform estimate for $H^{(k+1)}$ on the stated dyadic comparison
range. The definitions, hypotheses, and proofs are in the manuscript.

## Contents

| Path | Contents |
| --- | --- |
| `generalized-multiplication-tables.pdf` | Compiled paper |
| `main.tex` | LaTeX entry point |
| `sections/` | Main text |
| `appendices/` | Technical appendices in the same paper |
| `figures/` | TikZ source for the profile diagram |
| `references/references.bib` | All references cited in the paper |
| `main.bbl` | Generated bibliography for submission systems |
| `build.py` | Build helper that keeps temporary files outside the repository |
| `reproducibility/` | Worked-example verification code and its reproducible results |

## Build the paper

Run commands from the repository root. No files from another project or research
directory are required.

### Tectonic

Install Tectonic and Python 3, then run:

```sh
python3 build.py
```

This regenerates `generalized-multiplication-tables.pdf` and `main.bbl`.
Tectonic obtains standard TeX resources on its first run. Once these resources
are cached, `python3 build.py --offline` builds without network access.
If Tectonic is not on `PATH`, pass `--compiler /path/to/tectonic`.

### TeX Live or MiKTeX

A standard LaTeX installation with BibTeX can also compile the sources:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

These commands produce `main.pdf`. The included `main.bbl` is a convenience;
the full bibliography can be regenerated from `references/references.bib`.

The sources use the standard AMS packages, Latin Modern, geometry, microtype,
enumitem, booktabs, longtable, array, graphicx, TikZ and hyperref. No shell escape,
external images, proprietary fonts, custom document class, or custom bibliography
style is needed.

## Reproduce the worked-example checks

With Python 3, run:

```sh
python3 reproducibility/verify_examples.py
```

The script uses only the Python standard library, performs exact rational
interval calculations, prints `PASS`, and regenerates
`reproducibility/examples_results.json`. The checked-in result is provided for
comparison. Run Python normally, without `-O`, because assertions verify the
finite inequalities.

The checks cover logarithm and exponential enclosures, root isolation, KKT
signs, the 15 comparisons in the six-dimensional ray, and the 120 terminal
comparisons in the 18-dimensional family. The analytic derivations and exact
identities are given in Section 4 and Appendix C; the script supplements those
proofs.
