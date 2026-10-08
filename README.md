# Dyadic divisor products and rectangular multiplication tables

**Author: Tseng**

[Read the paper](dyadic-divisor-products-and-rectangular-multiplication-tables.pdf).
This repository contains the LaTeX source, compiled manuscript, and
reproducibility checks for the worked examples.

For every fixed integer $k\ge1$, the paper gives a uniform
order-of-magnitude formula for the dyadic divisor-product count

$$
H^{(k+1)}(x,\mathbf y,2\mathbf y),\qquad
3\le y_1\le\cdots\le y_k,\quad
x\ge2^k(y_1\cdots y_k)y_k.
$$

It also gives a formula for the number of distinct products in every
rectangular multiplication table with real ordered sides
$1\le N_1\le\cdots\le N_{k+1}$, including bounded sides and
coincident cutoffs. Comparison constants depend only on $k$.
The manuscript also develops conditional persistence estimates and a
scalar large-deviation principle.

## Contents

| File | Purpose |
| --- | --- |
| [main.tex](main.tex) | Complete manuscript, including appendices and bibliography |
| [dyadic-divisor-products-and-rectangular-multiplication-tables.pdf](dyadic-divisor-products-and-rectangular-multiplication-tables.pdf) | Compiled paper |
| [build.py](build.py) | Portable build helper |
| [reproducibility/verify_examples.py](reproducibility/verify_examples.py) | Exact rational checks for the two worked profiles |
| [CITATION.cff](CITATION.cff) | Citation metadata |

The source uses standard LaTeX packages. It has no separate figures,
source inputs, bibliography database, or generated data dependencies.

## Build

Requires Python 3.9 or later and [Tectonic](https://tectonic-typesetting.github.io/)
on `PATH`. From the repository root, run:

```sh
python3 build.py
```

This rebuilds the named PDF above. Intermediate files are created in a
temporary directory and removed automatically. Tectonic may download
standard TeX resources on the first build; an offline build requires
those resources to be cached:

```sh
python3 build.py --offline
python3 build.py --compiler /path/to/tectonic
```

A local Tectonic resource bundle and cache can be selected with
`--bundle /path/to/bundle --cache-dir /path/to/cache`.
Alternatively, with a TeX distribution and `latexmk` installed, run
`python3 build.py --compiler latexmk`. No BibTeX or Biber step is needed.

## Reproduce the worked examples

```sh
python3 reproducibility/verify_examples.py
```

The checks use only Python's standard library and exact rational interval
bounds with explicit series remainders. They check the scalar constants
and geometric comparisons in Section 3.7; the script states its precise
coverage. These checks do not evaluate the conditional persistence
exponent or certify the full theorems. They print `PASS` and leave no
generated files in the repository.

## Cite

Citation metadata is provided in [CITATION.cff](CITATION.cff).
This is an unpublished manuscript; no DOI is assigned here.

```bibtex
@unpublished{TsengDyadicDivisorProducts2026,
  author = {Tseng},
  title = {Dyadic divisor products and rectangular multiplication tables},
  year = {2026},
  note = {Manuscript},
  url = {https://github.com/Tseng-math/generalized-multiplication-tables-paper}
}
```

## AI disclosure

The author selected the problem and managed the project; the mathematical derivations were carried out by artificial intelligence.
