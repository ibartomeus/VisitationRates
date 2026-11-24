[![DOI](https://zenodo.org/badge/1050989794.svg)](https://doi.org/10.5281/zenodo.24499316)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-blue.svg)](LICENSE-data)

# Predicting optimal pollinator visitation rates to crops using simple mechanistic models

Ignasi Bartomeus^1, Maxime Eeraerts^1, Mike Garratt^2, Jake Bishop^2, Ellen Knight^2, Tom Breeze^2

^1 Estación Biológica de Doñana (EBD-CSIC). Sevilla, Spain.
^2 School of Agriculture, Policy and Development, University of Reading, Reading, Berkshire, RG6 6AR, UK.

## Overview

This repository contains data and code to reproduce the paper: 
Predicting optimal pollinator visitation rates to crops using simple mechanistic models
Ignasi Bartomeus, Maxime Eeraerts, Mike Garratt, Jake Bishop, Ellen Knight, Tom Breeze.

The model functions live in the companion package
[`VisitationR`](https://github.com/ibartomeus/VisitationR).

## Repository structure

```
data/
  VisitationRates_Data.csv          # observational visit-response data (5 studies)
  VisitationRates_ParameterData.csv # literature parameters for ~25 crops
  parameters_real_crops.csv         # fitted vs literature model parameters
  README.md                         # data dictionary (columns, units, sources)
  Watermelon_data/                  # raw watermelon records + extraction script
scripts/
  Real_crops.qmd                    # Quarto doc: initial exploration to fit models.
  Final_analysis.R                  # main analysis producing the figures
```

## Dependencies

- R (>= 4.4)
- [`VisitationR`](https://github.com/ibartomeus/VisitationR), installed with
  `remotes::install_github("ibartomeus/VisitationR")`.

## Data

The three top-level CSVs in `data/` are described field-by-field in
[`data/README.md`](data/README.md), including units, missing-value codes, and
the source study or DOI for each entry. The data are licensed under
[CC BY 4.0](LICENSE-data).

## Licence

- Code (`scripts/`): [MIT](LICENSE)
- Data (`data/`): [CC BY 4.0](LICENSE-data)
