# Welcome to StratPRS

This repo contains results and scripts for reproducing and repurposing our analyses in our work,

> Alan Aw, Jeremy McRae, Elior Rahmani and Yun Song (2023+) "Population stratification explains low sensitivities of overly parameterized polygenic risk scores"

# Directory Structure

## Results

Under `results`:
- `random_projections` contains results from **Predictability Inflation by rPRS**
- `effect_perturbation` contains results from **Sensitivity to pPRS and sPRS**
- `PGS_catalog` contains results from analyses of polygenic risk scores obtained from the [PGS Catalogue](https://www.pgscatalog.org/)

## Scripts

Under `scripts`, same directory organization as Results.

# Notes

1. The summary of non-zero variant counts of all PGS Catalogue scores is available under `results/PGS_catalog/all_polygenic_scores`. The associated data wrangling script is available under `scripts/PGS_catalog/all_polygenic_scores`.    