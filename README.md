# Welcome to StratPRS

This repo contains results and scripts for reproducing and repurposing our analyses in our work,

> Alan Aw, Jeremy McRae, Elior Rahmani and Yun Song (2024+) "Highly parameterized polygenic scores tend to overfit to population stratification via random effects"

We also present our results in an interactive dashboard, available [here](https://alan-aw.shinyapps.io/stratPGS_v0/).

# Directory Structure

## Results

Under `results`:
- `random_projections` contains results from **Performance Inflation by rPGS**
- `effect_perturbation` contains results from **Performance Relative to pPGS and sPGS**
- `PGS_catalog` contains results from analyses of polygenic risk scores obtained from the [PGS Catalogue](https://www.pgscatalog.org/) (**Evaluation of MCH PGSs** in our paper)

## Scripts

Under `scripts`, same directory organization as Results.

# Notes

1. The summary of non-zero variant counts of all PGS Catalogue scores is available under `results/PGS_catalog/all_polygenic_scores`. The associated data wrangling script is available under `scripts/PGS_catalog/all_polygenic_scores`.