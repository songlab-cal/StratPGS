# Welcome to StratPGS

This repo contains results and scripts for reproducing and repurposing our analyses in our work,

> Alan Aw, Jeremy McRae, Elior Rahmani and Yun Song (2024+) "Highly parameterized polygenic scores tend to overfit to population stratification via random effects"

## Dashboard

<p align="center">
<picture>
  <img src="images/dashboard.png" width="400"/>
</picture>
</p>

As the image above shows, we also present our results in an interactive dashboard, available [here](https://alan-aw.shinyapps.io/stratPGS_v0/).

# Directory Structure

## Results

Under `results`:
- `random_projections` contains results from **Performance Inflation by rPGS**
- `effect_perturbation` contains results from **Performance Relative to pPGS and sPGS**
- `PGS_catalog` contains results from analyses of polygenic risk scores obtained from the [PGS Catalogue](https://www.pgscatalog.org/) (**Evaluation of MCH PGSs** in our paper)
- `PGS_perf_vs_pval_thres` contains the data file and plot showing the performance, averaged across traits, of PGSs we trained on UKB phenotypes using various GWAS p-value thresholds  

## Scripts

Under `scripts`, similar directory organization as Results. We include one additional subdirectory:
- `angular_central_gaussian`, which contains a script to simulate random vectors under the Angular Central Gaussian distribution. This is mentioned briefly in our Main Text and discussed in our Supplementary Material

## Logs

Under `logs`, we provide log files that record the statistical tests we performed, as described in our paper.

# Notes

1. The summary of non-zero variant counts of all PGS Catalogue scores is available under `results/PGS_catalog/all_polygenic_scores`. The associated data wrangling script is available under `scripts/PGS_catalog/all_polygenic_scores`.