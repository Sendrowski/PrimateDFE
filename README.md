# PrimateDFE <img align="right" width="100" src="resources/logo/large.png">

This repository contains the analysis pipeline used in the study *Comparison of the Distribution of Fitness Effects Across Primates*. The project investigates how DFEs vary across **38 catarrhine primates** and to what extent interspecific differences are explained by variation in effective population size Nₑ.

## Overview

The workflow performs:

- **SFS construction** – VCF parsing, degeneracy annotation (0-fold / 4-fold), and ancestral allele inference using [`fastDFE`](https://github.com/Sendrowski/fastDFE)  
- **DFE inference** – estimation of DFEs with gamma–exponential and discrete models using [`fastDFE`](https://github.com/Sendrowski/fastDFE)  
- **Comparative analyses** – associations between DFE properties and Nₑ, including phylogenetic regressions  

## Workflow

The analysis is implemented as a **Snakemake workflow**.

Structure:
- [`workflow/`](workflow/)
  - [`Snakefile`](workflow/Snakefile) – main workflow entry point
  - [`scripts/`](workflow/scripts/) – Python scripts used by rules
  - [`envs/`](workflow/envs/) – conda environment definitions
- [`resources/`](resources/) – reference data
- [`results/`](results/) – workflow outputs