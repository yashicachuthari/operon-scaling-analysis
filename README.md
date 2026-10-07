# Operon Scaling Analysis

Computational analysis of operon scaling patterns across bacterial chromosomes and plasmids.

## Research Overview

This project investigates how operon organization varies across bacterial replicons, particularly chromosomes and plasmids.

The analysis uses publicly available bacterial genome annotations and UniOP to examine relationships between replicon characteristics and operon organization.

## Research

Research conducted by Yashica Chuthari under the mentorship of Professor Rohan Maddamsetti at Rutgers University.

## Repository Structure

- `src/` — Python scripts for genome retrieval, filtering, and operon analysis
- `data/` — Input genome metadata (not included in repository)
- `results/` — Generated datasets and analysis outputs (not included)
- `UniOP/` — External operon analysis software (not included)

## Main Scripts

- `filter-genome-reports.py` — Filters genome metadata
- `fetch-gbk-annotation.py` — Retrieves genome annotations
- `run-uniop-full.py` — Runs the operon analysis pipeline

## Reproducibility

Detailed installation instructions, software versions, data sources, and commands for reproducing the figures will be added as the pipeline is documented and validated.

## Status

Ongoing research. Results are being independently validated.

## Code Attribution

The genome retrieval and filtering scripts were originally
developed by Professor Rohan Maddamsetti and obtained from:

https://github.com/rohanmaddamsetti/plasmid-scaling-laws

These scripts are used and, where applicable, modified
as part of this research project.

The original repository includes the GNU General Public
License version 2. The relevant license and copyright
notices are retained.

The UniOP software is developed separately:
https://github.com/hongsua/UniOP

The operon scaling analysis and research-specific pipeline
are being developed by Yashica Chuthari under the
mentorship of Professor Rohan Maddamsetti.

## Reproducing the R Figure

### Requirements

- R 4.6.1
- tidyverse (including ggplot2)

### Install R packages

Run:

R -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'

### Generate the figure

From the repository's root directory, run:

Rscript src/plot-operon-scaling.R

The script reads:

processed_data/operon_scaling_results.csv

It generates:

- results/figures/operon_scaling_R.png
- results/figures/operon_scaling_R.pdf

The processed dataset contains 60,934 bacterial replicons,
including 17,444 chromosomes and 43,490 plasmids.

Missing operon counts are treated as zero for plotting.
Replicons with missing or nonpositive lengths are excluded.

### Raw Data and Analysis

Genome annotations were obtained from NCBI RefSeq.
Operon predictions were generated using UniOP.

The complete genome-processing workflow and independent
validation procedures are still being documented.
