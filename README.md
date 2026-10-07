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
