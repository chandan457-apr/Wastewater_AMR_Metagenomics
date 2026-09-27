METAGENOMIC PROFILING OF AMR IN WASTEWATER: TAXONOMIC DIVERSITY, ARGS AND MOBILE GENETIC ELEMENTS ACROSS FOUR INDIAN CITIES


## Content

1. Introduction
2. Code
3. Data

## Introduction

This project presents a computational analysis of publicly available shotgun
metagenomic wastewater samples from four Indian metropolitan cities: Chennai,
Delhi, Kolkata and Mumbai.

A selected set of 16 wastewater metagenomes was analyzed to characterize
microbial community composition and antimicrobial resistance using both
read based and genome resolved metagenomic approaches.

Read based analyses included species level taxonomic profiling using Kraken2
and Bracken, alpha diversity, Aitchison based beta diversity analysis, and
antimicrobial resistance profiling using RGI and CARD.

Genome resolved analysis included metagenomic assembly, binning and quality
assessment of metagenome assembled genomes (MAGs), followed by taxonomic
classification and analysis of antimicrobial resistance genes and mobile
genetic elements.

## Code

The R scripts in the `R` directory are used for processing, analysis and
visualization of the project results.

The scripts include normalization of species level abundance data, calculation
of alpha diversity, Aitchison PCoA of taxonomic and ARG profiles, analysis of
top bacterial species and ARG associated drug classes, and MAG level
ARG–MGE analysis.

The required R packages are mentioned in 00_setup.R.

## Data

All processed data required for the R analysis are deposited at
[Zenodo]https://doi.org/10.5281/zenodo.22943902.

The Zenodo archive contains the processed tabular data used for taxonomic,
antimicrobial resistance, MAG and mobile genetic element analyses.

The data should be downloaded and placed in the appropriate `data` directory
before running the R scripts.

The raw sequencing data were obtained from the publicly available source
dataset used for this computational analysis.