# Early Fusion Multi-Omics Machine Learning Pipeline

**Author:** Susanta Behura  
**Version:** v1.0.0  
**DOI:** (to be assigned via Zenodo)  

## Overview

This repository provides a reproducible pipeline for **integrated analysis of host transcriptomics and gut microbiome data** using:

- Early data fusion
- Random Forest machine learning

The goal is to identify **host genes and microbial taxa associated with susceptibility or resistance to *Eimeria* infection in chicken**.


## Key Features

- Multi-omics integration (gene expression + microbiome)
- Variance-based feature filtering
- CLR-like microbiome transformation
- Early feature fusion
- Random Forest classification
- Feature importance ranking
- PCA visualization

## Input Files

Place the following files in the working directory:

1. `gene_expression.csv`  
   - Rows = genes  
   - Columns = samples  

2. `microbiome_abundance.csv`  
   - Rows = taxa  
   - Columns = samples  

3. `metadata.csv`  
   - Must contain:
     - `SampleID`
     - `Group` (e.g., Resistant / Susceptible)

## How to Run

```bash
Rscript scripts/early_fusion_pipeline.R
