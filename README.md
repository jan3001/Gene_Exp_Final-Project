# Gene Expression Analysis Pipeline

An end-to-end R script performing comprehensive gene expression analysis on the GSE53890 dataset, including data acquisition, preprocessing, outlier detection, gene filtering, statistical testing, dimensionality reduction, classification, discriminant gene extraction, and functional enrichment.

---

## Table of Contents

* [Overview](#overview)
* [Features](#features)
* [Getting Started](#getting-started)

  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Pipeline Steps](#pipeline-steps)
* [Repository Structure](#repository-structure)
* [Example Outputs](#example-outputs)
* [Contributing](#contributing)
* [Contact](#contact)

---

## Overview

This repository contains a single R script, `analysis_pipeline.R`, that guides you through a full gene expression workflow:

1. **Data Acquisition & Preprocessing**: Download GSE53890, clean and align expression data with sample metadata.
2. **Outlier Detection & Removal**: PCA and hierarchical clustering to detect and remove anomalous samples.
3. **Gene Filtering**: Remove low-expressed genes based on a user-defined threshold.
4. **Feature Selection**: Perform t-tests between age groups and adjust p-values (Benjamini–Hochberg).
5. **Dimensionality Reduction**: PCA and heatmap visualization on statistically significant genes.
6. **Classification**: Train an SVM on PCA features, evaluate performance, and visualize decision boundaries.
7. **Discriminant Gene Extraction**: Identify top genes driving PCA components and calculate log fold-changes.
8. **Functional Enrichment**: Map probes to gene symbols and perform GO enrichment analysis using clusterProfiler.

These steps showcase reproducible, end-to-end bioinformatics analysis in R, following best practices for clarity, documentation, and dependency management.

---

## Features

* **Automated package management**: Installs missing CRAN and Bioconductor dependencies.
* **Reproducible**: Uses a fixed seed for PCA and clustering.
* **Comprehensive**: Covers data cleaning through to biological interpretation.
* **Single-script**: Simple `Rscript analysis_pipeline.R` execution.

---

## Getting Started

### Prerequisites

* R (version 4.0+)
* Internet connection (to download GEO data)

### Installation

1. Clone the repo:

   ```bash
   git clone https://github.com/jan3001/Gene_Exp_Final-Project.git
   cd Gene_Exp_Final-Project
   ```
2. Install dependencies in R:

   ```r
   install.packages(c(
     "ggplot2","e1071","pheatmap",
     "AnnotationDbi","clusterProfiler","ggrepel"
   ))
   if (!requireNamespace("BiocManager", quietly=TRUE))
     install.packages("BiocManager")
   BiocManager::install(c(
     "GEOquery","org.Hs.eg.db","hgu133plus2.db"
   ))
   ```

---

## Usage

Run the full analysis pipeline:

```bash
Rscript analysis_pipeline.R
```

Outputs include console logs, figures (PCA plots, heatmaps, venn diagrams, decision boundaries), and enrichment tables. Customize thresholds or output paths by editing the script parameters.

---

## Pipeline Steps

1. **Data Acquisition**: `getGEO("GSE53890")` and initial exploration.
2. **Class Definition**: Parse age groups (`<60`, `≥60`).
3. **Outlier Detection**: PCA scatter plots and dendrogram-based removal.
4. **Gene Filtering**: Remove genes expressed in <20% of samples.
5. **Feature Selection**: T-tests, p-value adjustment, Venn diagrams of significant genes.
6. **Dimensionality Reduction**: PCA and heatmap on filtered genes.
7. **Classification**: SVM classification on first two PCs with performance metrics and decision boundary plots.
8. **Discriminant Gene Extraction**: Calculate contributions, log fold-changes; plot top genes.
9. **Functional Enrichment**: Map probe IDs to symbols; GO enrichment and dotplots.

---

## Repository Structure

```
Gene_Exp_Final-Project/
├── analysis_pipeline.R  # Main R script
├── README.md           # Project documentation
├── .gitignore          # Ignored files (R history, outputs)
└── outputs/            # Auto-generated figures and tables
```

---

## Example Outputs

<p align="center">
  <img src="https://github.com/jan3001/Gene_Exp_Final-Project/raw/main/outputs/PCA_plot.png" alt="PCA Plot" width="45%">
  <img src="https://github.com/jan3001/Gene_Exp_Final-Project/raw/main/outputs/heatmap.png" alt="Heatmap" width="45%">
</p>

<p align="center">
  <img src="https://github.com/jan3001/Gene_Exp_Final-Project/raw/main/outputs/venn_diagram.png" alt="Venn Diagram" width="30%">
  <img src="https://github.com/jan3001/Gene_Exp_Final-Project/raw/main/outputs/decision_boundary.png" alt="Decision Boundary" width="30%">
</p>

---

## Contributing

Please open issues or pull requests for enhancements, bug fixes, or suggestions.

---

## Contact

**Anjana Suresh**
GitHub: [@jan3001](https://github.com/jan3001)
