# Gene Expression Outlier Detection & PCA Example

A self‐contained R script that:

1. Downloads GSE53890  
2. Runs PCA for outlier detection  
3. Plots & prints results  

## How to run

\`\`\`bash
Rscript pca_example.R
\`\`\`

## Dependencies

Install in R:

\`\`\`r
install.packages(c(
  "ggplot2","e1071","pheatmap",
  "AnnotationDbi","clusterProfiler","ggrepel"))
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery","org.Hs.eg.db","hgu133plus2.db"))
\`\`\`
