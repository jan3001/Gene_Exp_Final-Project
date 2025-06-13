# Gene Expression Data Analysis and Visualization
# Final Project

# Load necessary libraries
library(GEOquery)      # For accessing GEO datasets
library(ggplot2)       # For plotting
library(e1071)         # For SVM classification
library(pheatmap)      # For heatmaps
library(AnnotationDbi) # For mapping gene IDs
library(org.Hs.eg.db)  # Human genome annotation
library(clusterProfiler)# For functional enrichment analysis
library(ggrepel)

# Ensure required packages are installed
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  BiocManager::install("hgu133plus2.db")
}
library(hgu133plus2.db)

# Set seed for reproducibility
set.seed(123)

# 1. Data Acquisition and Preprocessing

## 1.1 Load the GEO dataset
gse <- getGEO("GSE53890", GSEMatrix = TRUE)
expr_data <- exprs(gse[[1]])      # Expression data
metadata <- pData(gse[[1]])       # Sample metadata

## 1.2 Explore the data
dim(expr_data)        # Dimensions of expression data (genes x samples)
head(expr_data[, 1:5])# Preview first 5 samples

dim(metadata)         # Dimensions of metadata
colnames(metadata)    # Check available metadata columns

# 2. Define Class Structure (e.g., Age Groups)

## 2.1 Extract age information from metadata
metadata$age <- as.numeric(gsub(" years old", "", metadata$`age:ch1`))

## 2.2 Define age groups
metadata$new_age_group <- ifelse(metadata$age < 60, "<60", "≥60")

## 2.3 Check the distribution of samples in each age group
table(metadata$new_age_group)

age_group_distribution <- table(metadata$new_age_group)

# Bar plot
barplot(age_group_distribution, 
        main = "Sample Distribution by Age Group", 
        xlab = "Age Groups", 
        ylab = "Number of Samples", 
        col = c("skyblue", "lightcoral"), 
        border = "black")


# 3. Align Expression Data with Metadata

## Ensure that the columns of expression data match the rows of metadata
expr_data <- expr_data[, rownames(metadata)]

## Verify alignment
if (!all(colnames(expr_data) == rownames(metadata))) {
  stop("Mismatch between expression data columns and metadata rows.")
}

# 4. Outlier Detection and Removal

## 4.1 Perform PCA for outlier detection
pca <- prcomp(t(expr_data), scale. = TRUE)

## Create a data frame for PCA results
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  AgeGroup = metadata$new_age_group
)

## 4.2 Plot PCA to visualize potential outliers
ggplot(pca_df, aes(x = PC1, y = PC2, color = AgeGroup)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of Samples (New Age Groups)",
    x = "PC1",
    y = "PC2"
  )

pca_df$SampleName <- rownames(pca_df)
ggplot(pca_df, aes(x = PC1, y = PC2, color = AgeGroup)) +
  geom_point(size = 3) +
  geom_text(aes(label = SampleName), vjust = -0.5, size = 3) +  # Adjust vjust and size as needed
  theme_minimal() +
  labs(
    title = "PCA of Samples (New Age Groups)",
    x = "PC1",
    y = "PC2"
  )

## 4.3 Perform Hierarchical Clustering
# Calculate distance matrix
dist_matrix <- dist(t(expr_data))

# Hierarchical clustering
hc_samples <- hclust(dist_matrix, method = "average")

# Label samples with ID and age group
hc_samples$labels <- paste(rownames(metadata), metadata$new_age_group, sep = " - ")

# Plot dendrogram
plot(hc_samples, main = "Hierarchical Clustering of Samples")

# 4.4 Identify and Remove Outliers
# Cut tree into clusters (adjust 'k' as needed)
k <- 2  # Number of clusters
rect.hclust(hc_samples, k = k, border = "red")

# Assign samples to clusters
cluster_assignments <- cutree(hc_samples, k = k)

# Compare clusters with age groups
table(cluster_assignments, metadata$new_age_group)



# Identify samples in the cluster considered as outliers (e.g., cluster 2)
outlier_cluster <- 2  # Adjust based on dendrogram
outlier_samples <- rownames(metadata[cluster_assignments == outlier_cluster, ])

# Print outlier samples
print("Outlier Samples:")
print(outlier_samples)

# Remove outliers from expression data and metadata
expr_data <- expr_data[, !colnames(expr_data) %in% outlier_samples]
metadata <- metadata[!rownames(metadata) %in% outlier_samples, ]

## Verify dimensions after removal
dim(expr_data)
dim(metadata)

# Re-run PCA after removing outliers
pca_clean <- prcomp(t(expr_data), scale. = TRUE)
pca_clean_df <- data.frame(pca_clean$x, AgeGroup = metadata$new_age_group)

# PCA plot after outlier removal
library(ggplot2)
ggplot(pca_clean_df, aes(x = PC1, y = PC2, color = AgeGroup)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA After Outlier Removal",
    x = "PC1",
    y = "PC2"
  )

# Recalculate distance matrix
dist_matrix_clean <- dist(t(expr_data))

# Perform hierarchical clustering
hc_clean <- hclust(dist_matrix_clean, method = "average")

# Dendrogram after outlier removal
plot(hc_clean, main = "Hierarchical Clustering After Outlier Removal", sub = "", xlab = "")

# Optionally add cluster boxes
k <- 2  # Adjust as needed
rect.hclust(hc_clean, k = k, border = "red")

# 5. Gene Filtering

## 5.1 Filter out genes with low expression
# Set a threshold for minimum expression (adjust as needed)
expression_threshold <- 5

# Calculate the number of samples where each gene's expression exceeds the threshold
expressed_samples <- rowSums(expr_data > expression_threshold)

# Set a cutoff 
sample_cutoff <- 0.2 * ncol(expr_data)

# Filter genes
expr_data_filtered <- expr_data[expressed_samples >= sample_cutoff, ]

# Number of genes after filtering
cat("Number of genes after filtering:", nrow(expr_data_filtered), "\n")

## 5.2 Plot histogram of mean expression per gene
gene_means <- rowMeans(expr_data_filtered)
hist(
  gene_means,
  breaks = 50,
  main = "Gene Expression Distribution After Filtering",
  xlab = "Mean Expression"
)

# 6. Feature Selection

## 6.1 Create a factor for age groups
age_group <- factor(metadata$new_age_group)

## 6.2 Perform statistical tests (t-test for two groups)
p_values <- apply(
  expr_data_filtered, 1,
  function(x) t.test(x ~ age_group)$p.value
)

## 6.3 Adjust p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")

## 6.4 Select significant genes (adjusted p-value < 0.05)
significant_genes <- names(adjusted_p_values[adjusted_p_values < 0.05])

# Number of significant genes identified
cat("Number of significant genes:", length(significant_genes), "\n")

## 6.5 Plot p-value distributions
# Raw p-values
hist(
  p_values,
  breaks = 50,
  main = "P-value Distribution",
  xlab = "P-value"
)

# Adjusted p-values
hist(
  adjusted_p_values,
  breaks = 50,
  main = "Adjusted P-value Distribution",
  xlab = "Adjusted P-value"
)

library(VennDiagram)

# Raw significant genes (p-value < 0.05)
raw_significant_genes <- names(p_values[p_values < 0.05])

# Adjusted significant genes (adjusted p-value < 0.05)
adjusted_significant_genes <- names(adjusted_p_values[adjusted_p_values < 0.05])

# Create Venn Diagram
venn.plot <- venn.diagram(
  x = list(
    Raw = raw_significant_genes,
    Adjusted = adjusted_significant_genes
  ),
  filename = NULL,
  main = "Overlap of Significant Genes",
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

# Plot Venn Diagram
dev.off()
grid.draw(venn.plot)



# 7. Dimensionality Reduction and Visualization

## 7.1 Subset expression data to significant genes
expr_data_significant <- expr_data_filtered[significant_genes, ]

## 7.2 Perform PCA on significant genes
pca_significant <- prcomp(t(expr_data_significant), scale. = TRUE)

## Create PCA dataframe for plotting
pca_df_significant <- data.frame(
  PC1 = pca_significant$x[, 1],
  PC2 = pca_significant$x[, 2],
  AgeGroup = metadata$new_age_group
)

## 7.3 Plot PCA of significant genes
ggplot(pca_df_significant, aes(x = PC1, y = PC2, color = AgeGroup)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of Significant Genes",
    x = "PC1",
    y = "PC2"
  )

## 7.4 Heatmap of significant genes
# Normalize data for heatmap
expr_data_norm <- t(scale(t(expr_data_significant)))

# Create annotation for samples
annotation_col <- data.frame(AgeGroup = metadata$new_age_group)
rownames(annotation_col) <- colnames(expr_data_norm)

# Plot heatmap
pheatmap(
  expr_data_norm,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = "Heatmap of Significant Genes"
)

# Add age group information to sample labels
labels_with_age <- paste(rownames(metadata), metadata$new_age_group, sep = " - ")

# Perform hierarchical clustering
hc <- hclust(dist(t(expr_data_significant)), method = "average")

# Update labels in dendrogram
plot(hc, labels = labels_with_age, main = "Dendrogram of Significant Genes with Age Groups")


# 8. Classification

## 8.1 Split data into training and test sets
train_indices <- sample(seq_len(nrow(metadata)), size = 0.7 * nrow(metadata))
test_indices <- setdiff(seq_len(nrow(metadata)), train_indices)

# Sizes of training and test sets
cat("Number of training samples:", length(train_indices), "\n")
cat("Number of test samples:", length(test_indices), "\n")

## 8.2 Prepare training and test data using PCA features
# Use PC1 and PC2 as features
train_pca <- pca_df_significant[train_indices, c("PC1", "PC2")]
test_pca <- pca_df_significant[test_indices, c("PC1", "PC2")]

# Labels for training and testing
train_labels <- age_group[train_indices]
test_labels <- age_group[test_indices]

## 8.3 Train SVM classifier
svm_pca_model <- svm(train_pca, as.factor(train_labels), kernel = "linear")

## 8.4 Predict on test data
pca_predictions <- predict(svm_pca_model, test_pca)

## 8.5 Evaluate classifier performance
# Confusion matrix
pca_conf_matrix <- table(Predicted = pca_predictions, Actual = test_labels)
print("Confusion Matrix:")
print(pca_conf_matrix)

# Calculate accuracy
pca_accuracy <- sum(diag(pca_conf_matrix)) / sum(pca_conf_matrix)
cat("Classification Accuracy (PCA features):", pca_accuracy, "\n")

## 8.6 Identify Important Genes Contributing to PC1 and PC2

# Get the loadings (rotation matrix) from PCA
loadings <- pca_significant$rotation

# Calculate the contribution of each gene to PC1 and PC2
# We can consider the absolute value of loadings as contribution
abs_loadings_PC1 <- abs(loadings[, 'PC1'])
abs_loadings_PC2 <- abs(loadings[, 'PC2'])

# Sum the contributions from PC1 and PC2 for each gene
total_contribution <- abs_loadings_PC1 + abs_loadings_PC2

# Create a data frame with gene names and their total contribution
gene_contributions <- data.frame(
  Gene = rownames(loadings),
  Contribution = total_contribution
)

# Order genes by their total contribution in decreasing order
gene_contributions_ordered <- gene_contributions[order(gene_contributions$Contribution, decreasing = TRUE), ]

# Select top 10 genes with the highest contributions
top_genes_pca <- head(gene_contributions_ordered, 10)

# Print top contributing genes
cat("Top 10 Genes Contributing to PC1 and PC2:\n")
print(top_genes_pca)

## 8.6 Plot Decision Boundary
# Create a grid of points to predict and plot decision boundary
x_min <- min(pca_df_significant$PC1) - 1
x_max <- max(pca_df_significant$PC1) + 1
y_min <- min(pca_df_significant$PC2) - 1
y_max <- max(pca_df_significant$PC2) + 1

# Create a grid of values for PC1 and PC2
grid_points <- expand.grid(
  PC1 = seq(x_min, x_max, length.out = 300),
  PC2 = seq(y_min, y_max, length.out = 300)
)

# Predict on the grid points
grid_predictions <- predict(svm_pca_model, grid_points)

# Combine grid points with predictions for plotting
grid_df <- cbind(grid_points, Predicted = as.factor(grid_predictions))

# Scatter plot with decision boundary
ggplot() +
  geom_tile(data = grid_df, aes(x = PC1, y = PC2, fill = Predicted), alpha = 0.3) +  # Decision boundary
  geom_point(data = pca_df_significant, aes(x = PC1, y = PC2, color = AgeGroup), size = 3) +  # Actual data points
  geom_text_repel(data = pca_df_significant, aes(x = PC1, y = PC2, label = rownames(pca_df_significant)), size = 3) +  # Add sample names
  scale_fill_manual(values = c("<60" = "lightblue", "≥60" = "pink"), name = "Predicted Class") +  # Custom fill colors
  scale_color_manual(values = c("<60" = "blue", "≥60" = "red"), name = "Actual Class") +  # Custom point colors
  labs(
    title = "SVM Decision Boundary (PCA Features)",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

## 8.7 Visualize classification results
# Add predicted and actual classes to test PCA data
test_pca$Predicted <- pca_predictions
test_pca$Actual <- test_labels

# Plot predicted vs. actual classes
ggplot(test_pca, aes(x = PC1, y = PC2, color = Predicted, shape = Actual)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA-Based Classification",
    x = "PC1",
    y = "PC2"
  )

# 9. Extract Top Discriminant Genes

## 9.1 Calculate mean expression for each gene in each group

# Ensure that 'age_group' is correctly assigned
age_group <- factor(metadata$new_age_group)

# Transpose the filtered expression data for aggregation
expr_data_transposed <- t(expr_data_filtered)

# Combine the expression data with the age group
expr_data_with_group <- data.frame(age_group, expr_data_transposed, check.names = FALSE)

# Calculate mean expression for each gene in each group
group_means <- aggregate(. ~ age_group, data = expr_data_with_group, FUN = mean)

# Transpose group_means to get genes as rows
group_means_t <- t(group_means[, -1])  # Exclude 'age_group' column

# Set column names as age groups
colnames(group_means_t) <- group_means$age_group

# Get gene names
gene_names <- rownames(group_means_t)

## 9.2 Calculate log fold changes

# Extract mean expressions for each group
mean_group1 <- group_means_t[, "<60"]
mean_group2 <- group_means_t[, "≥60"]

# Avoid division by zero by adding a small value
mean_group1[mean_group1 == 0] <- 1e-8
mean_group2[mean_group2 == 0] <- 1e-8

# Calculate log2 fold change
logFC <- log2(mean_group2 / mean_group1)

# Assign gene names to logFC
names(logFC) <- gene_names

# Handle infinite or NaN values
logFC[!is.finite(logFC)] <- NA

# Remove NA values
logFC <- na.omit(logFC)

## 9.3 Combine adjusted_p_values and logFC into results data frame

# Find common genes between adjusted_p_values and logFC
common_genes <- intersect(names(adjusted_p_values), names(logFC))

# Subset adjusted_p_values and logFC to common genes
adjusted_p_values_subset <- adjusted_p_values[common_genes]
logFC_subset <- logFC[common_genes]

# Create results data frame
results <- data.frame(
  Gene = common_genes,
  AdjustedPValue = adjusted_p_values_subset,
  logFC = logFC_subset,
  stringsAsFactors = FALSE
)

# Remove any rows with NA values
results <- na.omit(results)

## 9.4 Select Top 5 Positive and Negative Genes

# Ensure logFC is numeric
results$logFC <- as.numeric(results$logFC)

# Order results by logFC in decreasing order for positive logFC
results_ordered_pos <- results[order(results$logFC, decreasing = TRUE), ]

# Order results by logFC in increasing order for negative logFC
results_ordered_neg <- results[order(results$logFC, decreasing = FALSE), ]

# Select top 5 genes with the largest positive logFC
top_positive_genes <- head(results_ordered_pos, 5)

# Select top 5 genes with the largest negative logFC
top_negative_genes <- head(results_ordered_neg, 5)

# Combine top positive and negative genes
top_genes <- rbind(top_positive_genes, top_negative_genes)

# Print top genes
cat("Top 5 Genes with Positive LogFC:\n")
print(top_positive_genes)

cat("\nTop 5 Genes with Negative LogFC:\n")
print(top_negative_genes)

# Combine positive and negative top genes for plotting
top_genes$Direction <- ifelse(top_genes$logFC > 0, "Positive", "Negative")

# Create a bar plot
ggplot(top_genes, aes(x = reorder(Gene, logFC), y = logFC, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  coord_flip() +  # Flip the coordinates for better readability
  labs(
    title = "Top 5 Genes with Positive and Negative logFC",
    x = "Gene",
    y = "log2 Fold Change"
  ) +
  scale_fill_manual(values = c("Positive" = "skyblue", "Negative" = "salmon")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# 10. Functional Annotation and Enrichment Analysis

## 10.1 Map Affymetrix IDs to Gene Symbols
# Using hgu133plus2.db for mapping
library(AnnotationDbi)
library(hgu133plus2.db)

# Map probe IDs to gene symbols
mapped_genes <- mapIds(
  hgu133plus2.db,
  keys = top_genes$Gene,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Remove unmapped IDs
mapped_genes <- na.omit(mapped_genes)

# Print mapped gene symbols
cat("\nMapped Gene Symbols:\n")
print(mapped_genes)

## 10.2 Perform Functional Enrichment Analysis
# Using clusterProfiler for GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

go_results <- enrichGO(
  gene = mapped_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

## 10.3 View and Visualize GO enrichment results
if (!is.null(go_results) && nrow(go_results) > 0) {
  cat("\nGO Enrichment Results:\n")
  print(head(go_results))
  
  # Dotplot of top GO terms
  dotplot(go_results, showCategory = 10) +
    ggtitle("GO Enrichment of Top Genes")
} else {
  cat("\nNo significant GO terms found.\n")
}


# 11. Summary

# - Tested for outlier samples and provided visual proof (PCA plot and dendrogram).
# - Removed identified outliers from the dataset.
# - Filtered out genes with low expression based on defined criteria.
# - Conducted feature selection using appropriate statistical tests (t-test for two groups).
# - Adjusted for multiplicity using Benjamini-Hochberg method.
# - Provided the number of genes retained and thresholds used.
# - Plotted scores (p-values and adjusted p-values) in histograms.
# - Used PCA for dimensionality reduction and visualized samples in two-dimensional space.
# - Used a classification method (SVM) to classify samples into their respective classes.
# - Colored samples by predicted class membership and used different symbols for actual class memberships.
# - Extracted top 5 discriminant genes in positive and negative directions.
# - Provided gene names and functional information (associated pathways, GO terms) for these genes.
#
#
#
# End of Script
