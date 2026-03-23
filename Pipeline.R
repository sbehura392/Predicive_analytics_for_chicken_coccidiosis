############################################################
# Title: Early Fusion Multi-Omics Machine Learning Pipeline
# Author: Susanta Behura
# Description:
#   This pipeline performs integrated analysis of gut microbiome
#   and host gene expression data using early data fusion and
#   Random Forest machine learning to identify key features
#   associated with susceptibility or resistance to Eimeria infection
#   in chicken.
#
# Repository: [Add GitHub URL here]
# DOI: [Add Zenodo DOI after release]
#
# ----------------------------------------------------------
# INPUT FILES:
#   1. gene_expression.csv        (rows = genes, columns = samples)
#   2. microbiome_abundance.csv  (rows = taxa, columns = samples)
#   3. metadata.csv              (SampleID, Group)
#
# OUTPUT FILES:
#   - fused_feature_matrix.csv
#   - train_test_predictions.csv
#   - feature_importance.csv
#   - top20_feature_importance.png
#   - pca_early_fusion.png
#
# REQUIREMENTS:
#   R >= 4.0
#
# KEY METHODS:
#   - Log transformation
#   - Variance filtering
#   - CLR-like microbiome transformation
#   - Early data fusion
#   - Random Forest classification
#
# REPRODUCIBILITY:
#   Random seed fixed (set.seed = 123)
############################################################


# ==========================================================
# 1. LOAD REQUIRED PACKAGES
# ==========================================================
required_packages <- c("data.table", "dplyr", "caret", "randomForest", "ggplot2")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(required_packages, install_if_missing))

library(data.table)
library(dplyr)
library(caret)
library(randomForest)
library(ggplot2)


# ==========================================================
# 2. READ INPUT DATA
# ==========================================================
gene_df   <- fread("gene_expression.csv", data.table = FALSE)
microbe_df <- fread("microbiome_abundance.csv", data.table = FALSE)
meta_df   <- fread("metadata.csv", data.table = FALSE)


# ==========================================================
# 3. DATA CLEANING AND FORMATTING
# ==========================================================
# Assign rownames (feature IDs)
rownames(gene_df) <- gene_df[, 1]
gene_df <- gene_df[, -1, drop = FALSE]

rownames(microbe_df) <- microbe_df[, 1]
microbe_df <- microbe_df[, -1, drop = FALSE]

# Convert to numeric matrices
gene_df <- as.data.frame(lapply(gene_df, as.numeric), row.names = rownames(gene_df))
microbe_df <- as.data.frame(lapply(microbe_df, as.numeric), row.names = rownames(microbe_df))


# ==========================================================
# 4. ALIGN SHARED SAMPLES
# ==========================================================
common_samples <- Reduce(intersect,
                         list(colnames(gene_df),
                              colnames(microbe_df),
                              meta_df$SampleID))

if (length(common_samples) < 5) {
  stop("ERROR: Too few shared samples across datasets.")
}

gene_df   <- gene_df[, common_samples, drop = FALSE]
microbe_df <- microbe_df[, common_samples, drop = FALSE]
meta_df   <- meta_df %>% filter(SampleID %in% common_samples)

meta_df <- meta_df[match(common_samples, meta_df$SampleID), ]


# ==========================================================
# 5. GENE EXPRESSION PREPROCESSING
# ==========================================================
gene_log <- log2(gene_df + 1)

gene_var <- apply(gene_log, 1, var, na.rm = TRUE)
gene_filt <- gene_log[gene_var > quantile(gene_var, 0.75), ]


# ==========================================================
# 6. MICROBIOME PREPROCESSING
# ==========================================================
# Relative abundance
microbe_rel <- sweep(microbe_df, 2, colSums(microbe_df), "/")

# Handle zeros
microbe_rel[microbe_rel == 0] <- 1e-6

# Log transformation (CLR-like)
microbe_log <- log(microbe_rel)

# Variance filtering
microbe_var <- apply(microbe_log, 1, var)
microbe_filt <- microbe_log[microbe_var > quantile(microbe_var, 0.50), ]


# ==========================================================
# 7. TRANSPOSE TO SAMPLE × FEATURE
# ==========================================================
gene_mat   <- t(as.matrix(gene_filt))
microbe_mat <- t(as.matrix(microbe_filt))

colnames(gene_mat)   <- paste0("GENE_", make.names(colnames(gene_mat)))
colnames(microbe_mat) <- paste0("MICROBE_", make.names(colnames(microbe_mat)))


# ==========================================================
# 8. STANDARDIZATION
# ==========================================================
gene_scaled   <- scale(gene_mat)
microbe_scaled <- scale(microbe_mat)

gene_scaled[is.na(gene_scaled)] <- 0
microbe_scaled[is.na(microbe_scaled)] <- 0


# ==========================================================
# 9. EARLY DATA FUSION
# ==========================================================
fused_mat <- cbind(gene_scaled, microbe_scaled)
fused_df  <- as.data.frame(fused_mat)
fused_df$SampleID <- rownames(fused_df)

final_df <- meta_df %>%
  inner_join(fused_df, by = "SampleID")

write.csv(final_df, "fused_feature_matrix.csv", row.names = FALSE)


# ==========================================================
# 10. DEFINE OUTCOME
# ==========================================================
if (!"Group" %in% colnames(final_df)) {
  stop("ERROR: Metadata must contain 'Group' column.")
}

final_df$Group <- as.factor(final_df$Group)


# ==========================================================
# 11. TRAIN / TEST SPLIT
# ==========================================================
set.seed(123)
train_idx <- createDataPartition(final_df$Group, p = 0.8, list = FALSE)

train_df <- final_df[train_idx, ]
test_df  <- final_df[-train_idx, ]

x_train <- train_df %>% select(-SampleID, -Group)
y_train <- train_df$Group

x_test <- test_df %>% select(-SampleID, -Group)
y_test <- test_df$Group


# ==========================================================
# 12. RANDOM FOREST MODEL
# ==========================================================
set.seed(123)

rf_model <- randomForest(
  x = x_train,
  y = y_train,
  ntree = 500,
  importance = TRUE
)

print(rf_model)


# ==========================================================
# 13. MODEL EVALUATION
# ==========================================================
pred <- predict(rf_model, newdata = x_test)

conf_mat <- confusionMatrix(pred, y_test)
print(conf_mat)

pred_df <- data.frame(
  SampleID = test_df$SampleID,
  Actual = y_test,
  Predicted = pred
)

write.csv(pred_df, "train_test_predictions.csv", row.names = FALSE)


# ==========================================================
# 14. FEATURE IMPORTANCE
# ==========================================================
imp <- importance(rf_model)

imp_df <- data.frame(
  Feature = rownames(imp),
  MeanDecreaseGini = imp[, "MeanDecreaseGini"]
) %>%
  arrange(desc(MeanDecreaseGini))

write.csv(imp_df, "feature_importance.csv", row.names = FALSE)


# Plot top features
top_imp <- head(imp_df, 20)

ggplot(top_imp, aes(x = reorder(Feature, MeanDecreaseGini),
                    y = MeanDecreaseGini)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Top 20 Features (Random Forest)",
       x = "Feature",
       y = "Importance")

ggsave("top20_feature_importance.png", width = 8, height = 6)


# ==========================================================
# 15. PCA VISUALIZATION
# ==========================================================
pca_res <- prcomp(final_df %>% select(-SampleID, -Group), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = final_df$Group
)

ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA of Early-Fused Data")

ggsave("pca_early_fusion.png", width = 7, height = 5)


cat("✔ Pipeline completed successfully.\n")