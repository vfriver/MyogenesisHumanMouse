library(pheatmap)

# Compute log2(TPM + 1)
log_tpm_human <- log2(counts_to_tpm(human_counts, human_lengths) + 1)
log_tpm_mouse <- log2(counts_to_tpm(mouse_counts, mouse_lengths) + 1)

# Combine human and mouse matrices
log_tpm_combined <- cbind(log_tpm_human, log_tpm_mouse)

# Save matrix
write.csv(log_tpm_combined, "/home/log2TPM_combined.csv")


# Basic heatmap
pheatmap(
  log_tpm_combined,
  scale = "none",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",  
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 3,
  fontsize_col = 3,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)

# Column annotation for species
annotation_col <- data.frame(
  Species = rep(c("Human", "Mouse"), c(4,4))
)
rownames(annotation_col) <- colnames(log_tpm_combined)

# Heatmap color palette
heat_colors <- colorRampPalette(c("navy", "khaki1", "firebrick3"))(100)

# Final clustered heatmap
pheatmap(
  log_tpm_combined,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",  
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 3,
  fontsize_col = 3,
  color = heat_colors
)
