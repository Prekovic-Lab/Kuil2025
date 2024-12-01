result_supers <- list()
#Define all genes you need. There are negative controls in here you want to exclude since they increase padj 
all_genes <- rownames(Subset)
#Filter out names starting with "Negative"
genes <- all_genes[!grepl("^Negative", all_genes)]

Idents(Subset) <- Subset$Condition
# Iterate over selected_colnames
for (i in selected_colnames) {
  # Subset the Seurat object based on Supergroup
  subset_Subset <- subset(Subset, Cat == i)
  print(i)
  
  # Initialize an empty data frame to store results for the current iteration
  result_super <- data.frame()
  
  # Perform FindMarkers command with error handling
  tryCatch({
    DEG1 <- FindMarkers(subset_Subset, ident.1 = "y", ident.2 = "n", slot = "data", min.pct = 0, logfc.threshold = 0, features = genes, test.use = "MAST")
    DEG1$i <- i  # Add Supergroup information
    DEG1$ident1 <- "y"
    DEG1$ident2 <- "n"
    result_super <- bind_rows(result_super, DEG1)
  }, error = function(e) {
    message(paste("Skipping due to error in FindMarkers for", i, ":", e$message))
  })
  
  # Append result_super to result_supers if it's not empty
  if (nrow(result_super) > 0) {
    result_supers <- bind_rows(result_supers, result_super)
  }
}

# Extract gene names from row names
result_supers$Gene <- gsub("\\.\\.\\..*", "", rownames(result_supers))

DEG <- subset(result_supers, i == "Astrocytes")

#Plotting code is the same as DEG_ALL except for xlim
plot <- ggplot(DEG, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5), size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(0, 0.5, -0.5), linetype = "dashed", color = c("black", "gray", "gray")) +
  geom_text_repel(
    data = subset(DEG, abs(avg_log2FC) > 0.5 & p_val_adj < 0.05), 
    aes(label = Gene), 
    hjust = -0.1, vjust = 0
  ) +
  labs(x = "Average log2 Fold Change", y = "-log10(Adjusted p-value)",
       title = "Volcano Plot in Astrocytes irradiated over non-irradiated") +
  theme_minimal() +
  NoLegend()

ggsave("DEG_Astrocytes.svg", plot = last_plot())