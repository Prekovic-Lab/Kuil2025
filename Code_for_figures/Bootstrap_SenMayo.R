expression_data <- GetAssayData(Subset, slot = "data")

SenMayo <- genes <- c(
  "Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", "Bmp2", "Bmp6", 
  "C3", "Ccl1", "Ccl2", "Ccl20", "Ccl24", "Ccl26", "Ccl3", "Ccl4", "Ccl5", 
  "Ccl7", "Ccl8", "Cd55", "Cd9", "Csf1", "Csf2", "Csf2rb", "Cst10", "Ctnnb1", 
  "Ctsb", "Cxcl1", "Cxcl10", "Cxcl12", "Cxcl16", "Cxcl2", "Cxcl3", "Cxcr2", 
  "Dkk1", "Edn1", "Egf", "Egfr", "Ereg", "Esm1", "Ets2", "Fas", "Fgf1", 
  "Fgf2", "Fgf7", "Gdf15", "Gem", "Gmfg", "Hgf", "Hmgb1", "Icam1", "Icam5", 
  "Igf1", "Igfbp1", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6", "Igfbp7", 
  "Il10", "Il13", "Il15", "Il18", "Il1a", "Il1b", "Il2", "Il6", "Il6st", 
  "Il7", "Inha", "Iqgap2", "Itga2", "Itpka", "Jun", "Kitl", "Lcp1", "Mif", 
  "Mmp13", "Mmp10", "Mmp12", "Mmp14", "Mmp2", "Mmp3", "Mmp9", "Nap1l4", 
  "Nrg1", "Pappa", "Pecam1", "Pgf", "Pigf", "Plat", "Plau", "Plaur", "Ptbp1", 
  "Ptger2", "Ptges", "Rps6ka5", "Scamp4", "Selplg", "Sema3f", "Serpinb3a", 
  "Serpine1", "Serpine2", "Spp1", "Spx", "Timp2", "Tnf", "Tnfrsf11b", "Tnfrsf1a", 
  "Tnfrsf1b", "Tubgcp2", "Vegfa", "Vegfc", "Vgf", "Wnt16", "Wnt2")

#Senescence and quiescence genes as defined by Nanostring 
SenNANO <- c(
  "Ago3", "Anapc16", "Atm", "Atr", "Btrc", "Cacna1d", "Calm1", "Calm2", "Calm3",
  "Ccl2", "Cdc27", "Cdk6", "Cxcl2", "Fbxw11", "Foxo1", "Foxo3", "H3f3b", "Hipk2",
  "Id1", "Id2", "Igfbp7", "Il6", "Itpr1", "Itpr2", "Map2k4", "Map2k6", "Map3k5",
  "Map4k4", "Mapk14", "Mapk8", "Mcu", "Nfkb1", "Phc2", "Phc3", "Ppp3ca", "Pten",
  "Raf1", "Rb1", "Sirt1", "Slc25a4", "Smad2", "Smad3", "Stat3", "Terf2ip", "Tfdp2",
  "Tgfb1", "Tgfb2", "Tgfbr1", "Tnik", "Tnrc6a", "Tnrc6b", "Tnrc6c", "Trp53",
  "Trpm7", "Ube2e1"
)

pathways <- list("Senes" = SenMayo, "SenesNANO" = SenNANO)

filtered_genes <- all_genes[!grepl("Negative", all_genes, ignore.case = TRUE)]
set.seed(123)  # Set seed for reproducibility
sample_list <- list()

# Perform the sampling 1000 times

## 31 genes SenMayo, 55 genes SenNANO
for (i in 1:1000) {
  sampled_genes <- sample(filtered_genes, 31, replace = FALSE)  # Sample 31 genes without replacement
  sample_list[[paste0("sample_", i)]] <- sampled_genes  # Name each element as sample_1, sample_2, etc.
}
#print(sample_list)

pathways <- c(pathways, sample_list)

gsva_results <- gsva(expression_data, pathways, method = "ssgsea")

# Convert GSVA results into a data frame and add to Seurat metadata
gsva_results_transposed <- t(gsva_results)
gsva_results_df_Mayo <- as.data.frame(gsva_results_transposed)

gsvea_results_df_Mayo[,2] <- Subset$Condition

plot_data <- data.frame(
  Value = c(gsva_results_df_Mayo[, 1], average_distribution),
  Distribution = factor(c(as.character(gsva_results_df_Mayo[, 2]), 
                          rep("Average Distribution", length(average_distribution))))
)

# Plot the distributions
ggplot(plot_data, aes(x = Value, fill = Distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Column 1 (Split by Column 2) vs. Average Distribution",
       x = "Value",
       y = "Density") +
  scale_fill_manual(values = c("blue", "green", "red")) +  # Adjust colors as needed
  theme_minimal() +
  theme(legend.position = "top")

ggsave("Bootstrap_SenMayo_n31.svg",plot = last_plot())
write.csv(gsva_results_df_Mayo, "gsva_results_df_Mayo.csv")