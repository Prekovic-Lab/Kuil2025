Subset <- ScaleData(Subset, features = all_genes)


goi <- c("Pdgfra")
expression <- FetchData(Subset, vars = goi)
expression <- as.data.frame(expression)
expression$cell <- Subset$cell
expression$Condition <- Subset$Condition
expression$Type_sup_all <- Subset$Type_sup_all
expression$Type_prob_all <- Subset$Type_prob_all
expression$Cat <- Subset$Cat  

# Filter the dataframe for "Oligodendrocyte precursor cells"
filtered_df <- expression[expression$Type_sup_all == "Oligodendrocyte precursor cells", ]

#Does take the 0's into account. So the mean of y is lower, but that's due to more 0 values in there. See bootstrapping
means <- aggregate(Pdgfra ~ Condition, data = filtered_df, FUN = mean)

# Perform a t-test (or Wilcoxon test if data is not normally distributed) between the two conditions
stat_test <- wilcox.test(Pdgfra ~ Condition, data = filtered_df)

custom_colors <- c("n" = "lightblue", "y" = "salmon")  # Replace with your conditions and desired colors

# Generate the plot
ggplot(filtered_df, aes(x = Condition, y = Pdgfra)) +
  geom_violin(trim = FALSE, aes(fill = Condition)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Boxplot inside violin plot
  geom_jitter(width = 0.2, alpha = 0.5) +  # Show individual data points with jitter
  stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.25) +  # Show the mean as a red diamond
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), hjust = 1.8, vjust = -1.5) +  # Show the mean as text
  scale_y_continuous(expand = c(0, 0)) +  # Ensure the y-axis starts at 0
  scale_fill_manual(values = custom_colors) +  # Apply custom colors for the violins
  ggtitle("Gene Expression of Pdgfra per Condition in OPCs") +
  labs(y = "Scaled Pdgfra expression") +
  theme_minimal() +
  stat_compare_means(method = "wilcox.test", label.y = max(filtered_df$Pdgfra), hjust = -1)  # Add test result

