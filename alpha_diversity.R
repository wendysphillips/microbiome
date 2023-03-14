# Load these libraries needed to run this code
library(tidyverse)
library(vegan)
library(cowplot)
library(ggpubr)


# import a taxatable and a metadata table
taxatable <- read.delim("data/zackular_taxatable.csv", sep = ",", check.names = FALSE, row.names = 1)

metadata <- read.delim("data/zackular_metadata.csv", sep = ",", check.names = FALSE, row.names = 1)

# Quick check that samples are in columns and ballpark of total read counts 
colSums(taxatable)[1:5]

# Check what the minimim read count is in sample with fewest reads
# If > 10,000 & amplicon data, reasonable to use this value to rarefy
# If minimum < 10,000, consider dropping low count samples
min(colSums(taxatable))

# Transpose the dataframe for use with vegan
# I add the trailing _t to the name because I consider
# standard taxatable format to be samples in columns, taxa in rows
taxa_rarefied_t <- data.frame(vegan::rrarefy(t(taxatable), min(colSums(taxatable))), check.names = FALSE)

# Quick peak shows this has samples in rows and taxa in columns
taxa_rarefied_t[1:2,1:2]

# Check that all samples now have same read count
rowSums(taxa_rarefied_t)[1:5]

# Remove any taxa that might have all 0 read counts after rarefaction
taxa_rarefied_t <- taxa_rarefied_t[, colSums(taxa_rarefied_t) != 0]

# Get alpha diversity - number of taxa
num_taxa <- vegan::specnumber(taxa_rarefied_t)

# Get alpha diversity - Shannon
shannon <- vegan::diversity(taxa_rarefied_t, index = "shannon", MARGIN = 1)

# Get alpha diversity - Inverse Simpson
simpson <- vegan::diversity(taxa_rarefied_t, index = "invsimpson", MARGIN = 1)

# Get alpha diversity - Chao1
chao1 <- vegan::estimateR(taxa_rarefied_t)[2, ]

# Combine all alpha metrics into a single dataframe
alpha_values <- data.frame(Num_taxa = num_taxa, Chao1 = chao1, Shannon = shannon, Simpson = simpson)

# Merge alpha values with metadata
alpha_values <- merge(metadata, alpha_values, by = 0, )

# Rename the new Row.names column created by merge
alpha_values <- alpha_values %>% dplyr::rename("Sample" = "Row.names")

# Basic boxplot comparing Shannon diversity index values between disease states
shannon_plot <- ggplot(alpha_values, aes(x = DiseaseState, y = Shannon)) +
  geom_boxplot(outlier.color = NA, lwd = 0.5) +
  theme_cowplot() +
  geom_jitter(size = 3, width = 0.2, alpha = 0.4, aes(color = DiseaseState)) +
  scale_color_manual(values = c("#7153db", "#008eb0", "magenta")) +
  theme(legend.position = "none") +
  labs(x = "Group")

shannon_plot

# Create a function to make a similar plot with any of the metrics
alpha_plot <- function(input_table, metric, meta = "DiseaseState") {
  ggplot(input_table, aes_string(x = meta, y = metric)) +
    geom_boxplot(outlier.color = NA, lwd = 0.5) +
    theme_cowplot() +
    geom_jitter(size = 3, width = 0.2, alpha = 0.4, aes(color = DiseaseState)) +
    scale_color_manual(values = c("#7153db", "#008eb0", "#cd53db")) +
    theme(legend.position = "none") +
    labs(x = "Group")
}

# Now can easily create the same plot for alpha diversity
alpha_plot(alpha_values, "Chao1")

# Or number of taxa
alpha_plot(alpha_values, "Num_taxa") + labs(y = "Number of taxa")

# Perform a wilcoxon rank sum test between two disease states
wilcox.test(alpha_values$Shannon[alpha_values$DiseaseState == "H"], alpha_values$Shannon[alpha_values$DiseaseState == "CRC"])

# Create a function that will include the results of a wilcoxon rank sum test
# This function is for a metadata factor that has only two levels
alpha_plot_stat <- function(input_table, metric, meta = "DiseaseState") {
  # Do Wilcoxon rank sum test
  temp <- as.formula(paste0(metric, " ~ ", meta))
  wilc <- input_table %>% rstatix::wilcox_test(temp)

  # Plot with p-value
  g <- ggplot(input_table, aes_string(x = "DiseaseState", y = metric)) +
    geom_boxplot(outlier.color = NA, lwd = 0.5) +
    theme_cowplot() +
    geom_jitter(size = 3, width = 0.2, alpha = 0.4, aes(color = DiseaseState)) +
    scale_color_manual(values = c("#7153db", "#008eb0", "#cd53db")) +
    theme(legend.position = "none") +
    labs(x = "Group") 
  
  g <- g + ggpubr::stat_pvalue_manual(wilc,
      label = "p",
      y.position = 1.1 * max(alpha_values[, metric]),
      step.increase = 0.1
    )
  g
}

# Chao1 with wilcoxon-rank sum test p-values
alpha_plot_stat(alpha_values, "Chao1")

# Shannon with wilcoxon-rank sum test p-values
alpha_plot_stat(alpha_values, "Shannon")
