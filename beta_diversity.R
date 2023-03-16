# This is just the start of this code. More to come!

# Load libraries that will be used
library(vegan)
library(ape)

# First run the filter_hmp_data.R code
# Can be done by uncommenting next line
source("filter_hmp_data.R")

# Calculate Bray-Curtis dissimilarity matrix
BC <- vegan::vegdist(t(species), method = "bray")

# Get pcoa
BC_pcoa <- ape::pcoa(BC)

# Get coordinates first first two pcoa axes
BC_pcoa_df <- data.frame(
  pcoa1 = BC_pcoa$vectors[, 1],
  pcoa2 = BC_pcoa$vectors[, 2]
)

# Simple plot
pc_plot <- ggplot(BC_pcoa_df, aes(x = pcoa1, y = pcoa2)) +
  geom_point(alpha = 0.5) +
  theme_classic()

pc_plot

# Combine with metadata
BC_pcoa_df$External_ID <- colnames(species)
BC_pcoa_all <- left_join(BC_pcoa_df, metadata, by = "External_ID")

# Plot with metadata variable
pc_plot <- ggplot(BC_pcoa_all, aes(x = pcoa1, y = pcoa2)) +
  geom_point(alpha = 0.5, aes(color = diagnosis)) +
  theme_classic() +
  labs(x = paste0("PCoA1 ", round(BC_pcoa$values$Eigenvalues[1], 1), " %"), y = paste0("PCoA2 ", round(BC_pcoa$values$Eigenvalues[2], 1), " %"))

pc_plot

# Make the above into a function
#  to input different metadata
pcoa_plot <- function(meta = "Chemotherapy") {
  pc_plot <- ggplot(BC_pcoa_all, aes(x = pcoa1, y = pcoa2)) +
    geom_point(alpha = 0.5, aes_string(color = meta)) +
    theme_classic() +
    labs(x = paste0("PCoA1, ", round(BC_pcoa$values$Relative_eig[1], 2), "%"), y = paste0("PCoA2, ", round(BC_pcoa$values$Relative_eig[2], 2), " %"))

  pc_plot
}

# Plot with a categorical variable
pcoa_plot("sex")

# Plot with a continuous metadata variable
pcoa_plot("Age_at_diagnosis")

# Change legend title
pcoa_plot("site_name") + guides(color = guide_legend("Site"))
