# This is just the start of this code. More to come!

# Load libraries that will be used
library(vegan)
library(ecodist)

# First run the filter_hmp_data.R code
# Can be done by uncommenting next line
# source("filter_hmp_data.R")

BC <- vegan::vegdist(t(species), method = "bray")

# Get pcoa
BC_pcoa <- ecodist::pco(BC)

# Get info on first two pcoa axes
BC_pcoa_df <- data.frame(
  pcoa1 = BC_pcoa$vectors[,1],
  pcoa2 = BC_pcoa$vectors[,2]
)

pc_plot <- ggplot(BC_pcoa_df, aes(x = pcoa1, y = pcoa2)) +
  geom_point(alpha = 0.5) +
  theme_classic()

pc_plot
