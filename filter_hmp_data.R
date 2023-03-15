library(tidyverse)

# Loading data from HMP2 ----

# The metadata table was already partially filtered from full table
metadata <- read.delim("data/hmp2_metadata_partial.csv", sep = ",", check.names = FALSE)

# Get rid of spaces in column names
colnames(metadata) <- gsub(" ", "_", colnames(metadata))

# This is the MetaPhlAn3 taxatable
taxa <- read.delim("data/hmp_taxonomic_profiles_3.tsv", sep = "\t", check.names = FALSE)

# Remove superfluous part of sample names
colnames(taxa) <- gsub("_profile", "", colnames(taxa))

# Check the same samples are in metadata and taxatable
samps <- colnames(taxa)[2:1639]
all.equal(sort(samps), sort(metadata$`External ID`))

# Remove uninformative metadata -----

# Find how many unique values for each metadata column
temp <- c()
for (col in 1:ncol(metadata)){
      temp <- c(temp, length(unique(metadata[,col])))
}

# Filter out any metadata without any variation
meta <- metadata[,temp > 1]

# Create table with number of unique values
#  and ratio of smallest to largest value
temp_ratio <- c()
temp_n <- c()
for (col in 1:ncol(meta)) {
      temp_min <- min(table(meta[,col]))
      temp_max <- max(table(meta[,col]))
      temp_ratio <- c(temp_ratio, temp_min/temp_max)
      temp_n <- c(temp_n, nrow(table(meta[,col])))
}

# Put the above into a datafame with meta names
temp_df <- data.frame(meta_name = colnames(meta), n_vals = temp_n, ratio = temp_ratio)

# First filter, keep things with > 2 values
#  Of if just two values, the ratio indicates reasonable number in less abundant type
keep <- temp_df[(temp_df$n_vals > 2) | ((temp_df$n_vals == 2) & (temp_df$ratio > 0.03)),]
meta <- meta[,keep$meta_name]

# Out of curiousity, manually inspect the distribution of values in remaining columns
#  for those columns with no more than 10 levels
for (col in 1:ncol(meta)){
      z <- table(meta[,col])
      if (nrow(z) < 11) {
            print(z)
      }
}

# Taxatable manipulation ------

# Change the name of the taxa column to be more R friendly
colnames(taxa)[1] <- "taxa"

# Subset to species data
species <- taxa[grep("s__", taxa$taxa),]

# Move species strings to rownames
species <- species %>% remove_rownames() %>% column_to_rownames("taxa")

# Check that all samples have relative abundances that total ~100
sum( (colSums(species) > 98) & (colSums(species) < 102 ) ) == ncol(species)

# Some don't!
# Investigate those samples that don't sum to ~100
investigate <- species[,!((colSums(species) > 98) & (colSums(species) < 102 ) )]     
colSums(investigate)
# These samples have a total relative abundance of 0!

# Will drop samples without species abundances
species <- species[,((colSums(species) > 98) & (colSums(species) < 102 ) )] 

# Also drop those from metadata
meta <- meta[meta$External_ID %in% colnames(species),]

