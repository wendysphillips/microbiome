library(tidyverse)

# Loading data from HMP2 ----

# Data tables originated from https://ibdmdb.org/tunnel/public/summary.html

# The metadata table was already partially filtered from full table
metadata_full <- read.delim("data/hmp2_metadata_partial.csv", sep = ",", check.names = FALSE)

# Get rid of spaces in column names
colnames(metadata_full) <- gsub(" ", "_", colnames(metadata_full))

# This is the MetaPhlAn3 taxatable
taxa <- read.delim("data/hmp2_taxonomic_profiles_3.tsv", sep = "\t", check.names = FALSE)

# Remove superfluous part of sample names
colnames(taxa) <- gsub("_profile", "", colnames(taxa))

# Check the same samples are in metadata and taxatable
samps <- colnames(taxa)[2:1639]
all.equal(sort(samps), sort(metadata_full$`External ID`))

# Remove uninformative metadata -----

# Find how many unique values for each metadata column
temp <- c()
for (col in 1:ncol(metadata_full)){
      temp <- c(temp, length(unique(metadata_full[,col])))
}

# Filter out any metadata without any variation
metadata <- metadata_full[,temp > 1]
rm(temp)

# Create table with number of unique values
#  and ratio of smallest to largest value
temp_ratio <- c()
temp_n <- c()
for (col in 1:ncol(metadata)) {
      temp_min <- min(table(metadata[,col]))
      temp_max <- max(table(metadata[,col]))
      temp_ratio <- c(temp_ratio, temp_min/temp_max)
      temp_n <- c(temp_n, nrow(table(metadata[,col])))
      
      rm(temp_min, temp_max)
}

# Put the above into a datafame with metadata names
temp_meta_info <- data.frame(meta_name = colnames(metadata), n_vals = temp_n, ratio = temp_ratio)
rm(temp_n, temp_ratio)

# First filter, keep things with > 2 values
#  Of if just two values, the ratio indicates reasonable number in less abundant type
keep <- temp_meta_info[(temp_meta_info$n_vals > 2) | ((temp_meta_info$n_vals == 2) & (temp_meta_info$ratio > 0.03)),]
metadata <- metadata[,keep$meta_name]

# Out of curiousity, manually inspect the distribution of values in remaining columns
#  for those columns with no more than 10 levels
for (col in 1:ncol(metadata)){
      z <- table(metadata[,col])
      if (nrow(z) < 11) {
            print(z)
      }
      rm(z)
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
metadata <- metadata[metadata$External_ID %in% colnames(species),]

# Tidy up
rm(investigate, keep,temp_meta_info, col, samps)
