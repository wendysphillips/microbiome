library(tidyverse)

# First run the filter_hmp_data.R code
source("filter_hmp_data.R")

mycolors <- c(
  "#E4A6A6", "#D8908F", "#BF6C6B", "#DF6F6D", "#F96B68", "#DD3933", "#F3E1B4", "#DDC88A", "#AA944B", "#DCBC41", "#F3CF36", "#D5E6B9", "#AEC38A",
  "#8EA75F", "#D3F496", "#B6DC66", "#9DC71B", "#ABEE9A", "#78C960", "#50B61B", "#C3E3E1", "#8CBFBD", "#4AA19D",
  "#B2CFEB", "#82ADD1", "#3F82AE", "#B9BBE4", "#9C9FD2", "#7D80B0", "#C0B3D7", "#AC96D0", "#8A66BB", "#EAB5EF", "#D38DD9", "#CC6CD3", "#E696B9",
  "#DB75A4", "#D1398C"
)
set.seed(119)
my_colors_random <- sample(mycolors, length(mycolors), replace = FALSE)


# Create function to summarize taxa to different levels
get_specific_level <- function(input, level = "Genus", sep = ";") {
  temp <- input |> tibble::rownames_to_column("taxa")
  temp2 <- temp |>
    tidyr::separate(taxa,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
      extra = "merge",
      fill = "right",
      sep = sep
    )

  temp3 <- temp2 |> tidyr::unite("taxa", grep("Kingdom", colnames(temp2)):grep(level, colnames(temp2)), sep = ";")
  temp3 <- temp3 |> dplyr::select(-any_of(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")))
  temp3 <- temp3 |>
    dplyr::group_by(taxa) %>%
    summarize_all(list(sum))
  temp3 <- temp3 |> tibble::column_to_rownames("taxa")
}

# Get genus level abundances
genus_table <- get_specific_level(species, "Genus", sep = "\\|")

# Shorten name to just genus identifier
rownames(genus_table) <- gsub(".*;g__", "", rownames(genus_table))

# Get IDs for all subjects to use for subsetting
subject_list <- unique(metadata$Participant_ID)

# Subset to just samples from the first subject in the list
genus_s2 <- genus_table[, metadata$External_ID[metadata$Participant_ID == subject_list[4]]]

# Save only rows for taxa that have a minimum relative abundance of 2% in at least one sample
genus_s2_filt <- genus_s2[apply(genus_s2, 1, max) > 2, ]

# Save taxa names to factor later
taxa_names_kept <- rownames(genus_s2_filt)

# Get the total abundance of taxa not included
genus_s2_filt_other <- genus_s2[!rownames(genus_s2) %in% taxa_names_kept, ]

# Add that 'Other' value to plotting values so that they will total 100
genus_s2_filt <- rbind(genus_s2_filt, colSums(genus_s2_filt_other))
rownames(genus_s2_filt)[nrow(genus_s2_filt)] <- "Other"

# Transform table to long format required for plotting
genus_long <- genus_s2_filt %>%
  rownames_to_column("Genus") %>%
  pivot_longer(cols = 2:ncol(genus_long))

# Factor the genus names so Other is separated out
genus_long$Genus <- factor(genus_long$Genus, levels = c("Other", rev(taxa_names_kept)))

# Produce stacked bar plot
ggplot(genus_long, aes(x = name, y = value, fill = Genus)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_manual(values = sample(my_colors_random, length(unique(genus_long$Genus))))
