# Script to create normalized mean PoPS score across 111 UKB phenotypes
# Author: Nate
# Date: March 12 2025

# install libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)

# load libraries
ukb.phenos.pops = fread("~/Desktop/PoPS_FullResults.txt")

# Create an empty list to store trait-specific subsets
trait_list <- list()

# Split data by trait
trait_list <- split(ukb.phenos.pops, ukb.phenos.pops$trait)

# Apply normalization within each trait
trait_list <- lapply(trait_list, function(df) {
  df <- df %>%
    mutate(pops_norm = as.numeric(scale(pops_score, center = TRUE, scale = TRUE))) # Normalize
  return(df)
})

# Combine all trait-specific data frames into one long-format data frame
combined_df = bind_rows(trait_list, .id = "trait")

# Compute summary statistics per gene across all traits
gene_summary <- combined_df %>%
  group_by(gene) %>%
  summarise(
    mean_pops_norm = mean(pops_norm, na.rm = TRUE),
    median_pops_norm = median(pops_norm, na.rm = TRUE),
    sd_pops_norm = sd(pops_norm, na.rm = TRUE)
  )

# Compute summary statistics per gene across all traits
ensembl_summary <- combined_df %>%
  group_by(ensgid) %>%
  summarise(
    mean_pops_norm = mean(pops_norm, na.rm = TRUE),
    median_pops_norm = median(pops_norm, na.rm = TRUE),
    sd_pops_norm = sd(pops_norm, na.rm = TRUE)
  )

# Save using fwrite with tab separator
fwrite(ensembl_summary, "~/drugsets/DATA/COV/ensembl_mean_pops.cov", sep = "\t", quote = FALSE)
fwrite(gene_summary, "~/drugsets/DATA/COV/symbol_mean_pops.cov", sep = "\t", quote = FALSE)










