# Create drug gene sets
# May 28 2026
# Nate YB

# ---- (0) Load packages ----

# load bioconductor and database packages first
library(biomaRt)

# load data handling and generic packages
library(arrow)
library(data.table)

# load tidyverse last so its functions take absolute priority
library(purrr)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)

# ---- (1) Define functions ---- 

write.magma.geneset = function(df, set.id, gene.id, file.path) {
  
  # subset, merge duplicate sets, and deduplicate genes
  df.sub = df %>%
    dplyr::select(all_of(c(set.id, gene.id))) %>%
    # explicitly remove unmapped/empty set IDs before grouping
    filter(!is.na(.data[[set.id]]) & str_trim(.data[[set.id]]) != "") %>%
    group_by(.data[[set.id]]) %>%
    summarise(
      !!gene.id := list(unique(na.omit(unlist(.data[[gene.id]])))),
      .groups = "drop"
    ) %>%
    filter(lengths(.data[[gene.id]]) >= 5)
  
  # -- error checking -- 
  
  # duplicate set ids
  if (any(duplicated(df.sub[[set.id]]))) {
    stop("QC FAIL: Duplicate set IDs detected in '", set.id, "'.")
  }
  
  # duolicate genes inside sets
  duplicate_genes_check = map_lgl(df.sub[[gene.id]], ~ any(duplicated(.x)))
  if (any(duplicate_genes_check)) {
    stop("QC FAIL: Duplicate genes detected within individual sets for '", set.id, "'.")
  }
  
  # checking if data is left after filtering
  if (nrow(df.sub) == 0) {
    stop("QC FAIL: No valid gene sets with > 1 gene remain for '", set.id, "'.")
  }
  
  # collapse lists into tab-separated
  lines = map2_chr(df.sub[[set.id]], df.sub[[gene.id]], function(name, genes) {
    # remove blanks
    clean.genes = genes[genes != ""]
    # combine with set id using tabs for both the separator and the collapse
    paste(name, paste(clean.genes, collapse = "\t"), sep = "\t")
  })
  
  # compile background set
  all.genes = unique(unlist(df.sub[[gene.id]]))
  all.genes = all.genes[all.genes != ""]
  
  # format background set using tabs
  druggable.line = paste("druggable", paste(all.genes, collapse = "\t"), sep = "\t")
  
  # append background set
  lines = c(lines, druggable.line)
  
  # save
  writeLines(lines, con = file.path)
  
  # print
  message("Saved ", length(lines) - 1, " standard sets + 1 'druggable' background set to: ", file.path)
}

# helper to extract the highest clinical stage from a vector
get.max.stage = function(stages) {
  stage_levels = c(
    "UNKNOWN", "PRECLINICAL", "IND", "EARLY_PHASE_1",
    "PHASE_1", "PHASE_1_2", "PHASE_2", "PHASE_2_3",
    "PHASE_3", "PREAPPROVAL", "APPROVAL"
  )
  stages = na.omit(stages)
  if(length(stages) == 0) {
    return("UNKNOWN")
  }
  ordered_stages = factor(stages, levels = stage_levels, ordered = TRUE)
  as.character(max(ordered_stages))
}

# helper function to handle list summarization consistently
collapse_unique = function(x) {
  x = unique(na.omit(x))
  if (length(x) == 0) return(NULL) # NULLs are safer than empty lists in tibbles
  return(x)
}

# ---- (2) Load data ---- 

# load master data files
load("/Users/nyb/phd_code/drugsets-may2026/data/master-data-files-2026/master_long_2026.RData")

# ---- (3) Save DRUG gene sets ---- 

meta.list = stack.long.clean %>%
  group_by(master_id) %>%
  arrange(nchar(drug_name_raw), .by_group = TRUE) %>%
  summarise(
    # --- flat identifiers ---
    harmonized_chembl_id = first(na.omit(harmonized_chembl_id)),
    drug_name            = first(na.omit(drug_name)),
    drug_name_raw        = first(na.omit(drug_name_raw)),
    
    # --- clinical stage ---
    max_clinical_stage   = if(all(is.na(max_clinical_stage))) NA_character_ else get.max.stage(na.omit(max_clinical_stage)),
    
    # --- list-columns ---
    ensembl         = list(collapse_unique(ensembl)),
    entrez          = list(collapse_unique(entrez)),
    symbol          = list(collapse_unique(symbol)),
    moa_chembl      = list(collapse_unique(moa_chembl)),
    moa_targ_chembl = list(collapse_unique(moa_targ_chembl)),
    moa_clue        = list(collapse_unique(moa_clue)),
    moa_targ_clue   = list(collapse_unique(moa_targ_clue)),
    moa_ot          = list(collapse_unique(moa_ot)),
    moa_targ_ot     = list(collapse_unique(moa_targ_ot)),
    atc1            = list(collapse_unique(atc1)),
    atc2            = list(collapse_unique(atc2)),
    atc3            = list(collapse_unique(atc3)),
    atc4            = list(collapse_unique(atc4)),
    indication_efo  = list(collapse_unique(indication_efo)),
    
    .groups = "drop"
  ) %>%
  mutate(
    set.name = coalesce(drug_name, master_id)
  ) %>%
  dplyr::select(set.name, everything())

# check 
nrow(meta.list) == length(unique(stack.long.clean$master_id))

# subset to approved drugs
meta.list.apprv = meta.list %>%
  filter(lengths(ensembl) > 1) %>%
  filter(toupper(max_clinical_stage) == "APPROVAL")

# save
write.magma.geneset(
  df = meta.list.apprv, 
  set.id = "set.name", 
  gene.id = "ensembl", 
  file.path = "/Users/nyb/drugsets/DATA/GENESETS/solo_ensembl_cond_sets.txt"
)

# ---- (4) Save MOA gene sets ---- 

moa.targ.list = stack.long.clean %>%
  group_by(moa_targ_clue) %>%
  arrange(nchar(drug_name_raw), .by_group = TRUE) %>%
  summarise(
    # --- flat identifiers ---
    moa_targ_clue   = first(na.omit(moa_targ_clue)),
    
    # --- clinical stage ---
    max_clinical_stage   = if(all(is.na(max_clinical_stage))) NA_character_ else get.max.stage(na.omit(max_clinical_stage)),
    
    # --- list-columns ---
    harmonized_chembl_id = list(collapse_unique(harmonized_chembl_id)),
    drug_name            = list(collapse_unique(drug_name)),
    drug_name_raw        = list(collapse_unique(drug_name_raw)),
    master_id       = list(collapse_unique(master_id)),
    ensembl         = list(collapse_unique(ensembl)),
    entrez          = list(collapse_unique(entrez)),
    symbol          = list(collapse_unique(symbol)),
    moa_chembl      = list(collapse_unique(moa_chembl)),
    moa_targ_chembl = list(collapse_unique(moa_targ_chembl)),
    moa_clue        = list(collapse_unique(moa_clue)),
    moa_ot          = list(collapse_unique(moa_ot)),
    moa_targ_ot     = list(collapse_unique(moa_targ_ot)),
    atc1            = list(collapse_unique(atc1)),
    atc2            = list(collapse_unique(atc2)),
    atc3            = list(collapse_unique(atc3)),
    atc4            = list(collapse_unique(atc4)),
    indication_efo  = list(collapse_unique(indication_efo)),
    
    .groups = "drop"
  ) %>%
  mutate(
    moa_targ_clue_clean = str_replace_all(moa_targ_clue, " ", "-")
  ) %>%
  mutate(
    set.name = moa_targ_clue_clean
  ) %>%
  dplyr::select(set.name, everything())

# For the MoA data frame (one MoA per row)
write.magma.geneset(
  df = moa.targ.list, 
  set.id = "set.name", 
  gene.id = "ensembl", 
  file.path = "/Users/nyb/drugsets/DATA/GENESETS/moa_ensembl_cond_sets.txt"
)


# ---- (5) Save meta data for drug sets ---- 


#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

# Define the expected command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to input gene set file (e.g., .csv or .rds)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="filtered_sets.csv", 
              help="Path to output file [default= %default]", metavar="character"),
  make_option(c("-m", "--min_genes"), type="integer", default=10, 
              help="Minimum number of genes required to retain a drug set [default= %default]", metavar="integer"),
  make_option(c("-r", "--remove_hla"), type="logical", default=TRUE, 
              help="Boolean flag to remove HLA genes from sets [default= %default]", metavar="logical")
)

# Parse the arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Enforce required arguments
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Execution halted: --input argument must be supplied.", call.=FALSE)
}

# ---------------------------------------------------------
# Processing Logic
# ---------------------------------------------------------

# Read data
# Example assumes a long-format CSV with 'drug_id' and 'gene_symbol'
gene_data <- read.csv(opt$input)

# Apply dynamic filtering based on parsed arguments
if (opt$remove_hla) {
  gene_data <- gene_data %>%
    filter(!grepl("^HLA-", gene_symbol))
}

filtered_data <- gene_data %>%
  group_by(drug_id) %>%
  filter(n() >= opt$min_genes) %>%
  ungroup()

# Write output
write.csv(filtered_data, opt$output, row.names = FALSE)











