# Clear environment and load packages
rm(list = ls())  
library(vroom)    
library(tidyr)   
library(dplyr)   
library(data.table) 

# Set working directory (modify path accordingly)
setwd("Yueying Liu Lab")

# 1. Data Loading & Validation  download from https://www.ebi.ac.uk/gwas/ GWAS Caltalog
if(!file.exists('Your own GWAS FILE')) {
  stop("Input file not found! Please check file path.")
}

# Read compressed GWAS data with progress
Ruirui <- vroom('Your own GWAS FILE', 
                col_names = TRUE,
                show_col_types = FALSE,
                progress = TRUE)

# 2. Data Transformation Pipeline
processed_data <- Ruirui %>%
  # Column renaming and type conversion
  rename(
    SNP = variant_id,
    chr = chromosome,
    pos = base_pair_location,
    effect_allele = effect_allele,
    other_allele = other_allele,
    pval = p_value,
    eaf = effect_allele_frequency,
    beta = beta,
    se = standard_error
  ) %>%
  # Essential columns for SMR analysis
  select(chr, pos, other_allele, effect_allele, SNP, pval, eaf, beta, se) %>%
  mutate(
    chr = as.integer(chr),  # Ensure numeric chromosome
    pos = as.integer(pos),  # Ensure numeric position
    pval = as.numeric(pval), # Numeric p-values
    across(c(eaf, beta, se), as.numeric) # Ensure numeric types
  ) %>%
  # Filter valid genetic associations
  filter(pval > 0 & pval <= 1 & eaf > 0 & eaf < 1)

# 3. Save processed data
# Standard CSV format
data.table::fwrite(processed_data, "epilepsy.csv", row.names = FALSE)

