#---------------------------
# 0. Environment Setup
#---------------------------

# Install required packages with version checking
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_bioc <- c("GEOquery", "tidyverse")
BiocManager::install(required_bioc, update = FALSE)
suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

#---------------------------
# 1. Data Loading
#---------------------------
# Load expression data from series matrix
gse <- getGEO(filename = "GSEXXXXXX", #You own GSE file
              GSEMatrix = TRUE,
              getGPL = FALSE)  

#---------------------------
# 2. Expression Matrix Extraction
#---------------------------
expr_matrix <- exprs(gse) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Probe_ID")

#---------------------------
# 3. Annotation Processing
#---------------------------
#' Custom Platform Annotation Processor
#' @param gpl_number GEO platform number (e.g., "GPL16791")
process_platform_annotation <- function(gpl_number) {
  # Download platform data with cache
  gpl <- getGEO(gpl_number, destdir = ".", AnnotGPL = TRUE)
  
  # Extract feature data with key columns
  feature_data <- Table(gpl) %>% 
    dplyr::select(ID, Symbol, Entrez_Gene_ID, Chromosome) %>% 
    dplyr::distinct()  # Remove duplicate probes
  
  # Validate ID matching
  validate_id_match(expr_matrix$Probe_ID, feature_data$ID)
  
  return(feature_data)
}

# Execute annotation processing
annot_df <- process_platform_annotation(gse@annotation)

#---------------------------
# 4. Annotation Merging
#---------------------------
annotated_data <- expr_matrix %>% 
  inner_join(annot_df, by = c("Probe_ID" = "ID")) %>%
  dplyr::filter(!is.na(Symbol) & Symbol != "") %>%  # Remove empty symbols
  dplyr::group_by(Symbol) %>% 
  dplyr::summarise(across(where(is.numeric), mean), .groups = 'drop') %>%  # Average duplicates
  dplyr::rename_with(~gsub("GSM\\d+", "Sample_", .x))  # Standardize sample names
#---------------------------
# 5. Low Expression Filtering
#---------------------------
annotated_data_filtered <- annotated_data %>%
  dplyr::mutate(gene_median = matrixStats::rowMedians(as.matrix(select(., starts_with("Sample_"))))) %>%
  dplyr::filter(gene_median > quantile(gene_median, probs = 0.25)) %>%
  dplyr::select(-gene_median)
# Save annotated matrix

write_csv(annotated_data, "expression_matrix.csv")
