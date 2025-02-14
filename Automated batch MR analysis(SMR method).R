# -------------------------------
# SMR Data Processing Pipeline
# Please following the https://yanglab.westlake.edu.cn/software/smr/#Overview to get SMR file before using this R code
# -------------------------------

### SECTION 1: DATA PROCESSING ###

# 1.1 Import SMR results with data validation
input_file <- "You own SMR file"

if (!file.exists(input_file)) {
  stop("Input SMR file not found. Check file path: ", input_file)
}

data <- read.delim(input_file, 
                   header = TRUE, 
                   stringsAsFactors = FALSE,
                   na.strings = c("NA", "", "NaN"))

# 1.2 Multiple testing correction
FDR_threshold <- 0.05  # Standard threshold for FDR
results <- data %>%
  mutate(
    FDR = p.adjust(p_SMR, method = "BH"),
    Significant = ifelse(FDR < FDR_threshold, "Yes", "No")
  )

# 1.3 Gene ID conversion with error handling
library(clusterProfiler)
library(org.Hs.eg.db)

convert_gene_ids <- function(ensembl_ids) {
  tryCatch({
    df <- bitr(ensembl_ids,
               fromType = 'ENSEMBL',
               toType = 'SYMBOL',
               OrgDb = org.Hs.eg.db)
    
    # Handle unmapped IDs
    missing_ids <- setdiff(ensembl_ids, df$ENSEMBL)
    if (length(missing_ids) > 0) {
      warning(paste(length(missing_ids), "ENSEMBL IDs failed to map:", 
                    paste(head(missing_ids), collapse = ", ")))
    }
    
    return(df)
  }, error = function(e) {
    stop("Gene ID conversion failed: ", e$message)
  })
}

gene_mapping <- convert_gene_ids(results$Gene)
results <- results %>%
  left_join(gene_mapping, by = c("Gene" = "ENSEMBL")) %>%
  rename(Gene_Symbol = SYMBOL) %>%
  relocate(Gene_Symbol, .after = Gene)

### SECTION 2: OUTPUT GENERATION ###

# 2.1 Save results with multiple format options
output_basename <- "SMR_results"

write.csv(results, 
          file = paste0(output_basename, ".csv"), 
          row.names = FALSE, 
          na = "")

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "SMR Results")
writeData(wb, 1, results)
saveWorkbook(wb, paste0(output_basename, ".xlsx"), overwrite = TRUE)

### SECTION 3: VISUALIZATION ###

# 3.1 Enhanced plotting functions
plot_smr_results <- function(plot_file, output_name = "SMR_Plot") {
  # 3.1.1 Data import validation
  if (!file.exists(plot_file)) {
    stop("Plot data file not found: ", plot_file)
  }
  
  # 3.1.2 Data processing
  smrdata <- tryCatch({
    ReadSMRData(plot_file)
  }, error = function(e) {
    stop("Failed to read SMR plot data: ", e$message)
  })
  
  # 3.1.3 Locus plot with custom parameters
  locus_plot <- SMRLocusPlot(
    data = smrdata,
    smr_thresh = 8.4e-6,     
    heidi_thresh = 0.05,     
    plotWindow = 1000,       
    title = "SMR Locus Plot",
    point_size = 2,
    gene_label_size = 3
  ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  
  # 3.1.4 Effect plot with axis transformation
  effect_plot <- SMREffectPlot(
    data = smrdata,
    trait_name = "Skin Not Sun Exposed",
    xlab = "eQTL Effect (log scale)",
    ylab = "GWAS Effect"
  ) +
    scale_x_continuous(trans = "log10") +
    theme_bw(base_size = 12)
  
  # 3.1.5 Save plots
  pdf(paste0(output_name, ".pdf"), width = 8, height = 6)
  print(locus_plot)
  print(effect_plot)
  dev.off()
  
  cat("Saved plots to:", paste0(output_name, ".pdf"), "\n")
}

# Example usage:
# plot_smr_results("../plot/myplot.ENSG00000235821.txt", "Final_SMR_Results")

