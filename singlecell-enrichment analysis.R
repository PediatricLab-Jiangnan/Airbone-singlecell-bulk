# Load required packages
library(Seurat)          
library(tidyverse)        
library(EnhancedVolcano)  
library(clusterProfiler)  
library(org.Hs.eg.db)    
library(enrichplot)      
library(RColorBrewer)    

#--------------------------
# Differential Expression Analysis
#--------------------------

#' Perform differential expression analysis between conditions
#'
#' @param seurat_obj A Seurat object containing single-cell data
#' @param cond1 Name of experimental condition (e.g., "Seizure")
#' @param cond2 Name of control condition (e.g., "Control")
#' @return Data frame of DE results with gene annotations

run_differential_analysis <- function(seurat_obj, cond1, cond2) {
  FindMarkers(
    seurat_obj,
    ident.1 = cond1,
    ident.2 = cond2,
    group.by = "condition",
    min.pct = 0.1,          # Require genes detected in at least 10% of cells
    logfc.threshold = 0.25, # Minimum log2 fold change for preliminary filtering
    test.use = "wilcox"     # Non-parametric Wilcoxon rank sum test
  ) %>%
    rownames_to_column("gene") %>%
    mutate(
      significance = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5,
      direction = case_when(
        avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Up",
        avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Down",
        TRUE ~ "NS"
      )
    )
}

# Execute DE analysis
deg_results <- run_differential_analysis(merged_seurat, "Seizure", "Control")

#--------------------------
# Visualization: Volcano Plot
#--------------------------

# Configure color scheme
volcano_colors <- c(NS = "grey80", Up = "red3", Down = "blue4")

# Create enhanced volcano plot
volcano_plot <- EnhancedVolcano(
  deg_results,
  lab = deg_results$gene,
  x = "avg_log2FC",
  y = "p_val_adj",
  title = "Seizure vs Control",
  subtitle = "Differential Expression Analysis",
  pCutoff = 0.05,                   # Adjusted p-value cutoff
  FCcutoff = 0.5,                   # Log2 fold change cutoff
  colCustom = unname(volcano_colors[deg_results$direction]),
  colAlpha = 0.7,                   # 70% opacity for points
  pointSize = 2.5,
  labSize = 3.5,
  legendPosition = "bottom",
  drawConnectors = TRUE,            # Connect labels to points
  max.overlaps = 30,                # Max label overlaps allowed
  border = "full",                  # Add plot border
  borderWidth = 1.5,
  borderColour = "black"
)

# Save visualization
ggsave("volcano_plot.pdf", volcano_plot, width = 10, height = 8)

#--------------------------
# Result Summary and Export
#--------------------------

# Generate summary statistics
deg_summary <- list(
  total_deg = sum(deg_results$significance),
  up_regulated = sum(deg_results$direction == "Up"),
  down_regulated = sum(deg_results$direction == "Down"),
  top_genes = deg_results %>%
    filter(significance) %>%
    arrange(p_val_adj) %>%
    head(10) %>%
    select(gene, avg_log2FC, p_val_adj)
)

# Print summary to console
cat(str_glue(
  "Differential Expression Summary:
  Total DEGs: {deg_summary$total_deg}
  Upregulated: {deg_summary$up_regulated}
  Downregulated: {deg_summary$down_regulated}\n"
))

# Export full results
write_csv(deg_results, "full_de_results.csv")

# Export significant results
deg_results %>%
  filter(significance) %>%
  write_csv("significant_de_genes.csv")

#--------------------------
# Functional Enrichment Analysis
#--------------------------

#' Perform functional enrichment analysis
#'
#' @param deg_df Data frame containing DE results
#' @param p_cutoff Adjusted p-value cutoff
#' @param fc_cutoff Absolute log2 fold change cutoff
#' @return List containing GO and KEGG enrichment results

run_enrichment_analysis <- function(deg_df, p_cutoff = 0.05, fc_cutoff = 0.5) {
  # Convert gene symbols to Entrez IDs
  gene_symbols <- deg_df %>%
    filter(p_val_adj < p_cutoff & abs(avg_log2FC) > fc_cutoff) %>%
    pull(gene)
  
  gene_conversion <- bitr(gene_symbols, 
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
  
  # GO enrichment analysis
  go_results <- enrichGO(
    gene = gene_conversion$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",                  # Biological Process
    pAdjustMethod = "BH",        # Benjamini-Hochberg correction
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # KEGG pathway analysis
  kegg_results <- enrichKEGG(
    gene = gene_conversion$ENTREZID,
    organism = "hsa",            # Human KEGG organism code
    pvalueCutoff = 0.05,
    keyType = "ncbi-geneid"
  )
  
  list(go = go_results, kegg = kegg_results)
}

# Execute enrichment analysis
enrichment_results <- run_enrichment_analysis(deg_results)

#--------------------------
# Enrichment Visualization
#--------------------------

# Custom theme for enrichment plots
enrichment_theme <- theme(
  plot.title = element_text(face = "bold", size = 14),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  legend.position = "right",
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(color = "grey90")
)

# Generate GO enrichment plot
go_plot <- dotplot(enrichment_results$go, 
                   showCategory = 20,
                   font.size = 10) +
  ggtitle("GO Biological Processes") +
  enrichment_theme +
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "Spectral")),
    name = "-log10(p.adjust)"
  )

# Generate KEGG pathway plot
kegg_plot <- dotplot(enrichment_results$kegg, 
                     showCategory = 20,
                     font.size = 10) +
  ggtitle("KEGG Pathways") +
  enrichment_theme +
  scale_color_distiller(
    palette = "RdYlBu",
    direction = -1,
    name = "-log10(p.adjust)"
  )

# Save enrichment plots
ggsave("go_enrichment.pdf", go_plot, width = 12, height = 9)
ggsave("kegg_enrichment.pdf", kegg_plot, width = 12, height = 9)
