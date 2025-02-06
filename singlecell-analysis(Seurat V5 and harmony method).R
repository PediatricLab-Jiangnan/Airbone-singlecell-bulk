# Load required packages ------------------------------------------------------
library(Seurat)         
library(harmony)        
library(ggplot2)        
library(dplyr)          
library(patchwork)      
library(clustree)       

# Configuration parameters ----------------------------------------------------
config <- list(
  # Directory settings
  work_dir = "path/to/working/directory",
  
  # QC parameters
  qc = list(
    min_features = 200,     # Minimum genes per cell
    max_features = 6000,    # Maximum genes per cell
    max_mt = 20,            # Maximum mitochondrial percentage
    mito_pattern = "^MT-"   # Mitochondrial gene pattern (adjust for organism)
  ),
  
  # Processing parameters
  processing = list(
    n_var_features = 2000,  # Number of variable features
    npcs = 30,              # Principal components for PCA
    clust_dims = 1:20,      # Dimensions for clustering
    resolutions = seq(0.1, 1.0, 0.1) # Clustering resolutions to test
  ),
  
  # Visualization parameters
  visualization = list(
    theme_settings = theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold")),
    color_conditions = c("Control" = "#E41A1C", "Seizure" = "#377EB8")
  )
)

# Set working directory -------------------------------------------------------
setwd(config$work_dir)


create_seurat_object <- function(matrix_obj, sample_name, condition) {
  seurat_obj <- CreateSeuratObject(
    counts = matrix_obj,
    project = sample_name,
    min.cells = 3,    # Filter genes detected in <3 cells
    min.features = 200 # Filter cells with <200 genes
  )
  
  # Add metadata
  seurat_obj$batch <- sample_name
  seurat_obj$condition <- condition
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj, 
    pattern = config$qc$mito_pattern
  )
  
  return(seurat_obj)
}

#' Generate QC violin plots
#' @param seurat_obj Seurat object
#' @param file_name Output file name
generate_qc_plots <- function(seurat_obj, file_name) {
  qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  plots <- VlnPlot(
    seurat_obj,
    features = qc_features,
    group.by = "batch",
    pt.size = 0,
    ncol = 3
  )
  
  # Customize plot titles
  plots <- lapply(1:3, function(i) {
    plots[[i]] + ggtitle(qc_features[i]) + config$visualization$theme_settings
  })
  
  # Save combined plot
  pdf(file_name, width = 15, height = 6)
  print(wrap_plots(plots))
  dev.off()
}

# Data loading and preprocessing ----------------------------------------------
samples <- list()

for (i in 1:3) {
  # Process control samples
  con_name <- paste0("Con", i)
  samples[[con_name]] <- create_seurat_object(
    readRDS(paste0(con_name, ".rds")),
    con_name,
    "Control"
  )
  
  # Process seizure samples
  se_name <- paste0("Se", i)
  samples[[se_name]] <- create_seurat_object(
    readRDS(paste0(se_name, ".rds")),
    se_name,
    "Seizure"
  )
}

# Pre-QC visualization --------------------------------------------------------
merged_preqc <- merge(samples[[1]], y = samples[-1], add.cell.ids = names(samples))
generate_qc_plots(merged_preqc, "pre_qc_violins_combined.pdf")

# Quality control -------------------------------------------------------------
samples <- lapply(samples, function(obj) {
  subset(obj, subset = 
           nFeature_RNA > config$qc$min_features &
           nFeature_RNA < config$qc$max_features &
           percent.mt < config$qc$max_mt)
})

# Post-QC processing ----------------------------------------------------------
merged_seurat <- merge(samples[[1]], y = samples[-1], add.cell.ids = names(samples))
generate_qc_plots(merged_seurat, "post_qc_violins_combined.pdf")

# Clean memory
rm(samples, merged_preqc)
gc()

# Data normalization and integration ------------------------------------------
merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = config$processing$n_var_features) %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npcs = config$processing$npcs) %>%
  IntegrateLayers(
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE
  )

# Save intermediate results
saveRDS(merged_seurat, "merged_seurat_harmony.rds")

# Dimensionality assessment ---------------------------------------------------
pdf("pca_elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(merged_seurat, ndims = config$processing$npcs) + 
  ggtitle("PCA Elbow Plot") +
  config$visualization$theme_settings
dev.off()

# Clustering analysis ---------------------------------------------------------
merged_seurat <- merged_seurat %>%
  FindNeighbors(reduction = "harmony", dims = config$processing$clust_dims) %>%
  FindClusters(resolution = config$processing$resolutions) %>%
  RunUMAP(reduction = "harmony", dims = config$processing$clust_dims)

# Clustering resolution evaluation --------------------------------------------
pdf("clustree_analysis.pdf", width = 12, height = 16)
clustree(merged_seurat, prefix = "RNA_snn_res.") + 
  ggtitle("Cluster Resolution Evaluation") +
  theme(legend.position = "right")
dev.off()

# Generate UMAP plots for all resolutions -------------------------------------
plot_list <- lapply(config$processing$resolutions, function(res) {
  DimPlot(merged_seurat,
          reduction = "umap",
          group.by = paste0("RNA_snn_res.", res),
          label = TRUE) +
    ggtitle(paste("Resolution:", res)) +
    config$visualization$theme_settings
})

pdf("umap_all_resolutions.pdf", width = 20, height = 15)
wrap_plots(plot_list, ncol = 3)
dev.off()

# Cell annotation -------------------------------------------------------------
Idents(merged_seurat) <- "RNA_snn_res.0.2"

# Marker gene identification
markers <- FindAllMarkers(
  merged_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Save marker results
write.csv(markers, "all_markers_res0.2.csv", row.names = FALSE)

# Cluster annotation ####### Annoation using scGPT (the https://github.com/bowang-lab/scGPT)
cluster_ids <- c(
  "Neurons", "Fibroblasts", "Oligodendrocytes", "Astrocytes",
  "OPCs", "Interneurons", "GABAergic Interneurons", 
  "Endothelial Cells", "Endothelial Progenitors", "Microglia",
  "Microglia", "Schwann Cells", "Pericytes", "T cells",
  "Neurons", "Fibroblasts", "Fibroblasts"
)

# Validate cluster annotations
if(length(levels(merged_seurat)) != length(cluster_ids)) {
  stop("Cluster count mismatch! Verify annotation labels.")
}

names(cluster_ids) <- levels(merged_seurat)
merged_seurat <- RenameIdents(merged_seurat, cluster_ids)
merged_seurat$cell_type <- Idents(merged_seurat)

# Marker gene visualization ---------------------------------------------------
marker_genes <- list(
  "Neurons" = c("NEFM", "NEFL", "NPSR1"),
  "Fibroblasts" = c("ADAMTS16", "FAP", "COL15A1"),
  "Oligodendrocytes" = c("MOBP", "MBP"),
  "Astrocytes" = c("AQP4", "GFAP"),
  "OPCs" = c("PDGFRA", "CSPG4"),
  "Microglia" = c("CX3CR1", "P2RY12")
)

# Generate diagnostic plots
pdf("celltype_validation_plots.pdf", width = 15, height = 10)
DotPlot(merged_seurat, features = marker_genes) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(merged_seurat, 
            features = unlist(marker_genes)[1:6],
            ncol = 3,
            order = TRUE)
dev.off()

# Final visualizations --------------------------------------------------------
plot_umap <- function(group_var, title) {
  DimPlot(merged_seurat,
          reduction = "umap",
          group.by = group_var,
          cols = if(group_var == "condition") config$visualization$color_conditions,
          label = TRUE) +
    ggtitle(title) +
    config$visualization$theme_settings
}

pdf("final_umaps.pdf", width = 12, height = 10)
plot_umap("cell_type", "Cell Type Annotation")
plot_umap("condition", "Experimental Condition")
plot_umap("batch", "Batch Distribution")
dev.off()

# Save final annotated object -------------------------------------------------
saveRDS(merged_seurat, "annotated_seurat_object.rds")
