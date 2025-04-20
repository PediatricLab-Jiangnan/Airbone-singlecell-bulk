# Set working directory
setwd("Your own path")

# Load required libraries
library(dplyr)             # Data manipulation
library(corrr)             # Correlation calculation
library(knitr)             # Table formatting (not used here, can remove if unnecessary)
library(ComplexHeatmap)    # Advanced heatmap plotting
library(circlize)          # Color mapping utilities
library(RColorBrewer)      # Color palettes

# --- Data Preparation ---
# Read the dataset
file_path <- "combined.csv"
df <- read.csv(file_path)

# Select only numeric columns, excluding 'Group'
numeric_cols <- df %>%
  select(-Group) %>%
  select(where(is.numeric))

# Calculate Spearman correlation matrix, handling missing values
cor_matrix <- cor(
  numeric_cols,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# --- Color Mapping Optimization ---
# Create a color function for the heatmap using RdBu palette
col_fun <- colorRamp2(
  seq(-1, 1, length.out = 11),
  rev(brewer.pal(11, "RdBu"))
)

# --- Draw Heatmap ---
# Add a star '*' for correlation coefficients with absolute value > 0.9 (excluding diagonal)
ht <- Heatmap(
  cor_matrix,
  name = "Spearman Corr",
  col = col_fun,
  rect_gp = gpar(col = "grey90", lwd = 0.2),          # Set cell border color and width
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_title = "Spearman Correlation Heatmap",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = gpar(fontsize = 8, fontface = "plain"),
  column_names_gp = gpar(fontsize = 8, fontface = "plain"),
  heatmap_legend_param = list(
    title = "Spearman Corr",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    legend_height = unit(3, "cm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Add a star for highly correlated pairs (excluding the diagonal)
    if (i != j && abs(cor_matrix[i, j]) > 0.9) {
      grid.text("*", x, y, gp = gpar(fontsize = 18, col = "black", fontface = "bold"))
    }
  }
)

# Draw the heatmap with legends on the right and annotation legends at the bottom
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "bottom")
