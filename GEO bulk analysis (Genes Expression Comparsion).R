#---------------------------
# 0. Environment Setup
#---------------------------
# Check and install pacman if missing
if (!require("pacman")) install.packages("pacman")

# Load required packages with pacman (auto-installs missing packages)
pacman::p_load(tidyverse, ggpubr, rstatix, purrr,fs)

#---------------------------
# 1. Data Preparation
#---------------------------
# Define core analysis genes
CORE_GENES <- c("NFE2L2", "FTH1", "GPX4", "KEAP1", 
               "NFKB1","SLC40A1", "SLC11A2")

# Read and process expression matrix
expr_matrix <- read_csv("expression_matrix.csv") %>%
  column_to_rownames("Symbol") %>%  # Set gene symbols as row names
  t() %>%  # Transpose to samples-as-rows format
  as.data.frame() %>%
  rownames_to_column("ID")  # Convert row names to ID column

# Read and process clinical data
clinical <- read_csv("clinical.csv") %>%
  mutate(Group = factor(Group, levels = c("Control", "Epilepsy")))

# Merge datasets and select relevant columns
merged_data <- expr_matrix %>%
  inner_join(clinical, by = "ID") %>%
  select(ID, Group, all_of(CORE_GENES))  # Use predefined gene list

#---------------------------
# 2. Data Reshaping
#---------------------------
long_data <- merged_data %>%
  pivot_longer(
    cols = -c(ID, Group),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(Gene = factor(Gene, levels = CORE_GENES))  # Maintain gene order

#---------------------------
# 3. Statistical Analysis (Raw P-value Version)
#---------------------------
stat_test <- long_data %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Group) %>%
  add_significance("p") %>%  # Use raw p-values
  add_xy_position(x = "Group", dodge = 0.8) %>%
  mutate(
    xmin = as.numeric(factor(group1, levels = levels(clinical$Group))),
    xmax = as.numeric(factor(group2, levels = levels(clinical$Group)))
  )

#---------------------------
# 4. Visualization Functions
#---------------------------
# Define consistent plot aesthetics
PLOT_THEME <- theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.subtitle = element_text(color = "gray40", size = 10)
  )

COLOR_SCHEME <- c(Control = "#1b9e77", Epilepsy = "#d95f02")

# Boxplot generator function
create_boxplot <- function(gene) {
  plot_data <- filter(long_data, Gene == gene)
  
  ggplot(plot_data, aes(x = Group, y = Expression)) +
    geom_boxplot(aes(fill = Group), width = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
    stat_pvalue_manual(
      filter(stat_test, Gene == gene),
      label = "p.signif",
      tip.length = 0.01,
      y.position = "y.position"
    ) +
    scale_fill_manual(values = COLOR_SCHEME) +
    labs(
      title = paste(gene, "Expression"),
      y = "TPM Expression (log2)",
      subtitle = "Raw p-value significance"
    ) +
    PLOT_THEME
}

#---------------------------
# 5. Generate Expression Plots
#---------------------------
# Create output directory if missing
dir_create("plots")

# Generate and save individual plots
walk(CORE_GENES, ~{
  plot_path <- path("plots", paste0("Expression_", .x, ".pdf"))
  
  ggsave(
    plot_path,
    plot = create_boxplot(.x),
    device = cairo_pdf,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300
  )
})

#---------------------------
# 6. Spearman Correlation Analysis
#---------------------------
# Define gene pairs for correlation analysis
GENE_PAIRS <- list(
  c("NFE2L2", "SLC40A1"),
  c("NFE2L2", "FTH1")
)

# Calculate correlations with FDR adjustment
cor_results <- map_df(GENE_PAIRS, ~{
  test_res <- cor.test(
    pull(merged_data, .x[1]),
    pull(merged_data, .x[2]),
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    Gene1 = .x[1],
    Gene2 = .x[2],
    rho = test_res$estimate,
    p.value = test_res$p.value
  )
}) %>%
  adjust_pvalue(p.col = "p.value", method = "BH") %>%
  add_significance("p.value.adj")

#---------------------------
# 7. Correlation Visualization
#---------------------------
# Correlation plot generator function
create_corplot <- function(gene_pair) {
  ggplot(merged_data, aes(x = .data[[gene_pair[1]]], y = .data[[gene_pair[2]]])) +
    geom_point(color = "#2c7bb6", size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "#d7191c", se = FALSE) +
    stat_cor(
      method = "spearman",
      label.sep = "\n",
      aes(label = paste(..r.label.., ..p.label.., sep = "~`; `~"))
    ) +
    labs(
      title = paste(gene_pair[1], "vs", gene_pair[2]),
      x = paste(gene_pair[1], "Expression (TPM)"),
      y = paste(gene_pair[2], "Expression (TPM)")
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

# Generate and save correlation plots
walk(GENE_PAIRS, ~{
  plot_path <- path("plots", paste0("Correlation_", .x[1], "_vs_", .x[2], ".pdf"))
  
  ggsave(
    plot_path,
    plot = create_corplot(.x),
    device = cairo_pdf,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300
  )
})
