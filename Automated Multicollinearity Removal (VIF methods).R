# Load required packages  
library(car)  
library(ggplot2)  

# Set working directory  
setwd("Your own path")  

# Load dataset and handle missing values  
data <- read.csv("your own information")  
TS_Syndrome <- na.omit(data)  # Rename to TS_Syndrome  

# Create linear model and calculate VIF values  
model <- lm(Group ~ ., data = TS_Syndrome)  
vif_values <- vif(model)  

# Create dataframe and filter VIF < 10  
vif_df <- data.frame(  
    Variable = names(vif_values),  
    VIF = vif_values  
) %>%   
    subset(VIF < 10) %>%  
    arrange(desc(VIF))  # Sort by VIF in descending order  

# Create color gradient  
n_colors <- nrow(vif_df)  
colors <- colorRampPalette(c("#7CCD7C", "#6495ED", "#15317E"))(n_colors)  

# Create plot using ggplot2  
p <- ggplot(vif_df, aes(x = reorder(Variable, VIF), y = VIF)) +  
    geom_bar(stat = "identity", fill = colors) +  
    coord_flip() +  # Horizontal bars  
    theme_minimal() +  
    labs(  
        title = "Variance Inflation Factor (VIF) of Features",  
        x = "",  
        y = "VIF"  
    ) +  
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),  
        axis.text.y = element_text(size = 10),  
        axis.text.x = element_text(size = 10),  
        panel.grid.major.y = element_blank(),  
        panel.grid.minor = element_blank()  
    )  

# Save plots  
ggsave("VIF_plot_filtered.pdf", p, width = 10, height = 6)  
ggsave("VIF_plot_filtered.png", p, width = 10, height = 6, dpi = 300)  

# Print results  
print("Variables with VIF < 10:")  
print(vif_df)  

# Optional: Save detailed VIF results to CSV  
write.csv(vif_df, "VIF_results_filtered.csv", row.names = FALSE)
