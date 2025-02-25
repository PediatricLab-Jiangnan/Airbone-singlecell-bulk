# Load required packages  
library(Boruta)  

# Set working directory  
setwd("Your Path here")  

# Load dataset and handle missing values  
data <- read.csv("Your Clinical information")  
TS_Syndrome <- na.omit(data)  

# Ensure the dependent variable is a factor  
TS_Syndrome$Group <- as.factor(TS_Syndrome$Group)  

# Set random seed for reproducibility  
set.seed(1234567)   

# Run Boruta feature selection algorithm  
Boruta.TS <- Boruta(Group ~ ., data = TS_Syndrome, doTrace = 2, maxRuns = 500)  

# Save plot as PDF with larger dimensions  
pdf("Boruta_TS_plot.pdf", width = 15, height = 10)  
# Adjust margins to accommodate long variable names  
par(mar = c(12, 4, 4, 2))  
# Create Boruta algorithm result plot  
plot(Boruta.TS,   
     las = 2,            # Vertical x-axis labels  
     cex.axis = 0.8,     # Adjust axis label size  
     xlab = "",          # Remove x-axis title  
     ylab = "Importance",  
     main = "Variable Importance in TS Syndrome Classification",  
     cex.main = 1.2)     # Adjust main title size  
grid()                   # Add grid lines  
dev.off()  

# Get feature selection results  
boruta_signif <- attStats(Boruta.TS)  

# Save results to CSV file  
write.csv(boruta_signif, "TS_Syndrome_Boruta_results.csv")  

# Results output section  
cat("\nTic Syndrome Feature Selection Results:\n\n")  

cat("Confirmed Important Variables:\n")  
confirmed_vars <- rownames(boruta_signif)[boruta_signif$decision == "Confirmed"]  
print(confirmed_vars)  

cat("\nTentative Important Variables:\n")  
tentative_vars <- rownames(boruta_signif)[boruta_signif$decision == "Tentative"]  
print(tentative_vars)  

cat("\nRejected Variables:\n")  
rejected_vars <- rownames(boruta_signif)[boruta_signif$decision == "Rejected"]  
print(rejected_vars)  

cat("\nFeature Selection Summary:\n")  
cat("Total Variables:", nrow(boruta_signif), "\n")  
cat("Confirmed Important Variables:", sum(boruta_signif$decision == "Confirmed"), "\n")  
cat("Tentative Variables:", sum(boruta_signif$decision == "Tentative"), "\n")  
cat("Rejected Variables:", sum(boruta_signif$decision == "Rejected"), "\n")  

# Output detailed statistics for important variables  
cat("\nImportant Variables Statistics (Sorted by Importance):\n")  
important_vars <- subset(boruta_signif, decision == "Confirmed")  
important_vars <- important_vars[order(-important_vars$meanImp), ]  
print(important_vars)  

# Add calculation of mean and standard deviation for important variables  
if(length(confirmed_vars) > 0) {  
    cat("\nDescriptive Statistics for Important Variables:\n")  
    important_stats <- data.frame(  
        Mean_Importance = important_vars$meanImp,  
        SD_Importance = important_vars$sdImp,  
        row.names = rownames(important_vars)  
    )  
    print(important_stats)  
}  
