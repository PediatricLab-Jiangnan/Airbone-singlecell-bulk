###Create exposure.txt and added your own interest GWAS (exposure)####
# Set working directory (modify path as needed)
setwd("Your own path")

# Load required packages
library(TwoSampleMR)  
library(dplyr)       
library(readr)      
library(ggplot2)   
library(tidyr)       

# Custom function for error logging
log_message <- function(message, type = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0("[", timestamp, "] ", type, ": ", message, "\n"))
}

# Initialize empty list to store results 
results_list <- list()

# Create output directories with error handling
create_dir <- function(dir_name) {
  if (!dir.exists(dir_name)) {
    dir_result <- tryCatch(
      {
        dir.create(dir_name)
        log_message(paste("Directory created:", dir_name))
      },
      error = function(e) {
        log_message(paste("Failed to create directory:", dir_name, "| Error:", e$message), "ERROR")
        stop("Fatal directory creation error")  
      }
    )
  }
}

create_dir("Scatter_Plots")
create_dir("Harmonized_Data")

# Read exposure IDs with validation
if (!file.exists("exposure.txt")) {
  stop("Exposure ID file 'exposure.txt' not found in working directory")
}
exposure_ids <- read_lines("exposure.txt") %>% 
  na.omit() %>% 
  unique()

# Parameter configuration 
PARAMS <- list(
  p1_threshold = 1e-05,
  clump_r2 = 0.001,
  clump_kb = 10000,
  maf_threshold = 0.01,
  outcome_id = "Your own outcome gwas id"
)

# Main processing loop
for (i in seq_along(exposure_ids)) {
  exposure_id <- exposure_ids[i]
  log_message(paste("Processing", i, "/", length(exposure_ids), "| Exposure ID:", exposure_id))
  
  result <- tryCatch({
    # ----------------------------------------------------------
    # Step 1: Extract exposure instruments
    # ----------------------------------------------------------
    exposure_data <- extract_instruments(
      outcomes = exposure_id,
      p1 = PARAMS$p1_threshold,
      clump = TRUE,
      r2 = PARAMS$clump_r2,
      kb = PARAMS$clump_kb
    )
    
    # Validate exposure data
    if (is.null(exposure_data) || nrow(exposure_data) == 0) {
      msg <- paste("No valid SNPs found for", exposure_id)
      log_message(msg, "WARN")
      return(list(success = FALSE, message = msg))
    }
    
    # ----------------------------------------------------------
    # Step 2: Extract outcome data
    # ----------------------------------------------------------
    outcome_data <- extract_outcome_data(
      snps = exposure_data$SNP,
      outcomes = PARAMS$outcome_id,
      proxies = FALSE,
      maf_threshold = PARAMS$maf_threshold
    )
    
    # Validate outcome data
    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      msg <- paste("No outcome data for", exposure_id)
      log_message(msg, "WARN")
      return(list(success = FALSE, message = msg))
    }
    
    # ----------------------------------------------------------
    # Step 3: Harmonize data
    # ----------------------------------------------------------
    harmonized_data <- harmonise_data(exposure_data, outcome_data)
    
    if (nrow(harmonized_data) == 0) {
      msg <- paste("Harmonization failed for", exposure_id)
      log_message(msg, "WARN")
      return(list(success = FALSE, message = msg))
    }
    
    # ----------------------------------------------------------
    # Step 4: Perform MR analysis
    # ----------------------------------------------------------
    mr_results <- mr(harmonized_data)
    
    if (nrow(mr_results) == 0) {
      msg <- paste("MR analysis failed for", exposure_id)
      log_message(msg, "WARN")
      return(list(success = FALSE, message = msg))
    }
    
    # ----------------------------------------------------------
    # Step 5: Calculate ORs with 95% CI
    # ----------------------------------------------------------
    mr_results <- mr_results %>%
      mutate(
        OR = exp(b),
        OR_lower = exp(b - 1.96 * se),
        OR_upper = exp(b + 1.96 * se),
        exposure_id = exposure_id,
        analysis_timestamp = Sys.time()
      ) %>%
      relocate(exposure_id, .before = method)
    
    # ----------------------------------------------------------
    # Step 6: Save harmonized data
    # ----------------------------------------------------------
    safe_write <- function(data, path) {
      tryCatch(
        {
          write_csv(data, path)
          log_message(paste("File saved:", path))
        },
        error = function(e) {
          log_message(paste("Failed to save", path, "| Error:", e$message), "ERROR")
        }
      )
    }
    
    safe_write(
      harmonized_data,
      file.path("Harmonized_Data", paste0("harmonized_", exposure_id, ".csv"))
    )
    
    # ----------------------------------------------------------
    # Step 7: Generate and save scatter plot
    # ----------------------------------------------------------
    tryCatch(
      {
        p <- mr_scatter_plot(mr_results, harmonized_data) +
          labs(
            title = paste("Exposure ID:", exposure_id),
            subtitle = paste("Outcome:", PARAMS$outcome_id),
            caption = paste(nrow(harmonized_data), "SNPs |", Sys.Date())
          )
        
        ggsave(
          filename = file.path("Scatter_Plots", paste0("scatter_", exposure_id, ".pdf")),
          plot = p,
          width = 8,
          height = 6,
          device = "pdf"
        )
      },
      error = function(e) {
        log_message(paste("Plot generation failed for", exposure_id, "| Error:", e$message), "ERROR")
      }
    )
    
    # ----------------------------------------------------------
    # Step 8: Return success status
    # ----------------------------------------------------------
    list(success = TRUE, data = mr_results)
    
  }, error = function(e) {
    log_message(paste("Critical error in", exposure_id, "|", e$message), "ERROR")
    return(list(success = FALSE, message = e$message))
  })
  
  # Store results if successful
  if (result$success) {
    results_list[[exposure_id]] <- result$data
  }
  
  # Garbage collection for large datasets
  gc()
}

# ----------------------------------------------------------
# Post-processing
# ----------------------------------------------------------

# Combine results
final_results <- bind_rows(results_list) %>%
  select(
    exposure_id, method, nsnp,
    b, se, pval,
    OR, OR_lower, OR_upper,
    analysis_timestamp
  )

# Save final results
safe_write(final_results, "MR_results_compiled.csv")

# Session info logging
sink("analysis_session_info.txt")
print(sessionInfo())
sink()

log_message("Analysis pipeline completed successfully")
