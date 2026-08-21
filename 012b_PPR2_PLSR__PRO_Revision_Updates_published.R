# 04.08.2026
# Written by Nargiz Rüter
# 
# ============================================================
# PPR2 REVISION — WITHIN-YEAR METRICS FROM SAVED .RData FILES
# UPDATED VERSION
# ============================================================
#
# This script:
#   1. Searches all immediate folders beginning with "PPR2_"
#   2. Finds files ending with "_plot_data.RData"
#   3. Loads ensemble calibration and validation predictions
#   4. Calculates reviewer-requested model diagnostics
#   5. Overwrites one revision-metrics CSV in each model folder
#   6. Overwrites one combined master CSV in the root folder
#   7. Saves a validation-only summary CSV in the root folder
#
# ============================================================

# ------------------------------------------------------------
# LOAD LIBRARIES
# ------------------------------------------------------------

library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tibble)

# ------------------------------------------------------------
# ROOT DIRECTORY
# ------------------------------------------------------------

root <- paste0(
  "C:/Users/..."
)

setwd(root)
options(scipen = 999)

# ------------------------------------------------------------
# LIST CANDIDATE FOLDERS
# ------------------------------------------------------------

folders <- list.dirs(root, full.names = TRUE, recursive = FALSE)
folders <- folders[grepl("^PPR2_", basename(folders))]

folder_info <- tibble(
  folder = folders,
  folder_name = basename(folders)
) %>%
  mutate(
    Calibration_Years = str_match(folder_name, "^PPR2_([^_]+)_")[,2],
    Run_Date = as.Date(str_match(folder_name, "^PPR2_[^_]+_(\\d{4}-\\d{2}-\\d{2})_")[,2]),
    Trait = str_match(folder_name, "^PPR2_[^_]+_\\d{4}-\\d{2}-\\d{2}_(.+)$")[,2]
  ) %>%
  filter(!is.na(Calibration_Years), !is.na(Run_Date), !is.na(Trait)) %>%
  group_by(Calibration_Years, Trait) %>%
  slice_max(Run_Date, n = 1, with_ties = FALSE) %>%
  ungroup()

folders <- folder_info$folder

cat("Newest PPR2 folders retained:", length(folders), "\n")

print(
  folder_info %>% select(Calibration_Years, Trait, Run_Date, folder_name),
  n = Inf
)

extract_valid_iterations <- function(folder) {
  files <- list.files(
    folder,
    pattern = "_REF_DF_MODELs-[0-9]{1,4}-STATISTICS\\.csv$",
    full.names = FALSE
  )
  
  if (length(files) == 0) return(NA_integer_)
  
  as.integer(
    str_match(files[1], "_REF_DF_MODELs-([0-9]{1,4})-STATISTICS\\.csv$")[1,2]
  )
}

# ------------------------------------------------------------
# HELPER FUNCTION: POPULATION STANDARD DEVIATION
# ------------------------------------------------------------

population_sd <- function(x) {
  sqrt(mean((x - mean(x))^2))
}

# ------------------------------------------------------------
# HELPER FUNCTION: SAFE CORRELATION
# ------------------------------------------------------------

safe_correlation <- function(x, y, method = "pearson") {
  if (length(x) < 2 || length(y) < 2 || sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  
  cor(x, y, method = method, use = "complete.obs")
}

# ------------------------------------------------------------
# HELPER FUNCTION: CALCULATE REVISION METRICS
# ------------------------------------------------------------

calculate_revision_metrics <- function(
    prediction_df, dataset_type, calibration_years, trait,
    valid_iterations, source_file, source_folder
) {
  
  required_columns <- c("measured", "predicted")
  
  if (!all(required_columns %in% names(prediction_df))) {
    stop(
      "Prediction dataframe does not contain required columns: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  dat <- prediction_df %>%
    transmute(
      observed = as.numeric(measured),
      predicted = as.numeric(predicted)
    ) %>%
    filter(is.finite(observed), is.finite(predicted))
  
  n <- nrow(dat)
  
  if (n < 2) {
    stop("Fewer than two valid observed-predicted pairs.")
  }
  
  observed <- dat$observed
  predicted <- dat$predicted
  residuals <- predicted - observed
  
  # ----------------------------------------------------------
  # BASIC DESCRIPTIVE STATISTICS
  # ----------------------------------------------------------
  
  observed_mean <- mean(observed)
  predicted_mean <- mean(predicted)
  
  observed_sd_sample <- sd(observed)
  predicted_sd_sample <- sd(predicted)
  
  observed_sd_population <- population_sd(observed)
  predicted_sd_population <- population_sd(predicted)
  
  observed_min <- min(observed)
  observed_max <- max(observed)
  observed_range <- observed_max - observed_min
  
  # ----------------------------------------------------------
  # ERROR METRICS
  # ----------------------------------------------------------
  
  SSE_PRESS <- sum((observed - predicted)^2)
  SST <- sum((observed - observed_mean)^2)
  MSE <- mean((observed - predicted)^2)
  RMSEP <- sqrt(MSE)
  MAE <- mean(abs(observed - predicted))
  bias <- mean(predicted - observed)
  
  # ----------------------------------------------------------
  # RESIDUAL-BASED R²
  # ----------------------------------------------------------
  
  Residual_R2 <- if (is.finite(SST) && SST > 0) {
    1 - (SSE_PRESS / SST)
  } else {
    NA_real_
  }
  
  # ----------------------------------------------------------
  # CORRELATION METRICS
  # ----------------------------------------------------------
  
  Pearson_r <- safe_correlation(observed, predicted, method = "pearson")
  Spearman_rho <- safe_correlation(observed, predicted, method = "spearman")
  
  # ----------------------------------------------------------
  # REGRESSION DIAGNOSTICS
  # ----------------------------------------------------------
  
  regression_model <- lm(predicted ~ observed)
  
  Regression_Intercept <- unname(coef(regression_model)[1])
  Regression_Slope <- unname(coef(regression_model)[2])
  R2_regression <- summary(regression_model)$r.squared
  Adjusted_R2_regression <- summary(regression_model)$adj.r.squared
  
  # ----------------------------------------------------------
  # NORMALIZED RMSEP
  # ----------------------------------------------------------
  
  nRMSEP_SD_percent <- if (is.finite(observed_sd_sample) && observed_sd_sample > 0) {
    100 * RMSEP / observed_sd_sample
  } else {
    NA_real_
  }
  
  nRMSEP_Range_percent <- if (is.finite(observed_range) && observed_range > 0) {
    100 * RMSEP / observed_range
  } else {
    NA_real_
  }
  
  # ----------------------------------------------------------
  # MSE DECOMPOSITION
  # ----------------------------------------------------------
  
  Squared_Bias_Component <- bias^2
  
  Unequal_Variance_Component <- (
    predicted_sd_population - observed_sd_population
  )^2
  
  Incomplete_Correlation_Component <- if (is.finite(Pearson_r)) {
    2 * predicted_sd_population * observed_sd_population * (1 - Pearson_r)
  } else {
    NA_real_
  }
  
  Decomposition_Sum <- sum(
    Squared_Bias_Component,
    Unequal_Variance_Component,
    Incomplete_Correlation_Component,
    na.rm = FALSE
  )
  
  # ----------------------------------------------------------
  # RELATIVE CONTRIBUTION OF EACH COMPONENT TO MSE
  # ----------------------------------------------------------
  
  Squared_Bias_Percent_MSE <- if (is.finite(MSE) && MSE > 0) {
    100 * Squared_Bias_Component / MSE
  } else {
    NA_real_
  }
  
  Unequal_Variance_Percent_MSE <- if (is.finite(MSE) && MSE > 0) {
    100 * Unequal_Variance_Component / MSE
  } else {
    NA_real_
  }
  
  Incomplete_Correlation_Percent_MSE <- if (
    is.finite(MSE) && MSE > 0 && is.finite(Incomplete_Correlation_Component)
  ) {
    100 * Incomplete_Correlation_Component / MSE
  } else {
    NA_real_
  }
  
  # ----------------------------------------------------------
  # RETURN ONE OUTPUT ROW
  # ----------------------------------------------------------
  
  tibble(
    Calibration_Years = as.character(calibration_years),
    Trait = as.character(trait),
    Dataset = as.character(dataset_type),
    Valid_Iterations = valid_iterations,
    N = n,
    
    Observed_Mean = observed_mean,
    Predicted_Mean = predicted_mean,
    Observed_SD = observed_sd_sample,
    Predicted_SD = predicted_sd_sample,
    Observed_Min = observed_min,
    Observed_Max = observed_max,
    Observed_Range = observed_range,
    
    Residual_R2 = Residual_R2,
    
    Pearson_r = Pearson_r,
    Spearman_rho = Spearman_rho,
    
    R2_regression = R2_regression,
    Adjusted_R2_regression = Adjusted_R2_regression,
    
    Regression_Intercept = Regression_Intercept,
    Regression_Slope = Regression_Slope,
    
    Bias_Predicted_minus_Observed = bias,
    
    MAE = MAE,
    MSE = MSE,
    RMSEP = RMSEP,
    
    nRMSEP_SD_percent = nRMSEP_SD_percent,
    nRMSEP_Range_percent = nRMSEP_Range_percent,
    
    PRESS_SSE = SSE_PRESS,
    SST = SST,
    
    Squared_Bias_Component = Squared_Bias_Component,
    Unequal_Variance_Component = Unequal_Variance_Component,
    Incomplete_Correlation_Component = Incomplete_Correlation_Component,
    
    Squared_Bias_Percent_MSE = Squared_Bias_Percent_MSE,
    Unequal_Variance_Percent_MSE = Unequal_Variance_Percent_MSE,
    Incomplete_Correlation_Percent_MSE = Incomplete_Correlation_Percent_MSE,
    
    Decomposition_Sum = Decomposition_Sum,
    Decomposition_Difference_from_MSE = Decomposition_Sum - MSE,
    
    Source_Folder = basename(source_folder),
    Source_RData = basename(source_file)
  )
}

# ------------------------------------------------------------
# PROCESS ONE PPR2 FOLDER
# ------------------------------------------------------------

process_ppr2_folder <- function(folder) {
  
  cat("\nProcessing:", basename(folder), "\n")
  
  valid_iterations <- extract_valid_iterations(folder)
  
  rdata_files <- list.files(
    path = folder,
    pattern = "_plot_data\\.RData$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  if (length(rdata_files) == 0) {
    warning("No file ending with '_plot_data.RData' found in: ", basename(folder))
    return(NULL)
  }
  
  if (length(rdata_files) > 1) {
    warning(
      "More than one '_plot_data.RData' file found in ",
      basename(folder), ". All matching files will be processed."
    )
  }
  
  map_dfr(rdata_files, function(rdata_file) {
    
    rdata_env <- new.env(parent = emptyenv())
    loaded_objects <- load(rdata_file, envir = rdata_env)
    
    required_objects <- c("ensemble_cal_preds", "ensemble_val_preds")
    missing_objects <- setdiff(required_objects, loaded_objects)
    
    if (length(missing_objects) > 0) {
      warning(
        "Skipping file because required object(s) are missing: ",
        basename(rdata_file), " | Missing: ",
        paste(missing_objects, collapse = ", ")
      )
      return(NULL)
    }
    
    calibration_years <- if ("inYear" %in% loaded_objects) {
      get("inYear", envir = rdata_env)
    } else {
      str_extract(basename(rdata_file), "^[0-9]{4}(?:[-&][0-9]{4})*")
    }
    
    trait <- if ("inVar" %in% loaded_objects) {
      get("inVar", envir = rdata_env)
    } else {
      basename(rdata_file) %>%
        str_remove("^[0-9]{4}(?:[-&][0-9]{4})*_") %>%
        str_remove("_[0-9]{4}-[0-9]{2}-[0-9]{2}_plot_data\\.RData$")
    }
    
    calibration_metrics <- calculate_revision_metrics(
      prediction_df = get("ensemble_cal_preds", envir = rdata_env),
      dataset_type = "Calibration",
      calibration_years = calibration_years,
      trait = trait,
      valid_iterations = valid_iterations,
      source_file = rdata_file,
      source_folder = folder
    )
    
    validation_metrics <- calculate_revision_metrics(
      prediction_df = get("ensemble_val_preds", envir = rdata_env),
      dataset_type = "Within-year validation",
      calibration_years = calibration_years,
      trait = trait,
      valid_iterations = valid_iterations,
      source_file = rdata_file,
      source_folder = folder
    )
    
    file_results <- bind_rows(calibration_metrics, validation_metrics)
    
    file_results_to_save <- file_results %>%
      mutate(across(where(is.numeric), ~ round(.x, 6)))
    
    output_csv <- file.path(
      folder,
      paste0(
        tools::file_path_sans_ext(basename(rdata_file)),
        "_REVISION_METRICS.csv"
      )
    )
    
    write_csv(file_results_to_save, output_csv, na = "")
    cat("Saved/overwritten:", basename(output_csv), "\n")
    
    file_results
  })
}

# ------------------------------------------------------------
# RUN ALL PPR2 FOLDERS
# ------------------------------------------------------------

all_within_year_metrics <- map_dfr(folders, process_ppr2_folder)

# ------------------------------------------------------------
# SAVE COMBINED MASTER OUTPUT
# ------------------------------------------------------------

if (nrow(all_within_year_metrics) > 0) {
  
  all_within_year_metrics <- all_within_year_metrics %>%
    arrange(Calibration_Years, Trait, Dataset)
  
  master_output <- file.path(
    root,
    "PPR2_ALL_WITHIN_YEAR_REVISION_METRICS.csv"
  )
  
  write_csv(
    all_within_year_metrics %>%
      mutate(across(where(is.numeric), ~ round(.x, 6))),
    master_output,
    na = ""
  )
  
  validation_summary <- all_within_year_metrics %>%
    filter(Dataset == "Within-year validation") %>%
    select(
      Calibration_Years,
      Trait,
      Valid_Iterations,
      N,
      
      Residual_R2,
      
      Pearson_r,
      Spearman_rho,
      
      R2_regression,
      Adjusted_R2_regression,
      
      Regression_Intercept,
      Regression_Slope,
      
      Bias_Predicted_minus_Observed,
      
      MAE,
      MSE,
      RMSEP,
      
      nRMSEP_SD_percent,
      nRMSEP_Range_percent,
      
      PRESS_SSE,
      SST,
      
      Observed_Mean,
      Predicted_Mean,
      
      Observed_SD,
      Predicted_SD,
      
      Observed_Min,
      Observed_Max,
      Observed_Range,
      
      Squared_Bias_Component,
      Unequal_Variance_Component,
      Incomplete_Correlation_Component,
      
      Squared_Bias_Percent_MSE,
      Unequal_Variance_Percent_MSE,
      Incomplete_Correlation_Percent_MSE,
      
      Decomposition_Sum,
      Decomposition_Difference_from_MSE,
      
      Source_Folder,
      Source_RData
    ) %>%
    arrange(Calibration_Years, Trait)
  
  validation_summary_output <- file.path(
    root,
    "PPR2_WITHIN_YEAR_VALIDATION_REVISION_SUMMARY.csv"
  )
  
  write_csv(
    validation_summary %>%
      mutate(across(where(is.numeric), ~ round(.x, 6))),
    validation_summary_output,
    na = ""
  )
  
  cat("\n")
  cat("============================================\n")
  cat("PROCESSING COMPLETE\n")
  cat("============================================\n")
  cat("Folders searched:", length(folders), "\n")
  cat("Total metric rows:", nrow(all_within_year_metrics), "\n")
  cat("Validation-summary rows:", nrow(validation_summary), "\n")
  cat("Master metrics file:\n", master_output, "\n")
  cat("Validation summary file:\n", validation_summary_output, "\n")
  cat("============================================\n")
  
  print(validation_summary, n = Inf, width = Inf)
  
} else {
  warning("No valid revision metrics were produced.")
}



### END ###

