# 04.08.2026
# Written by Nargiz Rüter
# 
# ============================================================
# PPR2 REVISION — CROSS-YEAR METRICS FROM PREDICTION CSV FILES
# ============================================================
#
# Reads prediction CSVs only, excludes _STATS.csv and
# _REVISION_METRICS.csv files, calculates harmonized cross-year
# revision metrics, saves individual metrics beside prediction
# CSVs, and saves combined summaries in the root folder.
# ============================================================

library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

# ------------------------------------------------------------
# ROOT DIRECTORY
# ------------------------------------------------------------

root <- paste0(
  "C:/Users/.../")

options(scipen = 999)

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

OVERWRITE_REVISION_OUTPUTS <- TRUE
ROUND_DIGITS <- 6

valid_years <- c("2010", "2011", "2012", "2013")

# ------------------------------------------------------------
# FIND MODEL FOLDERS
# ------------------------------------------------------------

model_folders <- list.dirs(root, full.names = TRUE, recursive = FALSE)
model_folders <- model_folders[
  str_detect(basename(model_folders), "^Predict-From-Cumulative_")
]

cat("Cross-year model folders found:", length(model_folders), "\n")

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

extract_years <- function(x) {
  years <- str_extract_all(as.character(x), "2010|2011|2012|2013")[[1]]
  unique(years)
}

population_sd <- function(x) {
  sqrt(mean((x - mean(x))^2))
}

safe_correlation <- function(x, y, method = "pearson") {
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  
  if (length(x) < 2 || length(y) < 2 ||
      !is.finite(sd(x)) || !is.finite(sd(y)) ||
      sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  
  suppressWarnings(cor(x, y, method = method, use = "complete.obs"))
}

# ------------------------------------------------------------
# PARSE MODEL FOLDER NAME
# ------------------------------------------------------------

parse_model_folder <- function(folder) {
  folder_name <- basename(folder)
  
  pattern <- paste0(
    "^Predict-From-Cumulative_",
    "((?:2010|2011|2012|2013)(?:-(?:2010|2011|2012|2013))*)",
    "_(.+)$"
  )
  
  parsed <- str_match(folder_name, pattern)
  
  if (nrow(parsed) == 0 || is.na(parsed[1,2]) || is.na(parsed[1,3])) {
    return(list(
      Calibration_Years = NA_character_,
      Trait = NA_character_
    ))
  }
  
  list(
    Calibration_Years = parsed[1,2],
    Trait = parsed[1,3]
  )
}

# ------------------------------------------------------------
# EXTRACT TARGET YEARS FROM PREDICTION FILENAME
# ------------------------------------------------------------

extract_target_year_label <- function(filename) {
  filename_no_extension <- tools::file_path_sans_ext(basename(filename))
  
  pattern <- paste0(
    "_((?:2010|2011|2012|2013)",
    "(?:-(?:2010|2011|2012|2013))*)$"
  )
  
  target_match <- str_match(filename_no_extension, pattern)
  
  if (nrow(target_match) == 0 || is.na(target_match[1,2])) {
    return(NA_character_)
  }
  
  target_match[1,2]
}

# ------------------------------------------------------------
# CLASSIFY SINGLE- OR COMBINED-YEAR DATASET
# ------------------------------------------------------------

classify_year_set <- function(year_label) {
  years <- extract_years(year_label)
  
  if (length(years) == 1) return("Single-year")
  if (length(years) > 1) return("Combined-year")
  
  NA_character_
}

# ------------------------------------------------------------
# CALCULATE REVISION METRICS
# ------------------------------------------------------------

calculate_revision_metrics <- function(
    observed, predicted, calibration_years, target_years,
    trait, valid_iterations, source_folder, source_prediction_csv) {
  
  dat <- tibble(
    observed = suppressWarnings(as.numeric(observed)),
    predicted = suppressWarnings(as.numeric(predicted))
  ) %>%
    filter(is.finite(observed), is.finite(predicted))
  
  n <- nrow(dat)
  
  if (n < 2) stop("Fewer than two valid observed-predicted pairs.")
  
  observed <- dat$observed
  predicted <- dat$predicted
  
  # Descriptive statistics
  observed_mean <- mean(observed)
  predicted_mean <- mean(predicted)
  
  observed_sd_sample <- sd(observed)
  predicted_sd_sample <- sd(predicted)
  
  observed_sd_population <- population_sd(observed)
  predicted_sd_population <- population_sd(predicted)
  
  observed_min <- min(observed)
  observed_max <- max(observed)
  observed_range <- observed_max - observed_min
  
  predicted_min <- min(predicted)
  predicted_max <- max(predicted)
  predicted_range <- predicted_max - predicted_min
  
  # Error metrics
  residuals_predicted_minus_observed <- predicted - observed
  
  PRESS_SSE <- sum((observed - predicted)^2)
  SST <- sum((observed - observed_mean)^2)
  MSE <- mean((observed - predicted)^2)
  
  RMSE <- sqrt(MSE)
  RMSEP <- RMSE
  MAE <- mean(abs(observed - predicted))
  Bias_Predicted_minus_Observed <- mean(residuals_predicted_minus_observed)
  
  # MAPE
  nonzero_observed <- observed != 0
  
  MAPE_percent <- if (any(nonzero_observed)) {
    mean(abs(
      (observed[nonzero_observed] - predicted[nonzero_observed]) /
        observed[nonzero_observed]
    )) * 100
  } else {
    NA_real_
  }
  
  # Residual-based R²
  Residual_R2 <- if (is.finite(SST) && SST > 0) {
    1 - (PRESS_SSE / SST)
  } else {
    NA_real_
  }
  
  # Correlations
  Pearson_r <- safe_correlation(observed, predicted, method = "pearson")
  Spearman_rho <- safe_correlation(observed, predicted, method = "spearman")
  
  # Regression diagnostics
  regression_model <- lm(predicted ~ observed)
  
  Regression_Intercept <- unname(coef(regression_model)[1])
  Regression_Slope <- unname(coef(regression_model)[2])
  R2_regression <- summary(regression_model)$r.squared
  Adjusted_R2_regression <- summary(regression_model)$adj.r.squared
  
  # Normalized RMSEP
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
  
  # MSE decomposition
  Squared_Bias_Component <- Bias_Predicted_minus_Observed^2
  
  Unequal_Variance_Component <- (
    predicted_sd_population - observed_sd_population
  )^2
  
  Incomplete_Correlation_Component <- if (is.finite(Pearson_r)) {
    2 * predicted_sd_population * observed_sd_population * (1 - Pearson_r)
  } else {
    NA_real_
  }
  
  Decomposition_Sum <- if (all(is.finite(c(
    Squared_Bias_Component,
    Unequal_Variance_Component,
    Incomplete_Correlation_Component
  )))) {
    Squared_Bias_Component +
      Unequal_Variance_Component +
      Incomplete_Correlation_Component
  } else {
    NA_real_
  }
  
  # Decomposition percentages
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
  
  # Calibration/target year overlap
  calibration_year_vector <- extract_years(calibration_years)
  target_year_vector <- extract_years(target_years)
  overlap_years <- intersect(calibration_year_vector, target_year_vector)
  
  tibble(
    Calibration_Years = as.character(calibration_years),
    Calibration_Type = classify_year_set(calibration_years),
    Target_Years = as.character(target_years),
    Target_Type = classify_year_set(target_years),
    Evaluation_Type = "Cross-year",
    Trait = as.character(trait),
    Valid_Iterations = valid_iterations,
    N = n,
    
    Year_Overlap_Count = length(overlap_years),
    Overlap_Years = paste(overlap_years, collapse = "-"),
    
    Observed_Mean = observed_mean,
    Predicted_Mean = predicted_mean,
    Observed_SD = observed_sd_sample,
    Predicted_SD = predicted_sd_sample,
    Observed_Min = observed_min,
    Observed_Max = observed_max,
    Observed_Range = observed_range,
    Predicted_Min = predicted_min,
    Predicted_Max = predicted_max,
    Predicted_Range = predicted_range,
    
    Residual_R2 = Residual_R2,
    Pearson_r = Pearson_r,
    Spearman_rho = Spearman_rho,
    
    R2_regression = R2_regression,
    Adjusted_R2_regression = Adjusted_R2_regression,
    Regression_Intercept = Regression_Intercept,
    Regression_Slope = Regression_Slope,
    
    Bias_Predicted_minus_Observed = Bias_Predicted_minus_Observed,
    
    MAE = MAE,
    MSE = MSE,
    RMSE = RMSE,
    RMSEP = RMSEP,
    MAPE_percent = MAPE_percent,
    
    nRMSEP_SD_percent = nRMSEP_SD_percent,
    nRMSEP_Range_percent = nRMSEP_Range_percent,
    
    PRESS_SSE = PRESS_SSE,
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
    Source_Prediction_CSV = basename(source_prediction_csv)
  )
}

# ------------------------------------------------------------
# FIND PREDICTION CSV FILES
# ------------------------------------------------------------

find_prediction_csvs <- function(folder) {
  all_csvs <- list.files(
    folder,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  
  all_csvs[
    str_detect(basename(all_csvs), "^Predicted-From_") &
      !str_detect(basename(all_csvs), "_STATS\\.csv$") &
      !str_detect(basename(all_csvs), "_REVISION_METRICS\\.csv$")
  ]
}

# ------------------------------------------------------------
# BUILD PREDICTION INVENTORY
# ------------------------------------------------------------

prediction_inventory <- map_dfr(model_folders, function(folder) {
  parsed_folder <- parse_model_folder(folder)
  prediction_files <- find_prediction_csvs(folder)
  
  if (length(prediction_files) == 0) {
    return(tibble(
      Model_Folder_Path = character(0),
      Model_Folder = character(0),
      Calibration_Years = character(0),
      Trait = character(0),
      Prediction_CSV = character(0)
    ))
  }
  
  tibble(
    Model_Folder_Path = folder,
    Model_Folder = basename(folder),
    Calibration_Years = parsed_folder$Calibration_Years,
    Trait = parsed_folder$Trait,
    Prediction_CSV = prediction_files
  )
})

cat("Prediction CSV files found:", nrow(prediction_inventory), "\n")

# ------------------------------------------------------------
# PROCESS ONE PREDICTION CSV
# ------------------------------------------------------------

process_prediction_csv <- function(
    prediction_csv, model_folder, calibration_years, trait) {
  
  filename <- basename(prediction_csv)
  target_years <- extract_target_year_label(filename)
  
  stats_file <- file.path(
    model_folder,
    paste0(tools::file_path_sans_ext(filename), "_STATS.csv")
  )
  
  valid_iterations <- if (file.exists(stats_file)) {
    stats_dat <- read_csv(stats_file, show_col_types = FALSE)
    
    if ("Valid_Iterations" %in% names(stats_dat)) {
      stats_dat$Valid_Iterations[1]
    } else {
      NA_integer_
    }
  } else {
    NA_integer_
  }
  
  # Validate metadata
  if (is.na(calibration_years) || calibration_years == "") {
    stop("Could not identify calibration years from folder: ", basename(model_folder))
  }
  
  if (is.na(trait) || trait == "") {
    stop("Could not identify trait from folder: ", basename(model_folder))
  }
  
  if (is.na(target_years) || target_years == "") {
    stop("Could not identify target years from filename: ", filename)
  }
  
  # Read prediction file
  dat <- read_csv(prediction_csv, show_col_types = FALSE, progress = FALSE)
  
  observed_column <- trait
  expected_predicted_column <- paste0("Cumulative_Model_Preds_", trait)
  
  if (!(observed_column %in% names(dat))) {
    stop("Observed trait column '", observed_column, "' not found in file: ", filename)
  }
  
  # Identify predicted column
  if (expected_predicted_column %in% names(dat)) {
    predicted_column <- expected_predicted_column
    
  } else {
    candidate_predicted_columns <- names(dat)[
      str_detect(names(dat), "^Cumulative_Model_Preds_")
    ]
    
    if (length(candidate_predicted_columns) == 1) {
      predicted_column <- candidate_predicted_columns[1]
      
      warning(
        "Expected prediction column '", expected_predicted_column,
        "' not found in ", filename, ". Using '",
        predicted_column, "' instead."
      )
      
    } else {
      stop("Could not uniquely identify the predicted-values column in file: ", filename)
    }
  }
  
  # Calculate metrics
  metrics <- calculate_revision_metrics(
    observed = dat[[observed_column]],
    predicted = dat[[predicted_column]],
    calibration_years = calibration_years,
    target_years = target_years,
    trait = trait,
    valid_iterations = valid_iterations,
    source_folder = model_folder,
    source_prediction_csv = prediction_csv
  )
  
  # Output path
  output_filename <- paste0(
    tools::file_path_sans_ext(filename),
    "_REVISION_METRICS.csv"
  )
  
  output_path <- file.path(model_folder, output_filename)
  
  if (file.exists(output_path) && !OVERWRITE_REVISION_OUTPUTS) {
    stop("Revision output already exists and overwriting is disabled: ", output_path)
  }
  
  write_csv(
    metrics %>% mutate(across(where(is.numeric), ~ round(.x, ROUND_DIGITS))),
    output_path,
    na = ""
  )
  
  metrics
}

# ------------------------------------------------------------
# PROCESS ALL FILES WITH PROGRESS BAR
# ------------------------------------------------------------

total_files <- nrow(prediction_inventory)

cat("\nProcessing", total_files, "cross-year prediction CSV files...\n")

progress_bar <- txtProgressBar(min = 0, max = total_files, style = 3)

results_list <- vector(mode = "list", length = total_files)
log_list <- vector(mode = "list", length = total_files)

for (i in seq_len(total_files)) {
  inventory_row <- prediction_inventory[i, ]
  
  processing_result <- tryCatch(
    {
      metrics <- process_prediction_csv(
        prediction_csv = inventory_row$Prediction_CSV,
        model_folder = inventory_row$Model_Folder_Path,
        calibration_years = inventory_row$Calibration_Years,
        trait = inventory_row$Trait
      )
      
      list(
        metrics = metrics,
        log = tibble(
          Model_Folder = inventory_row$Model_Folder,
          Prediction_CSV = basename(inventory_row$Prediction_CSV),
          Calibration_Years = inventory_row$Calibration_Years,
          Target_Years = extract_target_year_label(inventory_row$Prediction_CSV),
          Trait = inventory_row$Trait,
          Processing_Success = TRUE,
          Status = "Revision metrics calculated and saved"
        )
      )
    },
    
    error = function(e) {
      list(
        metrics = NULL,
        log = tibble(
          Model_Folder = inventory_row$Model_Folder,
          Prediction_CSV = basename(inventory_row$Prediction_CSV),
          Calibration_Years = inventory_row$Calibration_Years,
          Target_Years = extract_target_year_label(inventory_row$Prediction_CSV),
          Trait = inventory_row$Trait,
          Processing_Success = FALSE,
          Status = conditionMessage(e)
        )
      )
    }
  )
  
  results_list[[i]] <- processing_result$metrics
  log_list[[i]] <- processing_result$log
  
  setTxtProgressBar(progress_bar, i)
}

close(progress_bar)
cat("\nCross-year revision-metric processing finished.\n")

# ------------------------------------------------------------
# COMBINE RESULTS AND LOG
# ------------------------------------------------------------

all_cross_year_metrics <- bind_rows(results_list)
processing_log <- bind_rows(log_list)

# ------------------------------------------------------------
# CHECK FOR YEAR OVERLAPS
# ------------------------------------------------------------

if (
  nrow(all_cross_year_metrics) > 0 &&
  any(all_cross_year_metrics$Year_Overlap_Count > 0, na.rm = TRUE)
) {
  warning(
    "At least one processed prediction still contains overlapping ",
    "calibration and target years. Inspect Year_Overlap_Count and Overlap_Years."
  )
}

# ------------------------------------------------------------
# SORT COMPLETE OUTPUT
# ------------------------------------------------------------

all_cross_year_metrics <- all_cross_year_metrics %>%
  arrange(Trait, Calibration_Type, Calibration_Years, Target_Years)

# ------------------------------------------------------------
# SAVE COMPLETE MASTER OUTPUT
# ------------------------------------------------------------

master_output_path <- file.path(
  root,
  "PPR2_ALL_CROSS_YEAR_REVISION_METRICS.csv"
)

write_csv(
  all_cross_year_metrics %>%
    mutate(across(where(is.numeric), ~ round(.x, ROUND_DIGITS))),
  master_output_path,
  na = ""
)

# ------------------------------------------------------------
# CREATE REVIEWER-FACING SUMMARY
# ------------------------------------------------------------

cross_year_revision_summary <- all_cross_year_metrics %>%
  select(
    Calibration_Years,
    Calibration_Type,
    Target_Years,
    Target_Type,
    Evaluation_Type,
    Trait,
    Valid_Iterations,
    N,
    
    Residual_R2,
    
    Pearson_r,
    Spearman_rho,
    
    Bias_Predicted_minus_Observed,
    
    Regression_Intercept,
    Regression_Slope,
    R2_regression,
    Adjusted_R2_regression,
    
    MAE,
    RMSE,
    RMSEP,
    nRMSEP_SD_percent,
    nRMSEP_Range_percent,
    
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
    
    Year_Overlap_Count,
    Overlap_Years,
    
    Source_Folder,
    Source_Prediction_CSV
  ) %>%
  arrange(Trait, Calibration_Type, Calibration_Years, Target_Years)

# ------------------------------------------------------------
# SAVE REVIEWER-FACING SUMMARY
# ------------------------------------------------------------

summary_output_path <- file.path(
  root,
  "PPR2_CROSS_YEAR_REVISION_SUMMARY.csv"
)

write_csv(
  cross_year_revision_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, ROUND_DIGITS))),
  summary_output_path,
  na = ""
)

# ------------------------------------------------------------
# SAVE PROCESSING LOG
# ------------------------------------------------------------

processing_log_path <- file.path(
  root,
  "PPR2_CROSS_YEAR_REVISION_PROCESSING_LOG.csv"
)

write_csv(processing_log, processing_log_path, na = "")

# ------------------------------------------------------------
# CREATE FAILED-FILE LOG
# ------------------------------------------------------------

failed_files <- processing_log %>%
  filter(!Processing_Success)

failed_files_path <- file.path(
  root,
  "PPR2_CROSS_YEAR_REVISION_FAILED_FILES.csv"
)

write_csv(failed_files, failed_files_path, na = "")

# ------------------------------------------------------------
# COMPLETION SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("CROSS-YEAR REVISION METRICS COMPLETE\n")
cat("============================================\n")
cat("Model folders searched:", length(model_folders), "\n")
cat("Prediction files identified:", nrow(prediction_inventory), "\n")
cat("Prediction files successfully processed:",
    sum(processing_log$Processing_Success), "\n")
cat("Prediction files that failed:",
    sum(!processing_log$Processing_Success), "\n")
cat("Individual revision-metrics files were saved inside their respective model folders.\n")

cat("\nComplete master output:\n", master_output_path, "\n")
cat("\nCompact revision summary:\n", summary_output_path, "\n")
cat("\nProcessing log:\n", processing_log_path, "\n")
cat("\nFailed-file log:\n", failed_files_path, "\n")
cat("============================================\n")

# ------------------------------------------------------------
# DISPLAY SUMMARY
# ------------------------------------------------------------

print(cross_year_revision_summary, n = Inf, width = Inf)

# ------------------------------------------------------------
# QUICK CHECKS
# ------------------------------------------------------------

overlap_check <- all(all_cross_year_metrics$Year_Overlap_Count == 0)

cat("\nNo overlapping calibration/target years:", overlap_check, "\n")

cat(
  "Cross-year Residual R² range:",
  paste(range(all_cross_year_metrics$Residual_R2, na.rm = TRUE), collapse = " to "),
  "\n"
)

residual_r2_counts <- all_cross_year_metrics %>%
  summarise(
    Total = n(),
    Positive_Residual_R2 = sum(Residual_R2 > 0, na.rm = TRUE),
    Nonpositive_Residual_R2 = sum(Residual_R2 <= 0, na.rm = TRUE),
    Missing_Residual_R2 = sum(is.na(Residual_R2))
  )

print(residual_r2_counts)



### END ###

