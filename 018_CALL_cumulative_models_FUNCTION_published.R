# 07.08.2026
# Written by Nargiz Rüter

# ============================================================
# PPR2 — GENERATE REVISED CROSS-YEAR TRANSFER PREDICTIONS
# ============================================================
#
# IMPORTANT:
# Cross-year predictions are generated ONLY when the target
# dataset shares NO year with the calibration dataset.
#
# Examples:
#
# Calibration       Target          Action
# ------------------------------------------------------------
# 2010              2013            KEEP
# 2010              2010-2012       SKIP
# 2010-2011         2012            KEEP
# 2010-2011         2011-2013       SKIP
# 2010-2011-2012    2013            KEEP
# 2010-2011-2012    2012-2013       SKIP
#
# Calibration models are read from:
#   ./r_model_output/PPR2_...
#
# Transfer predictions are written to:
#   ./r_model_output/Transfer_Models_Revised/
#
# Master statistics are written to:
#   ./r_model_output/
#
# ============================================================

# ============================================================
# PPR2 — REVISED CROSS-YEAR TRANSFER PREDICTIONS
# ============================================================

setwd("C:/Users/.../")

library(tidyverse)
library(stringr)
library(purrr)

study_traits <- c("ADF_perc","ADL_perc","Biom_dry_g","Biom_wet_g","C_perc","N_perc","NDF_perc")

base_output_dir <- "./r_model_output/"
input_data_dir <- "./input_data/"
transfer_output_dir <- file.path(base_output_dir, "Transfer_Models_Revised")
dir.create(transfer_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# YEAR HELPERS
# ------------------------------------------------------------

year_map <- c(
  "1011"="2010-2011", "1012"="2010-2012", "1013"="2010-2013",
  "1112"="2011-2012", "1113"="2011-2013", "1213"="2012-2013",
  "2110"="2010-2011-2012", "3110"="2010-2011-2013",
  "3210"="2010-2012-2013", "3211"="2011-2012-2013",
  "2010"="2010", "2011"="2011", "2012"="2012", "2013"="2013"
)

normalize_year_label <- function(x) {
  x <- as.character(x)
  if (length(x) == 1 && x %in% names(year_map)) return(unname(year_map[x]))
  x
}

extract_years <- function(x) {
  unique(str_extract_all(as.character(x), "2010|2011|2012|2013")[[1]])
}

years_overlap <- function(calibration_years, target_years) {
  length(intersect(extract_years(calibration_years), extract_years(target_years))) > 0
}

extract_input_year_label <- function(input_file) {
  filename <- basename(input_file)
  raw_code <- str_match(filename, "_input_([^_]+)_REF")[1,2]
  
  if (is.na(raw_code)) {
    raw_code <- str_extract(filename, paste(names(year_map), collapse="|"))
  }
  
  if (is.na(raw_code)) stop("Could not identify target year from: ", filename)
  
  normalize_year_label(raw_code)
}

# ------------------------------------------------------------
# FIND / PARSE MODEL FOLDERS
# ------------------------------------------------------------

parse_ppr2_folder <- function(folder) {
  parsed <- str_match(basename(folder), "^PPR2_([^_]+)_(\\d{4}-\\d{2}-\\d{2})_(.+)$")
  
  if (any(is.na(parsed[1,2:4]))) {
    return(tibble(full_path=folder, folder_name=basename(folder),
                  Calibration_Years=NA_character_, Run_Date=as.Date(NA), Trait=NA_character_))
  }
  
  tibble(
    full_path = folder,
    folder_name = basename(folder),
    Calibration_Years = normalize_year_label(parsed[1,2]),
    Run_Date = as.Date(parsed[1,3]),
    Trait = parsed[1,4]
  )
}

extract_valid_iterations <- function(folder) {
  files <- list.files(
    folder,
    pattern = "_REF_DF_MODELs-[0-9]{1,4}-STATISTICS\\.csv$",
    full.names = FALSE
  )
  
  if (length(files) == 0) {
    warning("No valid-iteration statistics file found in: ", basename(folder))
    return(NA_integer_)
  }
  
  if (length(files) > 1) {
    warning("More than one valid-iteration statistics file found in: ", basename(folder))
  }
  
  n_iter <- as.integer(
    str_match(files[1], "_REF_DF_MODELs-([0-9]{1,4})-STATISTICS\\.csv$")[1,2]
  )
  
  if (is.na(n_iter) || n_iter < 1 || n_iter > 1000) {
    warning("Invalid iteration count in: ", basename(folder))
    return(NA_integer_)
  }
  
  n_iter
}

all_model_folders <- list.dirs(base_output_dir, recursive=FALSE, full.names=TRUE)
all_model_folders <- all_model_folders[str_detect(basename(all_model_folders), "^PPR2_")]

trait_models_all <- map_dfr(all_model_folders, parse_ppr2_folder) %>%
  filter(!is.na(Calibration_Years), !is.na(Run_Date), !is.na(Trait), Trait %in% study_traits)

# Keep only newest rerun folder for each calibration-year × trait
trait_models <- trait_models_all %>%
  group_by(Calibration_Years, Trait) %>%
  slice_max(Run_Date, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  arrange(Calibration_Years, Trait)

trait_models <- trait_models %>%
  mutate(
    Valid_Iterations = map_int(full_path, extract_valid_iterations)
  )

excluded_models <- trait_models %>%
  filter(is.na(Valid_Iterations) | Valid_Iterations < 700)

write.csv(
  excluded_models,
  file.path(base_output_dir, "PPR2_TRANSFER_MODELS_EXCLUDED_VALID_ITERATIONS.csv"),
  row.names = FALSE
)

trait_models <- trait_models %>%
  filter(!is.na(Valid_Iterations), Valid_Iterations >= 700) %>%
  arrange(Calibration_Years, Trait)

cat("Calibration folders retained after Valid_Iterations >= 700 filter:",
    nrow(trait_models), "\n")

print(
  trait_models %>%
    select(Calibration_Years, Trait, Valid_Iterations, Run_Date, folder_name),
  n = Inf
)

# ------------------------------------------------------------
# TRANSFER FUNCTION
# ------------------------------------------------------------

# predict_with_cumulative_models_revised <- function(
    #     input_data_dir, coeffs_data_dir, output_base_dir,
#     trait_name, calibration_years, intercept_col=2) {

predict_with_cumulative_models_revised <- function(
    input_data_dir, coeffs_data_dir, output_base_dir,
    trait_name, calibration_years, valid_iterations, intercept_col=2) {
  
  input_files <- list.files(input_data_dir, pattern="^No_Duplicates_.*\\.csv$", full.names=TRUE)
  coeff_files <- list.files(coeffs_data_dir, pattern="_MAIN_MODELs-avg-COEFS\\.csv$", full.names=TRUE)
  
  if (length(coeff_files) == 0) {
    warning("No coefficient file found in: ", basename(coeffs_data_dir))
    return(tibble())
  }
  
  calibration_years <- normalize_year_label(calibration_years)
  all_metrics <- tibble()
  
  for (input_file in input_files) {
    
    target_years <- extract_input_year_label(input_file)
    
    # Skip if target shares ANY year with calibration
    if (years_overlap(calibration_years, target_years)) {
      cat("SKIP:", calibration_years, "->", target_years, "|", trait_name, "\n")
      next
    }
    
    cat("KEEP:", calibration_years, "->", target_years, "|", trait_name, "\n")
    
    year_data <- read.csv(input_file)
    
    if (!(trait_name %in% names(year_data)) || all(is.na(year_data[[trait_name]]))) {
      cat("Trait unavailable:", trait_name, "|", target_years, "\n")
      next
    }
    
    year_data <- year_data %>%
      mutate(
        nRow = row_number(),
        Unique_ID = paste0(nRow, "_", Plot_ID, "_", target_years),
        Year = target_years
      ) %>%
      relocate(Year, Unique_ID)
    
    year_data_filtered <- year_data[complete.cases(year_data[[trait_name]]), ]
    
    if (nrow(year_data_filtered) < 2) next
    
    wvl_cols <- grep("^X\\d+", names(year_data_filtered), value=TRUE)
    spectra_matrix <- as.matrix(year_data_filtered[, wvl_cols, drop=FALSE])
    
    for (coeff_file in coeff_files) {
      
      coeff_data <- read.csv(coeff_file)
      intercept <- coeff_data[1, intercept_col]
      coeffs <- coeff_data[-1, -c(1,3), drop = FALSE]
      
      if (nrow(coeffs) != ncol(spectra_matrix)) {
        stop("Coefficient/wavelength mismatch for ",
             calibration_years, " / ", trait_name, " -> ", target_years)
      }
      
      predictions <- intercept + spectra_matrix %*% as.matrix(coeffs)
      pred_col <- paste0("Cumulative_Model_Preds_", trait_name)
      year_data_filtered[[pred_col]] <- as.vector(predictions)
      
      output_dir <- file.path(
        output_base_dir,
        paste0("Predict-From-Cumulative_", calibration_years, "_", trait_name)
      )
      dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
      
      prediction_file <- file.path(
        output_dir,
        paste0("Predicted-From_", calibration_years,
               "_Model_", trait_name, "_", target_years, ".csv")
      )
      
      write.csv(year_data_filtered, prediction_file, row.names=FALSE)
      
      observed <- year_data_filtered[[trait_name]]
      predicted <- as.vector(predictions)
      
      mse <- mean((observed - predicted)^2)
      rmse <- sqrt(mse)
      mae <- mean(abs(observed - predicted))
      SSE <- sum((observed - predicted)^2)
      SST <- sum((observed - mean(observed))^2)
      R2 <- ifelse(SST > 0, 1 - SSE/SST, NA_real_)
      bias <- mean(predicted - observed)
      
      nonzero <- observed != 0
      mape <- ifelse(
        any(nonzero),
        mean(abs((observed[nonzero] - predicted[nonzero]) / observed[nonzero])) * 100,
        NA_real_
      )
      
      metrics <- tibble(
        Predicted_From = calibration_years,
        Input_Year = target_years,
        Trait = trait_name,
        Valid_Iterations = valid_iterations,
        N = length(observed),
        MAE = round(mae, 3),
        RMSE = round(rmse, 3),
        RMSEP = round(rmse, 3),
        MAPE_percent = round(mape, 2),
        Bias = round(bias, 3),
        PRESS = round(SSE, 1),
        R2 = round(R2, 6),
        MSE = round(mse, 3)
      )
      
      stats_file <- file.path(
        output_dir,
        paste0("Predicted-From_", calibration_years,
               "_Model_", trait_name, "_", target_years, "_STATS.csv")
      )
      
      write.csv(metrics, stats_file, row.names=FALSE)
      all_metrics <- bind_rows(all_metrics, metrics)
    }
  }
  
  all_metrics
}

# ------------------------------------------------------------
# RUN ALL VALID TRANSFERS
# ------------------------------------------------------------

all_results <- tibble()

for (i in seq_len(nrow(trait_models))) {
  
  coeffs_data_dir <- trait_models$full_path[i]
  trait_name <- trait_models$Trait[i]
  calibration_years <- trait_models$Calibration_Years[i]
  valid_iterations <- trait_models$Valid_Iterations[i]
  
  cat("\nCalibration:", calibration_years,
      "| Trait:", trait_name,
      "| Valid iterations:", valid_iterations,
      "| Folder:", basename(coeffs_data_dir), "\n")
  
  metrics_result <- predict_with_cumulative_models_revised(
    input_data_dir = input_data_dir,
    coeffs_data_dir = coeffs_data_dir,
    output_base_dir = transfer_output_dir,
    trait_name = trait_name,
    calibration_years = calibration_years,
    valid_iterations = valid_iterations
  )
  
  all_results <- bind_rows(all_results, metrics_result)
}

# ------------------------------------------------------------
# FINAL OVERLAP SAFETY CHECK
# ------------------------------------------------------------

if (nrow(all_results) > 0) {
  
  valid_transfer <- map2_lgl(
    all_results$Predicted_From,
    all_results$Input_Year,
    ~ !years_overlap(.x, .y)
  )
  
  if (!all(valid_transfer)) {
    stop("ERROR: At least one overlapping calibration-target case remains.")
  }
  
  cat("\nAll transfer predictions are temporally independent: TRUE\n")
}

# ------------------------------------------------------------
# FORMAT MASTER TABLE
# ------------------------------------------------------------

year_order <- c(
  "2010","2011","2012","2013",
  "2010-2011","2010-2012","2010-2013",
  "2011-2012","2011-2013","2012-2013",
  "2010-2011-2012","2010-2011-2013",
  "2010-2012-2013","2011-2012-2013"
)

all_results_filtered <- all_results %>%
  filter(Trait %in% study_traits) %>%
  mutate(
    Predicted_From = factor(Predicted_From, levels=year_order),
    Input_Year = factor(Input_Year, levels=year_order),
    Rounded_R2 = sprintf("%.2f", R2)
  ) %>%
  arrange(Trait, Predicted_From, Input_Year)

# ------------------------------------------------------------
# SAVE MASTER SUMMARY IN r_model_output
# ------------------------------------------------------------

master_summary_file <- file.path(
  base_output_dir,
  "PPR2_All_Traits_Master_Prediction_Statistics_ALL-MAIN_Revised.csv"
)

write.csv(all_results_filtered, master_summary_file, row.names=FALSE)

cat("\n============================================\n")
cat("TRANSFER RUN COMPLETE\n")
cat("Valid transfer predictions:", nrow(all_results_filtered), "\n")
cat("Transfer folders:", transfer_output_dir, "\n")
cat("Master summary:", master_summary_file, "\n")
cat("============================================\n")


### END ###

