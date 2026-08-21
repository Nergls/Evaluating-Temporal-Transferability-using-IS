### 24.04.2025 -> updated 05.08.2026
### Written by Nargiz Rüter

### PLSR for multi year IS data ###
### The code contains loops to build PLSR models in 2010, 2011, 2012, 2013 & year combinations for following traits:
                                                                                                                    # "ADF_perc", "ADL_perc",
                                                                                                                    # "Biom_dry_g", "Biom_wet_g",
                                                                                                                    # "C_perc", "N_perc",
                                                                                                                    # "NDF_perc"
## Traits' abbreviations:
# ADF_perc: Percentage of Acid Detergent Fiber
# ADL_perc: Percentage of Acid Detergent Lignin
# Biom_dry_g: Above-ground dry biomass (g m⁻²)
# Biom_wet_g: Above-ground fresh biomass (g m⁻²)
# C_perc: Percentage of Carbon
# N_perc: Percentage of Nitrogen
# NDF_perc: Percentage of Neutral Detergent Fiber

#### Load libraries needed to run the script #####
list.of.packages <- c("pls", "dplyr", "reshape2", "ggplot2", "spectratrait", "seecolor", "magrittr","openxlsx")
invisible(lapply(list.of.packages, library, character.only = TRUE))
# install.packages("remotes")
# remotes::install_github("TESTgroup-BNL/spectratrait")

# Avoid scientific notation globally
options(scipen = 999)

setwd("C:/Users/...")

# Input data reading
# REF_DF_LIST <- list.files("./input_data/" ,pattern = "No_Duplicates_Combined-Years_Model_input_") # combined years
REF_DF_LIST <- list.files("./input_data/" ,pattern = "No_Duplicates_Single-Year_Model_input_")   # single years
REF_DF_LIST_ALL <- data.frame(File_Name = REF_DF_LIST)
REF_DF_LIST_ALL$Year <- apply(REF_DF_LIST_ALL, 1, function(x) sub(".*_Model_input_([0-9]{4})_.*", "\\1", x[1]))

print(REF_DF_LIST_ALL)

# Define trait mappings for each file/year combination
study_traits <- c(
  "ADF_perc",
  "ADL_perc",
  "Biom_dry_g",
  "Biom_wet_g",
  "C_perc",
  "N_perc",
  "NDF_perc"
)

# Single-Year
year_trait_map <- list(
  "2010" = study_traits,
  "2011" = study_traits,
  "2012" = study_traits,
  "2013" = study_traits
)

# # Define trait mappings for each file/year combination
# year_trait_map <- list(
#   "1011" = study_traits,
#   "1012" = study_traits,
#   "1013" = study_traits,
#   "1112" = study_traits,
#   "1113" = study_traits,
#   "1213" = study_traits,
#   "2110" = study_traits,
#   "3110" = study_traits,
#   "3210" = study_traits,
#   "3211" = study_traits
#   ## "3333" = study_traits, # This combines 4 years, it has never been and will never be used in the study
# )

# Mapping from internal year code to readable string
year_label_map <- c(
  "2010" = "2010",
  "2011" = "2011",
  "2012" = "2012",
  "2013" = "2013"
)

# # # Mapping from internal year code to readable string
# year_label_map <- c(
#   "1011" = "2010-2011",
#   "1012" = "2010-2012",
#   "1013" = "2010-2013",
#   "1112" = "2011-2012",
#   "1113" = "2011-2013",
#   "1213" = "2012-2013",
#   "2110" = "2010-2011-2012",
#   "3110" = "2010-2011-2013",
#   "3210" = "2010-2012-2013",
#   "3211" = "2011-2012-2013"
#   # "3333" = "2010-2011-2012-2013"
# )

# [COLOR BLIND FRIENDLY PALETTE]
# palette_name <- paste("COLOR_BLIND_FRIENDLY_PALETTE.pdf", sep = "")
# pdf(palette_name,  width=15, height=10)
palette <- c( "#5566AA", "#117733", "#44AA66", "#55AA22", "#668822", "#99BB55", "#558877", "#88BBAA",
              "#AADDCC", "#44AA88", "#DDCC66", "#FFDD44", "#FFEE88", "#BB0011", "#000000", "#E69F00",
              "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

for (g in seq(nrow(REF_DF_LIST_ALL))) {
    
    REF_DF <- read.csv(
      file.path(
        "./input_data",
        REF_DF_LIST_ALL$File_Name[g]
      )
    )
    
    df_base <- as.data.frame(REF_DF)
    
    df_base$nRow <- seq_len(nrow(df_base))
    df_base <- df_base %>%
      relocate(nRow, .before = 1)
    
    df_base$Unique_ID <- paste(
      df_base$nRow,
      df_base$SAMPLE_ID,
      sep = "_"
    )
    
    df_base <- df_base %>%
      relocate(Unique_ID, .before = 1)
    
    year_code <- as.character(
      REF_DF_LIST_ALL$Year[g]
    )
    
    df_base$Year <- unname(
      year_label_map[year_code]
    )
    
    df_base <- df_base %>%
      relocate(Year, .before = 1)
    
    df_base$Unique_ID <- paste(
      df_base$Unique_ID,
      df_base$Year,
      sep = "_"
    )
    
    # # Assigning wavelength names into wvl variable 
    df_base <- df_base[
      ,
      seq_len(
        match("X2400", names(df_base))
      ),
      drop = FALSE
    ]
    
    # Check required wavelength endpoints
    if (!all(c("X410", "X2400") %in% names(df_base))) {
      stop(
        "Required wavelength columns X410 and/or X2400 were not found.\n",
        "First spectral-like column names are: ",
        paste(
          head(
            grep("^X?[0-9]+$", names(df_base), value = TRUE),
            20
          ),
          collapse = ", "
        )
      )
    }
    
    start_wvl <- match("X410", names(df_base))
    end_wvl   <- match("X2400", names(df_base))
    
    # Retain metadata, traits and spectra through X2400
    df_base <- df_base[
      ,
      seq_len(end_wvl),
      drop = FALSE
    ]
    
    # Identify wavelength columns directly
    wvl <- names(df_base)[
      start_wvl:end_wvl
    ]
    
    wvl <- wvl[
      grepl("^X[0-9]+$", wvl)
    ]
    
    # Optimal number of iterations per year is 25
    inIter <- 10 #UPDATE 10 to 1000 which is the actual number used for the paper analysis 
    original_inIter <- inIter
    
    inBands <- wvl
    
    # Naming the output folder according to the years
    inYear <- year_label_map[year_code]
    
    # Only traits mapped to this year combination
    if (year_code %in% names(year_trait_map)) {
      columns_of_interest <- year_trait_map[[year_code]]
    } else {
      message(paste("No trait mapping found for year:", year_code, "- Skipping."))
      next
    }
    
    for (col_name in columns_of_interest) {
      
      if (
        col_name %in% names(df_base) &&
        !all(is.na(df_base[[col_name]]))
      ) {
        
        inIter <- original_inIter
        inVar <- col_name
        
        # Start every trait from the untouched input dataframe
        df <- df_base %>%
          filter(
            !is.na(.data[[inVar]])
          )
        
        # ------------------------------------------------------------
        # TRAIT-INTEGRITY CHECK BEFORE MODELLING
        # ------------------------------------------------------------
        
        if (!(inVar %in% names(df))) {
          stop(
            "Trait column not found: ",
            inVar
          )
        }
        
        if (!is.numeric(df[[inVar]])) {
          stop(
            "Trait column is not numeric: ",
            inVar
          )
        }
        
        if (all(is.na(df[[inVar]]))) {
          stop(
            "Trait contains only missing values: ",
            inVar
          )
        }
        
        cat(
          "\n============================================\n",
          "INPUT FILE: ", REF_DF_LIST_ALL$File_Name[g], "\n",
          "YEAR CODE: ", year_code, "\n",
          "MODEL YEAR: ", inYear, "\n",
          "TRAIT: ", inVar, "\n",
          "N: ", nrow(df), "\n",
          "MIN: ", min(df[[inVar]], na.rm = TRUE), "\n",
          "MEAN: ", mean(df[[inVar]], na.rm = TRUE), "\n",
          "MAX: ", max(df[[inVar]], na.rm = TRUE), "\n",
          "============================================\n",
          sep = ""
        )
        
        # it creates a matrix to see which model had what specific ncomps, R2, rmsep etc.
        nComps_iter_df <- data.frame(
          Year = rep(NA_character_, inIter),
          nComps = rep(NA_real_, inIter),
          R2 = rep(NA_real_, inIter),
          RMSEP = rep(NA_real_, inIter),
          PRESS = rep(NA_real_, inIter)
        )
        
        # It creates an empty dataframe to save all coefficients from all inIter models
        coefs_iter_df <- data.frame(matrix(ncol = length(wvl)+1 , nrow = inIter))
        names(coefs_iter_df)[1] <- "Intercept"
        names(coefs_iter_df)[2:(length(wvl)+1)] <- as.numeric(gsub('X', '', wvl))
        
        # It creates an empty dataframe to save all vips from all inIter models
        vips_iter_df <- data.frame(matrix(ncol = length(wvl), nrow = inIter))
        names(vips_iter_df) <- as.numeric(gsub('X', '', wvl))
        
        # It creates an empty dataframe to save all X loadings from all inIter models' selected component
        loads_iter_df <- data.frame(matrix(ncol = length(wvl), nrow = inIter))
        names(loads_iter_df) <- as.numeric(gsub('X', '', wvl))
        
        # Create Output Directory
        dir.create(paste("./r_model_output/", "PPR2_", inYear, "_", Sys.Date(),"_" , inVar, sep = ""), recursive = T) #Single
        pati_main <- paste("./r_model_output/", "PPR2_", inYear, "_", Sys.Date(),"_" , inVar, sep = "")
        pati <- paste(pati_main, "/", "PPR2_", inYear, "_", Sys.Date(),"_" , inVar, sep = "")
        
        ##### Start of Calculations ####
        set.seed(6955866) 
        tmp_env <- new.env(parent = emptyenv())
        for (iter in seq(inIter)){ 
          
          inYear <- year_label_map[year_code]
          
          print(paste0("##################"))
          print(paste0("  Model --->  ", iter ,"  <---"))
          
          # Splitting data into Validation (test) and Calibration (train) data
          data.source <- spectratrait::create_data_split(dataset=df, approach="dplyr", split_seed=(iter+pi*(88033)), prop=0.8)
          cal.data <- data.source$cal_data # CAL calibration data
          val.data <- data.source$val_data # VAL external validation data
          
          # Structuring formula to be applied on all bands
          form <- paste(inBands,collapse="+")
          form <- as.formula(paste(inVar,"~",form))
          
          ##### Part 1: #####
          #Find number of components 
          set.seed(iter+pi*(69558))
          nComps <- if (0.8 * nrow(df) < 25) {
            floor(0.8 * 0.8 * nrow(df))
          } else {
            25
          }
          nCompsi <- if (0.8 * nrow(df) < 25) {
            floor(0.8 * 0.8 * nrow(df))
          } else {
            25
          }
          i.no <- 1 # number of iterations
          
          outMatR2 <- matrix(data=NA,nrow=i.no,ncol=nComps)
          outMatRMSEP <- matrix(data=NA,nrow=i.no,ncol=nComps)
          outMatPRESS <- matrix(data=NA,nrow=i.no,ncol=nComps)
          
          model.set <- plsr(form,data=cal.data,ncomp=nComps,
                            validation="LOO", 
                            method="oscorespls")
          resR2 <- pls::R2(model.set, intercept = F)[[1]]
          
          # Calculate the adjusted R²
          nw <- nrow(cal.data)  # Number of observations
          # pw <- nComps  # Number of components used in the model
          pw <- seq_along(resR2)  # 1..#components (matches resR2)
          
          adjusted_R2 <- 1 - ((1 - resR2) * (nw - 1) / pmax(1, (nw - pw - 1)))
          
          # Print the adjusted R²
          # adjusted_R2
          # Check if all R2 values are negative in this iteration
          if (all(resR2 < 0)) {
            print("All R2 values are negative in this iteration")
            #break  # Move to the next iteration of the outer loop # this will stop further iterations
            next # this will go back to the next iteration
          }else{
            
            outMatR2[1,seq(model.set$validation$ncomp)] <- resR2
            resRMSEP <- as.numeric(RMSEP(model.set,estimate="CV",intercept=F)$val)
            outMatRMSEP[1,seq(model.set$validation$ncomp)] <-resRMSEP
            # Lower values of PRESS indicate better predictive power
            resPRESS <- as.vector(model.set$validation$PRESS)
            outMatPRESS[1,seq(model.set$validation$ncomp)] <-resPRESS

            # summary(model.set)
            # explained variance per each component
            # explvar(model.set)
            
            # Set up the plot area
            ### Here we decide the number of components by plotting average of R2 and PRESS of inIter model iterations.
            # Plot all inIter models' R2 & PRESS
            outMatPRESS_df_t <- as.data.frame(t(outMatPRESS))
            outMatPRESS_df_t$comps <- 1:nrow(outMatPRESS_df_t)
            #outMatR2_df <- as.data.frame(outMatR2)
            outMatR2_df <- as.data.frame(t(outMatR2))
            outMatR2_df$comps <- 1:nrow(outMatR2_df)
            outMatPRESS_df_t <- merge(outMatPRESS_df_t,outMatR2_df, by="comps" )
            first_neg_index <- which(outMatPRESS_df_t$V1.y < 0)[1]
            
            # If there is at least one negative value in V1.y, filter out rows after the first negative value
            if (!is.na(first_neg_index)) {
              # outMatPRESS_df_t <- outMatPRESS_df_t[1:(first_neg_index - 1), ]
              outMatPRESS_df_t <- outMatPRESS_df_t[1:max(1, first_neg_index - 1), ]
            } else {
              # If there are no negative values in V1.y, do nothing
            }
            
            # Keep only components with valid PRESS and positive cross-validated R2
            positive_rows <- outMatPRESS_df_t %>%
              filter(
                is.finite(V1.x),
                is.finite(V1.y),
                V1.y > 0
              )
            
            # Skip the iteration if no eligible component exists
            if (nrow(positive_rows) == 0) {
              
              cat(
                "No positive cross-validated R2 values available in this iteration.\n"
              )
              
              next
            }
            
            # Among eligible components, choose the one with minimum PRESS
            min_V1x_row <- positive_rows %>%
              filter(
                V1.x == min(V1.x)
              )
            
            nComps <- min_V1x_row$comps[1]
            
            max_R2_index <- which.max(outMatR2_df$V1)
            max_comps <- outMatR2_df$comps[max_R2_index]
            #R2_axis_label <- c(min(outMatR2_df$V1)-0.2, max(outMatR2_df$V1)+0.2)
            
            # max_R2 <- round(as.numeric(max(outMatR2_df$V1)),2)
            max_R2 <- round(as.numeric(outMatR2_df$V1[outMatR2_df$comps == nComps]), 4)
            min_R2 <- round(as.numeric(min(outMatR2_df$V1)), 4)
            # min_PRESS <- round(as.numeric(min(outMatPRESS_df_t$V1.x)),2)
            min_PRESS <- round(as.numeric(outMatPRESS_df_t$V1.x[outMatPRESS_df_t$comps == nComps]), 2)
            max_PRESS <- round(as.numeric(max(outMatPRESS_df_t$V1.x)),2)
            
            # save each nComps in R session
            nComps <- nComps
            #assign(paste0("nComps", iter), nComps) #, envir = tmp_env
            
            # Save Some Statistics at each iteration into a dataframe
            nComps_iter_df[iter, "nComps"] <- nComps
            nComps_iter_df[iter, "R2"] <- max_R2
            # nComps_iter_df[iter, "RMSEP"] <- as.numeric(min(outMatRMSEP))
            nComps_iter_df[iter, "RMSEP"] <- as.numeric(outMatRMSEP[1, nComps])
            nComps_iter_df[iter, "PRESS"] <- min_PRESS
            nComps_iter_df[iter, "Valid_Iter"] <- iter
            nComps_iter_df[iter, "Year"] <- inYear
            nComps_iter_df[iter, "Processing"] <- inVar
            
            
            ##### Part 2: #####
            set.seed(6955866)
            
            model.pro <- plsr(form,data=cal.data,ncomp=nComps,
                              validation= "LOO", 
                              method="oscorespls")
            
            # summary(model.pro)
            # explained variance per each component
            # explvar(model.pro)
            
            
            ### Model Fit
            pred_cal_df <- data.frame(model.pro$fitted.values[,,nComps]) #11 #9
            pred_cal_df$measured <- model.pro$model[[inVar]]
            names(pred_cal_df)[1] <- "predicted"
            pred_cal_df <- pred_cal_df
            pred_cal_df$cal_data_Unique_ID <- cal.data$Unique_ID
            pred_cal_df$cal_data_measured <- cal.data[[inVar]]
            assign(paste0("pred_cal_df", iter), pred_cal_df, envir = tmp_env)
            
            
            ### Model applied on external Val Data
            pred_val <- predict(model.pro, newdata = val.data, ncomp = nComps)
            
            # Combine predicted and actual values in a data frame
            pred_val_df <- data.frame(pred_val)
            colnames(pred_val_df)[1] <- "predicted"
            pred_val_df$measured <- val.data[[inVar]]
            pred_val_df <- pred_val_df
            pred_val_df$val_data_Unique_ID <- val.data$Unique_ID
            pred_val_df$val_data_measured <- val.data[[inVar]]
            assign(paste0("pred_val_df", iter), pred_val_df, envir = tmp_env)
            
            
            # save each Coefficients in a dataframe
            coefs <- as.data.frame(t(as.data.frame(coef(model.pro, ncomp=nComps, intercept=TRUE))))
            # save each VIPs in a dataframe
            vips <- as.data.frame(t(as.data.frame(if (is.null(dim(spectratrait::VIP(model.pro)))) spectratrait::VIP(model.pro) else spectratrait::VIP(model.pro)[nComps, ] )))
            
            
            # save the most contributed bands to the nth component at each iteration
            x_loadings <- as.data.frame(model.pro$loadings[, nComps, drop = FALSE])
            names(x_loadings) <- "x_load"
            x_loadings <- as.data.frame(t(x_loadings))
            
          }
            
          # add variables into the dataframe
          coefs_iter_df[iter, ] <- coefs 
          vips_iter_df[iter, ] <- vips 
          loads_iter_df[iter, ] <- x_loadings 
            
        } 

        
        nComps_iter_df <- nComps_iter_df[complete.cases(nComps_iter_df[, c("nComps","R2","RMSEP","PRESS")]), , drop = FALSE]
        coefs_iter_df  <- coefs_iter_df[rownames(nComps_iter_df), , drop = FALSE]
        vips_iter_df   <- vips_iter_df[rownames(nComps_iter_df), , drop = FALSE]
        loads_iter_df  <- loads_iter_df[rownames(nComps_iter_df), , drop = FALSE]
        if (nrow(nComps_iter_df) == 0) next
        
        coefs_xlsx <- cbind(nComps_iter_df, coefs_iter_df)
        vips_xlsx  <- cbind(nComps_iter_df, vips_iter_df)
        loads_xlsx <- cbind(nComps_iter_df, loads_iter_df)
        
        inIterli <- nrow(nComps_iter_df)
        
        # save inIter models' statistics before averaging
        write.csv(coefs_iter_df, paste0(pati_main,"/", inYear, "_", inVar, "_REF_DF_MODELs-", inIterli , "-COEFS.csv"),row.names = FALSE)
        write.csv(vips_iter_df, paste0(pati_main,"/", inYear, "_", inVar, "_REF_DF_MODELs-", inIterli , "-VIPS.csv"),row.names = FALSE)
        write.csv(nComps_iter_df, paste0(pati_main,"/", inYear, "_", inVar, "_REF_DF_MODELs-", inIterli , "-STATISTICS.csv"),row.names = FALSE)
        write.csv(loads_iter_df, paste0(pati_main,"/", inYear, "_", inVar, "_REF_DF_MODELs-", inIterli , "-X_Loadings.csv"),row.names = FALSE)
        
        # added excel stats
        write.xlsx(
          coefs_xlsx,
          file = paste0(pati_main,"/", inYear, "_", inVar, "_REF_DF_MODELs-", inIterli , "-COEFS.xlsx"),
          rowNames = FALSE
        )
        
        write.xlsx(
          vips_xlsx,
          file = paste0(pati_main,"/", inYear, "_", inVar, "_REF_DF_MODELs-", inIterli , "-VIPS.xlsx"),
          rowNames = FALSE
        )
        
        write.xlsx(
          loads_xlsx,
          file = paste0(pati_main,"/", inYear, "_", inVar, "_REF_DF_MODELs-", inIterli , "-X_Loadings.xlsx"),
          rowNames = FALSE
        )
        
        
        
        ##### Part 3: #####
        #### Ensemble PLSR Models 
        
        # Find the most commonly selected nComps among inIter models
        # Create a frequency table of the values in the nComps column
        freq_table <- table(nComps_iter_df$nComps)
        # Find the nComps with the highest frequency
        common_nComp <- as.numeric(names(freq_table)[which.max(freq_table)])
        
        # Step 1: Find data frame names starting with "pred_cal_df"
        # Step 2: Count the number of data frames found after excluding "pred_cal_df" because "pred_cal_df"="pred_cal_df1"
        # Step 3: Check the count and take action accordingly

        if(length(ls(tmp_env, pattern = "^pred_cal_df\\d+$")) < 2){
          cat("Less than 2 data frames starting with 'pred_cal_df'. Breaking the loop or taking other action.\n")
          # You can break a loop or take another action here
        } else {
          cat("2 or more data frames starting with 'pred_cal_df' found.\n")
          # Continue with your loop or other operations
          
          #### THE HOLY MODEL ####
          ### Final Model Fit -> Averages of all predicted and measured of all models'
          
          ### Take average of predicted and measured from cal and cal datasets
          set.seed(6955866)
          
          # This ensembles iterations that were not empty
          ensemble_cal_list <- lapply(ls(tmp_env, pattern = "^pred_cal_df\\d+$"), function(x) get(x, envir = tmp_env))
          ensemble_val_list <- lapply(ls(tmp_env, pattern = "^pred_val_df\\d+$"), function(x) get(x, envir = tmp_env))
          
          ensemble_cal_pre <- do.call(rbind, ensemble_cal_list)
          ensemble_val_pre <- do.call(rbind, ensemble_val_list)
          
          ensemble_cal_preds <- ensemble_cal_pre %>%
            group_by(cal_data_Unique_ID) %>%
            summarize(
              predicted = mean(predicted),
              measured = mean(measured),
              cal_data_measured = mean(cal_data_measured))
          
          ensemble_val_preds <- ensemble_val_pre %>%
            group_by(val_data_Unique_ID) %>%
            summarize(
              predicted = mean(predicted),
              measured = mean(measured),
              val_data_measured = mean(val_data_measured))    
          
          # ------------------------------------------------------------
          # VERIFY THAT SAVED MEASURED VALUES MATCH THE TRAIT INPUT
          # ------------------------------------------------------------
          
          source_values_by_id <- df %>%
            select(
              Unique_ID,
              source_measured = all_of(inVar)
            )
          
          validation_check <- ensemble_val_preds %>%
            left_join(
              source_values_by_id,
              by = c(
                "val_data_Unique_ID" = "Unique_ID"
              )
            )
          
          if (any(is.na(validation_check$source_measured))) {
            stop(
              "Some validation IDs could not be matched back to ",
              "the source dataframe for trait: ",
              inVar
            )
          }
          
          max_measurement_difference <- max(
            abs(
              validation_check$measured -
                validation_check$source_measured
            ),
            na.rm = TRUE
          )
          
          if (
            !is.finite(max_measurement_difference) ||
            max_measurement_difference > 1e-10
          ) {
            stop(
              "Measured-value mismatch detected for ",
              inYear,
              " / ",
              inVar,
              ". Maximum difference = ",
              max_measurement_difference
            )
          }
          
          cat(
            "Measured-value integrity check passed for ",
            inYear,
            " / ",
            inVar,
            "\n",
            sep = ""
          )
          
          calibration_check <- ensemble_cal_preds %>%
            left_join(
              source_values_by_id,
              by = c(
                "cal_data_Unique_ID" = "Unique_ID"
              )
            )
          
          max_calibration_difference <- max(
            abs(
              calibration_check$measured -
                calibration_check$source_measured
            ),
            na.rm = TRUE
          )
          
          if (
            !is.finite(max_calibration_difference) ||
            max_calibration_difference > 1e-10
          ) {
            stop(
              "Calibration measured-value mismatch detected for ",
              inYear,
              " / ",
              inVar,
              ". Maximum difference = ",
              max_calibration_difference
            )
          }
          
          
          # THE HOLY MODEL #
          mod_fit <- lm(predicted~measured, ensemble_cal_preds)
          
          # Print the summary of the linear regression model
          R2_plot_cal <- round(summary(mod_fit)$adj.r.squared, 2)
          
          # Calculate RMSEP for calibration data
          residuals_cal <- ensemble_cal_preds$predicted - ensemble_cal_preds$measured
          mse_cal <- mean(residuals_cal^2)
          rmsep_cal <- round(sqrt(mse_cal),1)
          
          
          # Apply the model on val.data
          # Predict validation data -> take average of predicted and measured
          mod_fit_val <- lm(predicted~measured, data = ensemble_val_preds)
          # Print the summary of the linear regression model
          # summary(mod_fit_val)
          # summary(mod_fit_val)$adj.r.squared
          R2_plot_val <- round(summary(mod_fit_val)$adj.r.squared, 2)
          
          # Calculate RMSEP for validation data
          residuals_val <- ensemble_val_preds$predicted - ensemble_val_preds$measured
          mse_val <- mean(residuals_val^2)
          rmsep_val <- round(sqrt(mse_val),1)
          
          ########################################################
          
          
          ##### Part 4: #####
          # Calculated out of inIter models' average
          # Avg_Coefficents
          coefs_iter_avg_df <- as.data.frame(colMeans(na.omit(coefs_iter_df)))
          names(coefs_iter_avg_df)[1] <- "Avg_Coefficents"
          coefs_iter_avg_df$Wavelength <- rownames(coefs_iter_avg_df) #there is intercept here, so it can't be as.numeric
          
          # Avg_VIPs
          vips_iter_avg_df <- as.data.frame(colMeans(na.omit(vips_iter_df))) 
          names(vips_iter_avg_df)[1] <- "Avg_VIPs"
          vips_iter_avg_df$Wavelength <- as.numeric(rownames(vips_iter_avg_df))
          
          # Avg_X_Loadings
          loads_iter_avg_df <- as.data.frame(colMeans(na.omit(loads_iter_df)))  
          names(loads_iter_avg_df)[1] <- "Avg_X_Loadings"
          loads_iter_avg_df$Wavelength <- as.numeric(rownames(loads_iter_avg_df))
          
          # Avg_Models_Statistics
          #na.omit(nComps_iter_df)
          nComps_iter_avg_df <- colMeans(na.omit(nComps_iter_df[ ,3:5]))
          nComps_iter_avg_df <- as.data.frame(t(nComps_iter_avg_df))
          nComps_iter_avg_df$Year <- inYear
          nComps_iter_avg_df$inVar <- inVar
          nComps_iter_avg_df$Valid_Iter <- inIterli
          nComps_iter_avg_df$nComps <- round(mean(nComps_iter_df[complete.cases(nComps_iter_df), 2]),0)
          nComps_iter_avg_df$PRESS <- round(nComps_iter_avg_df$PRESS,4)
          nComps_iter_avg_df$RMSEP <- round(nComps_iter_avg_df$RMSEP,4)
          nComps_iter_avg_df$R2 <- round(nComps_iter_avg_df$R2,4)
          nComps_iter_avg_df <- nComps_iter_avg_df[ ,c("Year", "inVar", "Valid_Iter", "nComps", "PRESS", "RMSEP", "R2")]
          nComps_iter_avg_df$Cal_R2 <- round(summary(mod_fit)$adj.r.squared, 4)
          nComps_iter_avg_df$Val_R2 <- round(summary(mod_fit_val)$adj.r.squared, 4)
          nComps_iter_avg_df$Cal_RMSEP <- rmsep_cal
          nComps_iter_avg_df$Val_RMSEP <- rmsep_val
          nComps_iter_avg_df
          
          
          ##### Part 5: #####
          # The Holy model results as graphs etc.
          
          # Define the file path where variables will be saved
          save_file_path <- paste0(pati_main, "/", inYear, "_", inVar, "_", Sys.Date(), "_plot_data.RData")
          # Save all relevant variables to an RData file
          save(ensemble_cal_preds, mod_fit, R2_plot_cal, rmsep_cal,
               ensemble_val_preds, mod_fit_val, R2_plot_val, rmsep_val, pati_main, 
               inYear, inVar, file = save_file_path)
          
          # Plot predicted vs actual values
          png(file = paste0(pati_main,"/", inYear, "_", inVar, "_", "cal_prediction.png"), width = 1250, height = 850, res = 300)
          par(mar = c(3, 3, 0.4, 0.25)+0.1)#par(mar = c(bottom, left, top, right) + additional)
          plot(predicted ~ measured, ensemble_cal_preds, pch=15, cex=1, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
               xlim=c(0,2500), ylim=c(0,2000)) 
          axis(1, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
          axis(2, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
          mtext(expression("Measured [g" ~ m^{-2} ~ "]"), side=1, line=2.1, cex=1.25, font=2)
          mtext(expression("Predicted [g" ~ m^{-2} ~ "]"), side=2, line=1.5, cex=1.25, font=2)
          abline(mod_fit)
          # add expressions to top left corner
          mtext(paste(" R² =", sprintf(as.character(R2_plot_cal))), side = 3, line = -1.6, at = par("usr")[1], adj = 0, cex = 1.8)
          mtext(paste(" RMSEP =", sprintf(as.character(rmsep_cal))), side = 3, line = -3, at = par("usr")[1], adj = 0, cex = 1.8)
          #mtext(paste(" nComps =", sprintf(as.character(common_nComp))), side = 3, line = -3, at = par("usr")[1], adj = 0, cex = 1.8)
          graphics.off()
          
          # Plot predicted vs actual values
          png(file = paste0(pati_main,"/", inYear, "_", inVar, "_", "val_prediction.png"), width = 1250, height = 850, res = 300)
          par(mar = c(3, 3, 0.4, 0.25)+0.1)#par(mar = c(bottom, left, top, right) + additional)
          plot(predicted ~ measured, data = ensemble_val_preds, pch=15, cex=1, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
               xlim=c(0,2500), ylim=c(0,2000))  
          axis(1, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
          axis(2, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
          mtext(expression("Measured [g" ~ m^{-2} ~ "]"), side=1, line=2.1, cex=1.25, font=2)
          mtext(expression("Predicted [g" ~ m^{-2} ~ "]"), side=2, line=1.5, cex=1.25, font=2)
          abline(mod_fit_val)
          # add expressions to top left corner
          mtext(paste(" R² =", sprintf(as.character(R2_plot_val))), side = 3, line = -1.6, at = par("usr")[1], adj = 0, cex = 1.8)
          mtext(paste(" RMSEP =", sprintf(as.character(rmsep_val))), side = 3, line = -3, at = par("usr")[1], adj = 0, cex = 1.8)
          #mtext(paste(" nComps =", sprintf(as.character(common_nComp))), side = 3, line = -3, at = par("usr")[1], adj = 0, cex = 1.8)
          graphics.off()
          
          
          ### Plot
          # get the x value of the highest point
          
          plot_x_loadings <- as.data.frame(t(loads_iter_avg_df))
          plot_x_loadings <- plot_x_loadings[1,]
          plot_x_loadings <- melt(plot_x_loadings)
          
          top_15_positive <- plot_x_loadings %>%
            top_n(n = 30, wt = value)
          
          top_15_negative <- plot_x_loadings %>%
            arrange(desc(value)) %>%
            tail(n = 30)
          
          
          # Save positively and negatively correlated bands of the final model
          # which is the average of 25 models' x loadings. All top bands were selected as important
          # 25 times in different models
          Top_bands_pos <- as.data.frame(top_15_positive)
          names(Top_bands_pos) <- c("Wavelength", "top_30_positive")
          Top_bands_pos$Year <- inYear
          Top_bands_pos$inVar <- inVar
          Top_bands_pos$Model <- "Average"
          
          Top_bands_neg <- as.data.frame(top_15_negative)
          names(Top_bands_neg) <- c("Wavelength", "top_30_negative")
          Top_bands_neg$Year <- inYear
          Top_bands_neg$inVar <- inVar
          Top_bands_neg$Model <- "Average"
          
          # Define the file path where variables will be saved
          save_file_path <- paste0(pati_main, "/", inYear, "_", Sys.Date(),"_" , inVar, "_plot_wvl_contrib_data.RData")
          # Save all relevant variables to an RData file
          save(top_15_positive, top_15_negative, plot_x_loadings, inVar, pati_main, inYear, file = save_file_path)
          
          wvl_contrib <- ggplot() +
            geom_vline(data = top_15_positive, aes(xintercept = as.numeric(variable)), color = "#AADDCC", linetype = "solid", alpha = 0.85, size = 0.9) +
            geom_vline(data = top_15_negative, aes(xintercept = as.numeric(variable)), color = "#E69F00", linetype = "solid", alpha = 0.7, size = 0.9) +
            geom_line(data = plot_x_loadings, aes(x = factor(variable), y = `value`, group = inVar, color = inVar), size = 2.5, alpha = 0.9) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "#BB0011", size=1.5, alpha = 0.8) + # add a dashed red line at 0 on the y-axis
            labs( #subtitle = paste0(inYear," ", inVar, " ", inVar," PLSR models' top 30 bands"),
              x = "Wavelength [nm]", y = "Loadings") +
            scale_y_continuous(breaks = seq(round(min(top_15_negative$value),2), round(max(top_15_positive$value),2), 0.02), 
                               limits = c(min(top_15_negative$value),max(top_15_positive$value))) +
            scale_x_discrete(breaks = levels(factor(plot_x_loadings$variable))[seq(1, length(levels(factor(plot_x_loadings$variable))), 25)]
                             ,labels = c(410,600,680,775,920,1175,1415,1650,1860,2060,2235,2400)) +
            theme(axis.text.y = element_text(hjust = 1, size=18, face = "bold"),
                  axis.text.x = element_text(angle = 45,hjust = 1, size=18, face = "bold"),
                  #plot.subtitle=element_text(size=14,face="bold", color="black", hjust = 0.5),
                  axis.title.x = element_text(size = 24, face = "bold"),
                  axis.title.y = element_text(size = 24, face = "bold"),
                  legend.position = "top", # move the legend to the top
                  #legend.title = element_text(size = 18),
                  #legend.text = element_text(size = 18),
                  legend.margin = margin(t = 0, unit = "cm"), # adjust the top margin of the legend
                  plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")) +
            theme(panel.background = element_blank(),
                  #panel.grid.major = element_line(color = "gray"),
                  #panel.grid.minor = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 1),
                  plot.margin = unit(c(0.15, 0.25, 0.15, 0.15), "cm"))+ # top, right, bottom, left
            scale_color_manual(name = "", #VIPs:
                               values = c( inVar = "black"))# +
          ggsave(filename = file.path(paste0(pati_main,"/", inYear, "_", inVar, "_", "common-comps-", "REF_DF_MODEL_X_Loadings.pdf")), plot=wvl_contrib,
                 device="pdf", width = 24, height = 12, units = "cm", dpi = 300)
          ggsave(filename = file.path(paste0(pati_main,"/", inYear, "_", inVar, "_", "common-comps-", "REF_DF_MODEL_X_Loadings.png")), plot=wvl_contrib,
                 device="png", width = 24, height = 12, units = "cm", dpi = 300)
          wvl_contrib
          graphics.off()
          ###
          
          nComps_iter_avg_df$Valid_Iter <- inIterli
          
          ##### Save Coefficients, VIPs, X_Loadings #####
          # save averaged model's statistics
          write.csv(coefs_iter_avg_df, paste0(pati_main,"/", inYear, "_", inVar, "_", "MAIN_MODELs-avg-COEFS.csv"))
          write.csv(vips_iter_avg_df, paste0(pati_main,"/", inYear, "_", inVar, "_", "MAIN_MODELs-avg-VIPS.csv"))
          write.csv(nComps_iter_avg_df, paste0(pati_main,"/", inYear, "_", inVar, "_", "MAIN_MODELs-avg-STATISTICS.csv"))
          write.csv(loads_iter_avg_df, paste0(pati_main,"/", inYear, "_", inVar, "_", "MAIN_MODEL-avg-X_Loadings.csv"))
          write.csv(Top_bands_pos, paste0(pati_main,"/", inYear, "_", inVar, "_", "MAIN_MODEL-avg-top_30_POS_bands.csv"))
          write.csv(Top_bands_neg, paste0(pati_main,"/", inYear, "_", inVar, "_", "MAIN_MODEL-avg-top_30_NEG_bands.csv"))
          
          rm(list = ls(tmp_env), envir = tmp_env)
          gc()

          
        }
        
        
      } else {
        print(paste("Column not found:", col_name))
        # Handle case where column is not found
        # For example, continue to next iteration of columns_of_interest
        next  # This will continue to the next iteration of 'col_name'
      }
      

      
    }
    

    
    # --- CLEANUP MEMORY AFTER EACH TRAIT ---
    rm(list = ls(pattern = "^cal\\.data\\d+$"))
    rm(list = ls(pattern = "^val\\.data\\d+$"))
    rm(list = ls(pattern = "^model\\.pro\\d+$"))
    rm(list = ls(pattern = "^pred_cal_df\\d+$"))
    rm(list = ls(pattern = "^pred_val_df\\d+$"))
    gc()
    # ----------------------------------------
    
        
}

# View(cbind(cal.data1[,c("SAMPLE_ID", inVar)],pred_cal_df1))
# View(cbind(cal.data1[,c("SAMPLE_ID", inVar)],ensemble_cal_preds))
# View(cbind(pred_cal_df1,ensemble_cal_preds))

# View(cbind(val.data1[,c("SAMPLE_ID", inVar)],pred_val_df1))
# View(cbind(val.data1[,c("SAMPLE_ID", inVar)],ensemble_val_preds))
# View(cbind(pred_val_df1,ensemble_val_preds))



### END ###


