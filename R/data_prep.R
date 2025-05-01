#' Load clinical trial data from a CSV file
#' 
#' @param file_path Path to CSV file
#' @param data_format Format of the data (default, sdtm, adam)
#' @param ... Additional arguments passed to read.csv
#' @return A data frame containing clinical trial data
#' @export
load_clinical_trial <- function(file_path, data_format = c("default", "sdtm", "adam"), ...) {
  data_format <- match.arg(data_format)
  
  # Check if file exists
  check_file(file_path, "csv")
  
  # Read the data
  data <- utils::read.csv(file_path, stringsAsFactors = FALSE, ...)
  
  # Process based on format
  if (data_format == "sdtm") {
    data <- process_sdtm_data(data)
  } else if (data_format == "adam") {
    data <- process_adam_data(data)
  }
  
  # Add class
  class(data) <- c("bayestrials_data", class(data))
  
  return(data)
}

#' Process SDTM format data
#' 
#' @param data Data frame in SDTM format
#' @return Processed data frame
#' @keywords internal
process_sdtm_data <- function(data) {
  # Check for required SDTM variables
  required_vars <- c("USUBJID", "ARMCD", "ARM")
  
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required SDTM variables:", paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  
  # Add standardized variable names as attributes
  attr(data, "id_var") <- "USUBJID"
  attr(data, "treatment_var") <- "ARM"
  attr(data, "treatment_code_var") <- "ARMCD"
  
  # Look for common outcome variables
  outcome_candidates <- c("AVAL", "AVALC", "CNSR", "DTHFL")
  found_outcomes <- intersect(outcome_candidates, names(data))
  
  if (length(found_outcomes) > 0) {
    attr(data, "outcome_vars") <- found_outcomes
  }
  
  # Look for common demographic variables
  demo_candidates <- c("AGE", "SEX", "RACE", "ETHNIC", "COUNTRY", "SITEID")
  found_demos <- intersect(demo_candidates, names(data))
  
  if (length(found_demos) > 0) {
    attr(data, "demographic_vars") <- found_demos
  }
  
  return(data)
}

#' Process ADaM format data
#' 
#' @param data Data frame in ADaM format
#' @return Processed data frame
#' @keywords internal
process_adam_data <- function(data) {
  # Check for required ADaM variables
  required_vars <- c("USUBJID", "TRT01P", "TRT01PN")
  
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required ADaM variables:", paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  
  # Add standardized variable names as attributes
  attr(data, "id_var") <- "USUBJID"
  attr(data, "treatment_var") <- "TRT01P"
  attr(data, "treatment_code_var") <- "TRT01PN"
  
  # Look for common outcome variables
  outcome_candidates <- c("AVAL", "AVALC", "CHG", "PCHG", "CNS", "EVNT", "CNSR")
  found_outcomes <- intersect(outcome_candidates, names(data))
  
  if (length(found_outcomes) > 0) {
    attr(data, "outcome_vars") <- found_outcomes
  }
  
  # Look for common demographic variables
  demo_candidates <- c("AGE", "SEX", "RACE", "ETHNIC", "COUNTRY", "SITEID")
  found_demos <- intersect(demo_candidates, names(data))
  
  if (length(found_demos) > 0) {
    attr(data, "demographic_vars") <- found_demos
  }
  
  return(data)
}

#' Load multiple clinical trials
#' 
#' @param file_paths Vector of paths to CSV files
#' @param trial_names Vector of names for the trials (optional)
#' @param data_format Format of the data (default, sdtm, adam)
#' @param ... Additional arguments passed to read.csv
#' @return A list of data frames containing clinical trial data
#' @export
load_trials <- function(file_paths, trial_names = NULL, data_format = c("default", "sdtm", "adam"), ...) {
  data_format <- match.arg(data_format)
  
  # Check if all files exist
  sapply(file_paths, check_file, extension = "csv")
  
  # Generate trial names if not provided
  if (is.null(trial_names)) {
    trial_names <- paste0("trial", seq_along(file_paths))
  } else if (length(trial_names) != length(file_paths)) {
    stop("Length of trial_names must match length of file_paths", call. = FALSE)
  }
  
  # Load each trial
  trials <- lapply(file_paths, load_clinical_trial, data_format = data_format, ...)
  
  # Name the trials
  names(trials) <- trial_names
  
  # Add a class to the list
  class(trials) <- c("bayestrials_trials", class(trials))
  
  return(trials)
}

#' Harmonize variables across trials
#' 
#' @param trial_list List of trial data frames
#' @param mapping_file Path to variable mapping CSV file (optional)
#' @param auto_detect Logical indicating whether to attempt automatic variable detection
#' @param standard Data standard to use (none, sdtm, adam)
#' @return List of harmonized trial data frames
#' @export
harmonize_trials <- function(trial_list, 
                            mapping_file = NULL, 
                            auto_detect = TRUE,
                            standard = c("none", "sdtm", "adam")) {
  
  standard <- match.arg(standard)
  
  # Check if trial_list is valid
  if (!is.list(trial_list) || length(trial_list) == 0) {
    stop("trial_list must be a non-empty list of data frames", call. = FALSE)
  }
  
  # Read mapping file if provided
  mapping <- NULL
  if (!is.null(mapping_file)) {
    check_file(mapping_file, "csv")
    mapping <- utils::read.csv(mapping_file, stringsAsFactors = FALSE)
    
    # Check required column
    if (!"standard_name" %in% names(mapping)) {
      stop("Mapping file must contain a 'standard_name' column", call. = FALSE)
    }
  }
  
  # Auto-detect mappings if requested
  if (auto_detect) {
    mapping <- auto_detect_variables(trial_list, mapping, standard)
  }
  
  # If no mapping after auto-detection, create minimal mapping
  if (is.null(mapping)) {
    stop("No variable mapping available. Please provide a mapping file or enable auto_detect", call. = FALSE)
  }
  
  # Apply the mappings to each trial
  harmonized_trials <- lapply(seq_along(trial_list), function(i) {
    trial_name <- names(trial_list)[i]
    if (is.null(trial_name)) trial_name <- paste0("trial", i)
    
    trial_data <- trial_list[[i]]
    
    # Get the column name for this trial
    trial_col <- paste0(trial_name, "_name")
    
    # If trial column doesn't exist in mapping, skip this trial
    if (!trial_col %in% names(mapping)) {
      warning(paste("No mapping column for trial", trial_name, "- skipping harmonization"), call. = FALSE)
      return(trial_data)
    }
    
    # Create a named vector of mappings for this trial
    trial_mapping <- setNames(
      mapping$standard_name,
      mapping[[trial_col]]
    )
    
    # Remove NA mappings
    trial_mapping <- trial_mapping[!is.na(names(trial_mapping)) & names(trial_mapping) != ""]
    
    # Apply the mapping to the trial
    harmonized_data <- rename_variables(trial_data, trial_mapping)
    
    # Add attributes
    attr(harmonized_data, "trial_name") <- trial_name
    attr(harmonized_data, "original_names") <- names(trial_list[[i]])
    attr(harmonized_data, "mapping") <- trial_mapping
    
    return(harmonized_data)
  })
  
  # Preserve the original names
  names(harmonized_trials) <- names(trial_list)
  
  # Add full mapping as attribute
  attr(harmonized_trials, "full_mapping") <- mapping
  
  # Add harmonized class
  class(harmonized_trials) <- c("harmonized_trials", class(harmonized_trials))
  
  return(harmonized_trials)
}

#' Auto-detect variables across trials
#' 
#' @param trial_list List of trial data frames
#' @param existing_mapping Existing mapping data frame (optional)
#' @param standard Data standard to use (none, sdtm, adam)
#' @return Data frame with variable mappings
#' @keywords internal
auto_detect_variables <- function(trial_list, existing_mapping = NULL, standard = "none") {
  trial_names <- names(trial_list)
  if (is.null(trial_names)) {
    trial_names <- paste0("trial", seq_along(trial_list))
  }
  
  # Initialize mapping if not provided
  if (is.null(existing_mapping)) {
    existing_mapping <- data.frame(standard_name = character(0), stringsAsFactors = FALSE)
  }
  
  # Add columns for each trial if they don't exist
  for (trial_name in trial_names) {
    trial_col <- paste0(trial_name, "_name")
    if (!trial_col %in% names(existing_mapping)) {
      existing_mapping[[trial_col]] <- NA_character_
    }
  }
  
  # Common variable patterns to detect
  variable_patterns <- list(
    # ID variables
    id = c("id", "subject", "usubjid", "patid", "patnum", "ptid", "subjid"),
    
    # Treatment variables
    treatment = c("treatment", "arm", "trt", "trtp", "trt01p", "trtgrp", "armcd", "group"),
    
    # Outcome variables
    outcome = c("outcome", "response", "resp", "aval", "avalc", "cnsr", "event", "death", "dthfl", "chg"),
    
    # Demographic variables
    age = c("age", "age_years", "pt_age", "ageyears"),
    sex = c("sex", "gender"),
    race = c("race", "ethnic", "ethnicity"),
    
    # Site variables
    site = c("site", "center", "siteid", "clinic", "facility"),
    
    # Time variables
    time = c("time", "day", "week", "month", "visit", "visnum", "visitdy")
  )
  
  # Apply standard-specific patterns
  if (standard == "sdtm") {
    # Add SDTM-specific variable patterns
    sdtm_patterns <- list(
      id = "USUBJID",
      treatment = c("ARM", "ARMCD"),
      outcome = c("AVAL", "AVALC", "CNSR", "DTHFL"),
      age = "AGE",
      sex = "SEX",
      race = c("RACE", "ETHNIC"),
      site = "SITEID"
    )
    variable_patterns <- c(variable_patterns, sdtm_patterns)
  } else if (standard == "adam") {
    # Add ADaM-specific variable patterns
    adam_patterns <- list(
      id = "USUBJID",
      treatment = c("TRT01P", "TRT01PN"),
      outcome = c("AVAL", "AVALC", "CHG", "PCHG", "CNS", "EVNT", "CNSR"),
      age = "AGE",
      sex = "SEX",
      race = c("RACE", "ETHNIC"),
      site = "SITEID"
    )
    variable_patterns <- c(variable_patterns, adam_patterns)
  }
  
  # Find variables in each trial that match patterns
  for (std_var in names(variable_patterns)) {
    patterns <- variable_patterns[[std_var]]
    
    # If standard variable not in mapping, add it
    if (!std_var %in% existing_mapping$standard_name) {
      # Create a new row with the same columns as existing_mapping
      new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(existing_mapping)))
      names(new_row) <- names(existing_mapping)
      new_row$standard_name <- std_var
      
      # Add the new row
      existing_mapping <- rbind(existing_mapping, new_row)
    }
    
    # For each trial, find variables matching patterns
    for (i in seq_along(trial_list)) {
      trial_name <- trial_names[i]
      trial_col <- paste0(trial_name, "_name")
      
      # Get variable names in the trial
      var_names <- names(trial_list[[i]])
      
      # Convert all names to lowercase for matching
      var_names_lower <- tolower(var_names)
      
      # Match variables with patterns
      for (pattern in patterns) {
        pattern_lower <- tolower(pattern)
        matches <- var_names[var_names_lower == pattern_lower | 
                           grepl(paste0("^", pattern_lower, "$"), var_names_lower)]
        
        if (length(matches) > 0) {
          # Get row index for standard variable
          std_idx <- which(existing_mapping$standard_name == std_var)
          
          # If there's already a mapping for this trial, skip
          if (!is.na(existing_mapping[std_idx, trial_col])) {
            next
          }
          
          # Set the mapping
          existing_mapping[std_idx, trial_col] <- matches[1]
          break  # Stop after first match
        }
      }
    }
  }
  
  return(existing_mapping)
}

#' Rename variables in a data frame
#' 
#' @param data Data frame to rename variables in
#' @param mapping Named vector of new names, where names are old variable names
#' @return Data frame with renamed variables
#' @keywords internal
rename_variables <- function(data, mapping) {
  # Get valid mappings (old names that actually exist in the data)
  valid_mappings <- mapping[names(mapping) %in% names(data)]
  
  if (length(valid_mappings) == 0) {
    warning("No valid mappings found for this data frame", call. = FALSE)
    return(data)
  }
  
  # Create renamed data frame
  renamed_data <- data
  
  # Rename columns
  for (old_name in names(valid_mappings)) {
    new_name <- valid_mappings[old_name]
    
    # Skip if new name already exists and isn't the current column
    if (new_name %in% names(renamed_data) && !new_name == old_name) {
      warning(paste0("Column '", new_name, "' already exists in data frame. ", 
                    "Skipping renaming of '", old_name, "'"), call. = FALSE)
      next
    }
    
    # Rename the column
    names(renamed_data)[names(renamed_data) == old_name] <- new_name
  }
  
  return(renamed_data)
}

#' Preprocess data based on specifications
#' 
#' @param data Data frame to preprocess
#' @param specs Data frame or path to CSV file with covariate specifications
#' @param steps Character vector of preprocessing steps to perform
#' @return Preprocessed data frame
#' @export
preprocess_data <- function(data, 
                           specs = NULL, 
                           steps = c("transform", "recode", "convert", "center", "scale")) {
  
  # Check input data
  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }
  
  # Load specs if file path
  if (is.character(specs) && length(specs) == 1) {
    check_file(specs, "csv")
    specs <- utils::read.csv(specs, stringsAsFactors = FALSE)
  }
  
  # If no specs, return data as is
  if (is.null(specs)) {
    warning("No covariate specifications provided, returning data as is", call. = FALSE)
    return(data)
  }
  
  # Check specs format
  required_cols <- c("variable", "type")
  missing_cols <- setdiff(required_cols, names(specs))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in specs:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  
  # Prep the output data frame
  processed <- data
  
  # Track changes made
  changes <- list()
  
  # Process each variable
  for (i in 1:nrow(specs)) {
    var_name <- specs$variable[i]
    var_type <- specs$type[i]
    
    # Skip if variable not in data
    if (!var_name %in% names(processed)) {
      warning(paste("Variable", var_name, "not found in data, skipping"), call. = FALSE)
      next
    }
    
    var_changes <- list()
    
    # Apply transformations if requested
    if ("transform" %in% steps && "transformation" %in% names(specs) && !is.na(specs$transformation[i])) {
      transformation <- specs$transformation[i]
      
      if (transformation == "log") {
        # Check for non-positive values
        if (any(processed[[var_name]] <= 0, na.rm = TRUE)) {
          warning(paste("Variable", var_name, "contains non-positive values, cannot apply log transformation"),
                  call. = FALSE)
        } else {
          processed[[var_name]] <- log(processed[[var_name]])
          var_changes$transformation <- "log"
        }
      } else if (transformation == "logit") {
        # Check for values outside (0,1)
        if (any(processed[[var_name]] <= 0 | processed[[var_name]] >= 1, na.rm = TRUE)) {
          warning(paste("Variable", var_name, "contains values outside (0,1), cannot apply logit transformation"),
                  call. = FALSE)
        } else {
          processed[[var_name]] <- log(processed[[var_name]] / (1 - processed[[var_name]]))
          var_changes$transformation <- "logit"
        }
      } else if (transformation == "sqrt") {
        # Check for negative values
        if (any(processed[[var_name]] < 0, na.rm = TRUE)) {
          warning(paste("Variable", var_name, "contains negative values, cannot apply sqrt transformation"),
                  call. = FALSE)
        } else {
          processed[[var_name]] <- sqrt(processed[[var_name]])
          var_changes$transformation <- "sqrt"
        }
      } else if (transformation == "scale") {
        # Scale the variable
        if (!is.numeric(processed[[var_name]])) {
          warning(paste("Variable", var_name, "is not numeric, cannot apply scale transformation"),
                  call. = FALSE)
        } else {
          var_mean <- mean(processed[[var_name]], na.rm = TRUE)
          var_sd <- sd(processed[[var_name]], na.rm = TRUE)
          processed[[var_name]] <- (processed[[var_name]] - var_mean) / var_sd
          var_changes$transformation <- "scale"
          var_changes$mean <- var_mean
          var_changes$sd <- var_sd
        }
      }
    }
    
    # Convert variable types
    if ("convert" %in% steps) {
      if (var_type == "factor" && !is.factor(processed[[var_name]])) {
        processed[[var_name]] <- as.factor(processed[[var_name]])
        var_changes$conversion <- "to factor"
      } else if (var_type == "numeric" && !is.numeric(processed[[var_name]])) {
        # Only convert if it makes sense
        if (is.character(processed[[var_name]]) || is.factor(processed[[var_name]])) {
          # Try to convert to numeric
          new_var <- suppressWarnings(as.numeric(as.character(processed[[var_name]])))
          if (sum(is.na(new_var)) > sum(is.na(processed[[var_name]]))) {
            warning(paste("Converting", var_name, "to numeric resulted in NAs, skipping conversion"), call. = FALSE)
          } else {
            processed[[var_name]] <- new_var
            var_changes$conversion <- "to numeric"
          }
        }
      } else if (var_type == "character" && !is.character(processed[[var_name]])) {
        processed[[var_name]] <- as.character(processed[[var_name]])
        var_changes$conversion <- "to character"
      } else if (var_type == "logical" && !is.logical(processed[[var_name]])) {
        # Try to convert to logical
        if (all(unique(na.omit(processed[[var_name]])) %in% c(0, 1, "0", "1", "TRUE", "FALSE", "true", "false", "T", "F"))) {
          processed[[var_name]] <- as.logical(as.character(processed[[var_name]]))
          var_changes$conversion <- "to logical"
        } else {
          warning(paste("Variable", var_name, "contains values other than 0/1/TRUE/FALSE, cannot convert to logical"),
                  call. = FALSE)
        }
      }
    }
    
    # Recode factor variables
    if ("recode" %in% steps && var_type == "factor" && 
        "coding" %in% names(specs) && !is.na(specs$coding[i])) {
      
      coding <- specs$coding[i]
      reference <- if ("reference" %in% names(specs)) specs$reference[i] else NA
      
      # First ensure the variable is a factor
      if (!is.factor(processed[[var_name]])) {
        processed[[var_name]] <- as.factor(processed[[var_name]])
      }
      
      if (coding == "treatment") {
        # Set reference level if provided
        if (!is.na(reference) && reference %in% levels(processed[[var_name]])) {
          processed[[var_name]] <- relevel(processed[[var_name]], ref = reference)
          var_changes$reference <- reference
        }
        
        # Set contrasts but don't replace the variable with the contrast matrix
        contrasts(processed[[var_name]]) <- contr.treatment(nlevels(processed[[var_name]]))
        var_changes$coding <- "treatment contrasts"
      } else if (coding == "sum") {
        contrasts(processed[[var_name]]) <- contr.sum(nlevels(processed[[var_name]]))
        var_changes$coding <- "sum contrasts"
      } else if (coding == "helmert") {
        contrasts(processed[[var_name]]) <- contr.helmert(nlevels(processed[[var_name]]))
        var_changes$coding <- "helmert contrasts"
      } else if (coding == "dummy") {
        # Create dummy variables
        dummies <- model.matrix(~ processed[[var_name]] - 1)
        # Use proper names
        colnames(dummies) <- paste0(var_name, "_", levels(processed[[var_name]]))
        # Convert to data frame before binding
        dummies <- as.data.frame(dummies)
        processed <- cbind(processed, dummies)
        var_changes$coding <- "dummy variables"
      }
    }
    
    # Center numeric variables
    if ("center" %in% steps && var_type == "numeric" && is.numeric(processed[[var_name]])) {
      if ("transformation" %in% names(specs) && specs$transformation[i] == "scale") {
        # Already centered and scaled
        next
      }
      
      var_mean <- mean(processed[[var_name]], na.rm = TRUE)
      processed[[var_name]] <- processed[[var_name]] - var_mean
      var_changes$centering <- "mean centered"
      var_changes$mean <- var_mean
    }
    
    # Scale numeric variables
    if ("scale" %in% steps && var_type == "numeric" && is.numeric(processed[[var_name]])) {
      if ("transformation" %in% names(specs) && specs$transformation[i] == "scale") {
        # Already centered and scaled
        next
      }
      
      var_sd <- sd(processed[[var_name]], na.rm = TRUE)
      processed[[var_name]] <- processed[[var_name]] / var_sd
      var_changes$scaling <- "sd scaled"
      var_changes$sd <- var_sd
    }
    
    # Store changes for this variable
    if (length(var_changes) > 0) {
      changes[[var_name]] <- var_changes
    }
  }
  
  # Store preprocessing info as attribute
  attr(processed, "preprocessing") <- list(
    specs = specs,
    steps = steps,
    changes = changes
  )
  
  # Add class
  class(processed) <- c("preprocessed_data", class(processed))
  
  return(processed)
}

#' Create standardized trial dataset
#' 
#' @param data Data frame with trial data
#' @param id_var Name of ID variable
#' @param treatment_var Name of treatment variable
#' @param outcome_var Name of outcome variable
#' @param covariates Vector of covariate variable names
#' @param time_var Name of time variable (optional)
#' @return Standardized trial dataset
#' @export
create_trial_dataset <- function(data, 
                               id_var,
                               treatment_var,
                               outcome_var,
                               covariates = NULL,
                               time_var = NULL) {
  
  # Check input data
  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }
  
  # Check required variables exist
  required_vars <- c(id_var, treatment_var, outcome_var)
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste("Missing required variables:", paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  
  # Check covariates if provided
  if (!is.null(covariates)) {
    missing_covs <- setdiff(covariates, names(data))
    if (length(missing_covs) > 0) {
      warning(paste("Missing covariates:", paste(missing_covs, collapse = ", "), 
                   "- they will be excluded"), call. = FALSE)
      covariates <- setdiff(covariates, missing_covs)
    }
  }
  
  # Check time variable if provided
  if (!is.null(time_var) && !time_var %in% names(data)) {
    warning(paste("Time variable", time_var, "not found in data - ignoring"), call. = FALSE)
    time_var <- NULL
  }
  
  # Select variables to include
  vars_to_include <- c(id_var, treatment_var, outcome_var, covariates)
  if (!is.null(time_var)) {
    vars_to_include <- c(vars_to_include, time_var)
  }
  
  # Create standardized dataset
  std_data <- data[, vars_to_include, drop = FALSE]
  
  # Rename variables to standard names
  var_mapping <- c(
    id = id_var,
    treatment = treatment_var,
    outcome = outcome_var
  )
  
  if (!is.null(time_var)) {
    var_mapping <- c(var_mapping, time = time_var)
  }
  
  # Create reverse mapping (old -> new)
  reverse_mapping <- setNames(names(var_mapping), var_mapping)
  
  # Apply renaming
  std_data <- rename_variables(std_data, reverse_mapping)
  
  # Store original variable names
  attr(std_data, "original_names") <- setNames(vars_to_include, 
                                              ifelse(vars_to_include %in% names(reverse_mapping),
                                                    reverse_mapping[vars_to_include],
                                                    vars_to_include))
  
  # Store standardized variable info
  attr(std_data, "id_var") <- "id"
  attr(std_data, "treatment_var") <- "treatment"
  attr(std_data, "outcome_var") <- "outcome"
  attr(std_data, "covariate_vars") <- setdiff(names(std_data), c("id", "treatment", "outcome", "time"))
  
  if (!is.null(time_var)) {
    attr(std_data, "time_var") <- "time"
    attr(std_data, "data_structure") <- "longitudinal"
  } else {
    attr(std_data, "data_structure") <- "cross-sectional"
  }
  
  # Add class
  class(std_data) <- c("standardized_trial", class(std_data))
  
  return(std_data)
}