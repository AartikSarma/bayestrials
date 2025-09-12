#' PriorSpecification class
#' 
#' @description S3 class for storing and validating prior specifications
#' 
#' @export
PriorSpecification <- function() {
  structure(
    list(
      priors = data.frame(
        parameter = character(),
        distribution = character(),
        location = numeric(),
        scale = numeric(),
        lb = numeric(),
        ub = numeric(),
        stringsAsFactors = FALSE
      )
    ),
    class = "PriorSpecification"
  )
}

#' Print method for PriorSpecification
#' 
#' @param x A PriorSpecification object
#' @param ... Additional arguments passed to print
#' @return x invisibly
#' @export
print.PriorSpecification <- function(x, ...) {
  if (nrow(x$priors) == 0) {
    cat("PriorSpecification: [Empty]\n")
  } else {
    cat("PriorSpecification:\n")
    for (i in 1:nrow(x$priors)) {
      row <- x$priors[i, ]
      cat("  ", row$parameter, " ~ ", row$distribution, "(", sep = "")
      
      params <- c()
      if (!is.na(row$location)) params <- c(params, paste0("location = ", row$location))
      if (!is.na(row$scale)) params <- c(params, paste0("scale = ", row$scale))
      if (!is.na(row$lb)) params <- c(params, paste0("lb = ", row$lb))
      if (!is.na(row$ub)) params <- c(params, paste0("ub = ", row$ub))
      
      cat(paste(params, collapse = ", "), ")\n", sep = "")
    }
  }
  
  invisible(x)
}

#' Add a prior to a PriorSpecification
#' 
#' @param prior_spec A PriorSpecification object
#' @param parameter Name of the parameter
#' @param distribution Name of the distribution
#' @param location Location parameter
#' @param scale Scale parameter
#' @param lb Lower bound
#' @param ub Upper bound
#' @return Updated PriorSpecification object
#' @export
add_prior <- function(prior_spec, parameter, distribution, location = NA, scale = NA, lb = NA, ub = NA) {
  if (!inherits(prior_spec, "PriorSpecification")) {
    stop("prior_spec must be a PriorSpecification object", call. = FALSE)
  }
  
  # Check if parameter already exists
  if (parameter %in% prior_spec$priors$parameter) {
    warning(paste("Prior for parameter", parameter, "already exists and will be overwritten"), call. = FALSE)
    prior_spec$priors <- prior_spec$priors[prior_spec$priors$parameter != parameter, ]
  }
  
  # Validate distribution
  valid_distributions <- c("normal", "student_t", "cauchy", "gamma", "exponential", "beta", "uniform",
                          "half_normal", "half_cauchy", "half_student_t", "lognormal")
  
  if (!distribution %in% valid_distributions) {
    stop(paste("Invalid distribution:", distribution, 
              "Valid distributions are:", paste(valid_distributions, collapse = ", ")), call. = FALSE)
  }
  
  # Create new prior row
  new_prior <- data.frame(
    parameter = parameter,
    distribution = distribution,
    location = as.numeric(location),
    scale = as.numeric(scale),
    lb = as.numeric(lb),
    ub = as.numeric(ub),
    stringsAsFactors = FALSE
  )
  
  # Add prior to the specification
  prior_spec$priors <- rbind(prior_spec$priors, new_prior)
  
  return(prior_spec)
}

#' Remove a prior from a PriorSpecification
#' 
#' @param prior_spec A PriorSpecification object
#' @param parameter Name of the parameter to remove
#' @return Updated PriorSpecification object
#' @export
remove_prior <- function(prior_spec, parameter) {
  if (!inherits(prior_spec, "PriorSpecification")) {
    stop("prior_spec must be a PriorSpecification object", call. = FALSE)
  }
  
  # Check if parameter exists
  if (!parameter %in% prior_spec$priors$parameter) {
    warning(paste("Prior for parameter", parameter, "does not exist"), call. = FALSE)
    return(prior_spec)
  }
  
  # Remove prior
  prior_spec$priors <- prior_spec$priors[prior_spec$priors$parameter != parameter, ]
  
  return(prior_spec)
}

#' Read prior specifications from CSV
#' 
#' @param file_path Path to CSV file with prior specifications
#' @return A PriorSpecification object
#' @export
read_prior_spec <- function(file_path) {
  # Check if file exists
  check_file(file_path, "csv")
  
  # Read specifications
  specs <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  
  # Check required columns
  required_cols <- c("parameter", "distribution")
  missing_cols <- setdiff(required_cols, names(specs))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in prior specifications:", 
              paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  
  # Create PriorSpecification object
  prior_spec <- PriorSpecification()
  
  # Add each prior
  for (i in 1:nrow(specs)) {
    row <- specs[i, ]
    
    # Get parameters with default NA
    location <- if ("location" %in% names(row)) row$location else NA
    scale <- if ("scale" %in% names(row)) row$scale else NA
    lb <- if ("lb" %in% names(row)) row$lb else NA
    ub <- if ("ub" %in% names(row)) row$ub else NA
    
    # Add the prior
    prior_spec <- add_prior(
      prior_spec = prior_spec,
      parameter = row$parameter,
      distribution = row$distribution,
      location = location,
      scale = scale,
      lb = lb,
      ub = ub
    )
  }
  
  return(prior_spec)
}

#' Save prior specifications to CSV
#' 
#' @param prior_spec A PriorSpecification object
#' @param file_path Path to save prior specifications
#' @return Invisibly returns TRUE if successful
#' @export
save_prior_spec <- function(prior_spec, file_path) {
  if (!inherits(prior_spec, "PriorSpecification")) {
    stop("prior_spec must be a PriorSpecification object", call. = FALSE)
  }
  
  # Write to file
  utils::write.csv(prior_spec$priors, file_path, row.names = FALSE)
  
  invisible(TRUE)
}

#' Convert PriorSpecification to brms priors
#' 
#' @param prior_spec A PriorSpecification object
#' @param formula A model formula or ModelSpecification
#' @return A brms prior object
#' @export
to_brms_prior <- function(prior_spec, formula = NULL) {
  if (!inherits(prior_spec, "PriorSpecification")) {
    stop("prior_spec must be a PriorSpecification object", call. = FALSE)
  }
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Convert ModelSpecification to formula if needed
  if (inherits(formula, "ModelSpecification")) {
    formula <- build_formula(formula)
  }
  
  # Create brms priors directly without intermediate list
  combined_priors <- NULL
  
  for (i in 1:nrow(prior_spec$priors)) {
    row <- prior_spec$priors[i, ]
    
    # Extract parameter components
    param_info <- parse_parameter_name(row$parameter)
    
    # Create prior string
    prior_string <- paste0(row$distribution, "(")
    
    # Add parameters based on distribution
    if (row$distribution %in% c("normal", "student_t", "cauchy")) {
      if (!is.na(row$location)) prior_string <- paste0(prior_string, row$location)
      if (!is.na(row$scale)) prior_string <- paste0(prior_string, ", ", row$scale)
    } else if (row$distribution %in% c("gamma", "exponential", "beta")) {
      if (!is.na(row$location)) prior_string <- paste0(prior_string, row$location)
      if (!is.na(row$scale)) prior_string <- paste0(prior_string, ", ", row$scale)
    } else if (row$distribution == "uniform") {
      if (!is.na(row$lb)) prior_string <- paste0(prior_string, row$lb)
      if (!is.na(row$ub)) prior_string <- paste0(prior_string, ", ", row$ub)
    } else if (row$distribution == "half_normal") {
      if (!is.na(row$scale)) prior_string <- paste0(prior_string, "0, ", row$scale)
    } else if (row$distribution == "half_cauchy") {
      # brms uses student_t(1, 0, scale) for half-Cauchy
      prior_string <- "student_t(3, 0, "
      if (!is.na(row$scale)) prior_string <- paste0(prior_string, row$scale)
    } else if (row$distribution == "half_student_t") {
      if (!is.na(row$scale)) prior_string <- paste0(prior_string, "0, ", row$scale)
    } else if (row$distribution == "lognormal") {
      if (!is.na(row$location)) prior_string <- paste0(prior_string, row$location)
      if (!is.na(row$scale)) prior_string <- paste0(prior_string, ", ", row$scale)
    }
    
    prior_string <- paste0(prior_string, ")")
    
    # Skip adding bounds to prior string - half_cauchy and similar already have implicit bounds
    
    # Create brms prior object directly using eval() to avoid NSE issues
    # Handle different parameter classes
    if (!is.na(param_info$class) && param_info$class == "b" && !is.na(param_info$coef)) {
      prior_obj <- eval(call("prior", prior_string, class = "b", coef = param_info$coef), 
                       envir = asNamespace("brms"))
    } else if (!is.na(param_info$class) && param_info$class == "Intercept") {
      prior_obj <- eval(call("prior", prior_string, class = "Intercept"), 
                       envir = asNamespace("brms"))
    } else if (!is.na(param_info$class) && param_info$class %in% c("sigma", "shape", "nu", "phi")) {
      prior_obj <- eval(call("prior", prior_string, class = param_info$class), 
                       envir = asNamespace("brms"))
    } else if (!is.na(param_info$class) && param_info$class == "sd" && !is.na(param_info$group)) {
      prior_obj <- eval(call("prior", prior_string, class = "sd", group = param_info$group), 
                       envir = asNamespace("brms"))
    } else {
      # Default case - use class if available, otherwise default to "b"
      class_to_use <- if (!is.na(param_info$class)) param_info$class else "b"
      prior_obj <- eval(call("prior", prior_string, class = class_to_use), 
                       envir = asNamespace("brms"))
    }
    
    # Combine priors
    if (is.null(combined_priors)) {
      combined_priors <- prior_obj
    } else {
      combined_priors <- combined_priors + prior_obj
    }
  }
  
  return(combined_priors)
}

#' Parse parameter name into components
#' 
#' @param param_name Name of the parameter
#' @return List with class, coef, and group components
#' @keywords internal
parse_parameter_name <- function(param_name) {
  # Initialize components
  class <- NA
  coef <- NA
  group <- NA
  
  # Define patterns for common parameter types
  if (grepl("^b_", param_name)) {
    class <- "b"
    coef_part <- sub("^b_", "", param_name)
    if (coef_part != "") coef <- coef_part
  } else if (grepl("^sd_", param_name)) {
    class <- "sd"
    parts <- strsplit(sub("^sd_", "", param_name), "_")[[1]]
    if (length(parts) > 0) group <- parts[1]
    if (length(parts) > 1) coef <- paste(parts[-1], collapse = "_")
  } else if (grepl("^cor_", param_name)) {
    class <- "cor"
    parts <- strsplit(sub("^cor_", "", param_name), "_")[[1]]
    if (length(parts) > 0) group <- parts[1]
  } else if (param_name == "Intercept") {
    class <- "Intercept"
  } else if (param_name == "sigma") {
    class <- "sigma"
  } else if (param_name == "shape") {
    class <- "shape"
  } else if (param_name == "nu") {
    class <- "nu"
  } else if (param_name == "phi") {
    class <- "phi"
  } else if (param_name == "zi") {
    class <- "zi"
  } else if (param_name == "hu") {
    class <- "hu"
  } else if (param_name == "zeta") {
    class <- "zeta"
  } else if (param_name == "kappa") {
    class <- "kappa"
  } else {
    # Default to setting the whole name as the coefficient
    coef <- param_name
  }
  
  return(list(class = class, coef = coef, group = group))
}

#' Specify priors for a model
#' 
#' @param model_formula Model formula or ModelSpecification
#' @param prior_specs PriorSpecification object or NULL
#' @param prior_type Type of priors to use (default, weakly_informative, skeptical, optimistic, reference)
#' @param auto_scale Logical indicating whether to automatically scale priors
#' @param data Optional data for automatic scaling
#' @return A PriorSpecification object
#' @export
specify_priors <- function(model_formula,
                          prior_specs = NULL,
                          prior_type = c("default", "weakly_informative", "skeptical", "optimistic", "reference"),
                          auto_scale = TRUE,
                          data = NULL) {
  
  prior_type <- match.arg(prior_type)
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Convert ModelSpecification to formula if needed
  if (inherits(model_formula, "ModelSpecification")) {
    formula_obj <- build_formula(model_formula)
    family <- model_formula$family
    link <- model_formula$link
  } else {
    # Use brms to get formula info
    formula_obj <- model_formula
    
    # Extract family and link if possible
    if (inherits(formula_obj, "brmsformula")) {
      family <- formula_obj$family$family
      link <- formula_obj$family$link
    } else {
      family <- "gaussian"
      link <- "identity"
    }
  }
  
  # Get default priors from brms
  brms_priors <- brms::get_prior(formula_obj, data = data, family = paste0(family, "(link = '", link, "')"))
  
  # Create PriorSpecification if not provided
  if (is.null(prior_specs)) {
    prior_specs <- PriorSpecification()
  } else if (!inherits(prior_specs, "PriorSpecification")) {
    stop("prior_specs must be a PriorSpecification object", call. = FALSE)
  }
  
  # Apply prior type
  if (prior_type == "default") {
    # Default uses brms defaults
    return(prior_specs)
  } else {
    # Apply non-default priors based on parameter class
    for (i in 1:nrow(brms_priors)) {
      row <- brms_priors[i, ]
      
      # Skip if there's already a prior for this parameter
      param_name <- construct_parameter_name(row)
      if (param_name %in% prior_specs$priors$parameter) {
        next
      }
      
      # Set prior based on parameter class and prior type
      if (row$class == "b") {
        # Fixed effect prior
        if (prior_type == "weakly_informative") {
          # Weakly informative priors for fixed effects
          if (row$coef == "Intercept" || is.na(row$coef)) {
            prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 10)
          } else {
            prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 2.5)
          }
        } else if (prior_type == "skeptical") {
          # Skeptical priors (centered at 0 with narrow scale)
          if (row$coef == "Intercept" || is.na(row$coef)) {
            prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 5)
          } else {
            prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 1)
          }
        } else if (prior_type == "optimistic") {
          # Optimistic priors (for treatment effects, centered at positive value)
          if (row$coef == "treatment" || grepl("trt|arm", row$coef, ignore.case = TRUE)) {
            if (family == "gaussian") {
              prior_specs <- add_prior(prior_specs, param_name, "normal", 0.5, 1)
            } else if (family %in% c("bernoulli", "binomial")) {
              if (link == "logit") {
                prior_specs <- add_prior(prior_specs, param_name, "normal", 0.8, 1)
              } else {
                prior_specs <- add_prior(prior_specs, param_name, "normal", 0.4, 0.5)
              }
            } else {
              prior_specs <- add_prior(prior_specs, param_name, "normal", 0.5, 1)
            }
          } else if (row$coef == "Intercept" || is.na(row$coef)) {
            prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 10)
          } else {
            prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 2.5)
          }
        } else if (prior_type == "reference") {
          # Reference priors (wider to be less informative)
          if (row$coef == "Intercept" || is.na(row$coef)) {
            prior_specs <- add_prior(prior_specs, param_name, "student_t", 0, 10, 3)
          } else {
            prior_specs <- add_prior(prior_specs, param_name, "student_t", 0, 5, 3)
          }
        }
      } else if (row$class == "sd") {
        # Random effects standard deviation
        if (prior_type == "weakly_informative") {
          prior_specs <- add_prior(prior_specs, param_name, "half_cauchy", NA, 5, 0)
        } else if (prior_type == "skeptical") {
          prior_specs <- add_prior(prior_specs, param_name, "half_normal", NA, 2, 0)
        } else if (prior_type == "optimistic") {
          prior_specs <- add_prior(prior_specs, param_name, "half_cauchy", NA, 10, 0)
        } else if (prior_type == "reference") {
          prior_specs <- add_prior(prior_specs, param_name, "half_student_t", NA, 10, 0, 3)
        }
      } else if (row$class == "sigma") {
        # Residual standard deviation
        if (prior_type == "weakly_informative") {
          prior_specs <- add_prior(prior_specs, param_name, "half_cauchy", NA, 5, 0)
        } else if (prior_type == "skeptical") {
          prior_specs <- add_prior(prior_specs, param_name, "half_normal", NA, 2, 0)
        } else if (prior_type == "optimistic") {
          prior_specs <- add_prior(prior_specs, param_name, "half_cauchy", NA, 10, 0)
        } else if (prior_type == "reference") {
          prior_specs <- add_prior(prior_specs, param_name, "half_student_t", NA, 10, 0, 3)
        }
      } else if (row$class == "Intercept") {
        # Global intercept
        if (prior_type == "weakly_informative") {
          prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 10)
        } else if (prior_type == "skeptical") {
          prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 5)
        } else if (prior_type == "optimistic") {
          prior_specs <- add_prior(prior_specs, param_name, "normal", 0, 10)
        } else if (prior_type == "reference") {
          prior_specs <- add_prior(prior_specs, param_name, "student_t", 0, 10, 3)
        }
      }
      # Add more cases for other parameter classes as needed
    }
  }
  
  # Auto-scale priors if requested and data is provided
  if (auto_scale && !is.null(data)) {
    prior_specs <- auto_scale_priors(prior_specs, formula_obj, data)
  }
  
  return(prior_specs)
}

#' Construct a parameter name from brms prior components
#' 
#' @param row A row from brms get_prior output
#' @return Character string with parameter name
#' @keywords internal
construct_parameter_name <- function(row) {
  if (row$class == "b") {
    if (is.na(row$coef)) {
      return("b")
    } else {
      return(paste0("b_", row$coef))
    }
  } else if (row$class == "sd") {
    if (is.na(row$coef) && !is.na(row$group)) {
      return(paste0("sd_", row$group))
    } else if (!is.na(row$coef) && !is.na(row$group)) {
      return(paste0("sd_", row$group, "_", row$coef))
    } else {
      return("sd")
    }
  } else if (row$class == "cor") {
    if (!is.na(row$group)) {
      return(paste0("cor_", row$group))
    } else {
      return("cor")
    }
  } else {
    return(row$class)
  }
}

#' Auto-scale priors based on data
#' 
#' @param prior_specs PriorSpecification object
#' @param formula Model formula
#' @param data Data frame
#' @return Updated PriorSpecification object
#' @keywords internal
auto_scale_priors <- function(prior_specs, formula, data) {
  # Extract response variable from formula
  response_var <- tryCatch({
    if (inherits(formula, "brmsformula")) {
      all.vars(formula$formula[[2]])
    } else {
      all.vars(stats::as.formula(formula)[[2]])
    }
  }, error = function(e) {
    cli::cli_warn("Could not extract response variable from formula: {e$message}")
    return(NULL)
  })
  
  # If response variable extraction failed, return unchanged priors
  if (is.null(response_var)) {
    return(prior_specs)
  }
  
  # Calculate response scale
  if (response_var %in% names(data)) {
    response_data <- data[[response_var]]
    response_sd <- stats::sd(response_data, na.rm = TRUE)
    response_mean <- mean(response_data, na.rm = TRUE)
    
    # Adjust intercept priors
    for (i in 1:nrow(prior_specs$priors)) {
      row <- prior_specs$priors[i, ]
      
      # Adjust intercept priors based on response scale
      if (row$parameter == "b_Intercept" || row$parameter == "Intercept") {
        # For intercept, scale location to mean and scale parameter to response SD
        if (!is.na(row$location) && row$location == 0) {
          # If location is currently 0, center at response mean
          prior_specs$priors$location[i] <- response_mean
        }
        
        if (!is.na(row$scale)) {
          # Scale parameter should be wider than response SD
          prior_specs$priors$scale[i] <- max(prior_specs$priors$scale[i], 2 * response_sd)
        }
      } 
      # Adjust treatment effect priors
      else if (grepl("^b_treatment|^b_trt|^b_arm", row$parameter, ignore.case = TRUE)) {
        # For treatment effects, scale should be related to response SD
        if (!is.na(row$scale)) {
          # Reasonable prior scale for treatment effects
          prior_specs$priors$scale[i] <- min(prior_specs$priors$scale[i], response_sd)
        }
      }
      # Adjust sigma/sd priors
      else if (row$parameter == "sigma" || grepl("^sd_", row$parameter)) {
        # For variance parameters, scale to response SD
        if (!is.na(row$scale)) {
          prior_specs$priors$scale[i] <- max(prior_specs$priors$scale[i], response_sd)
        }
      }
    }
  }
  
  return(prior_specs)
}

#' Update priors with different levels of skepticism
#' 
#' @param model_spec ModelSpecification object
#' @param skepticism Level of skepticism (high, medium, low)
#' @param treatment_var Name of treatment variable
#' @return Updated ModelSpecification object
#' @export
update_priors <- function(model_spec, 
                         skepticism = c("high", "medium", "low"),
                         treatment_var = "treatment") {
  
  skepticism <- match.arg(skepticism)
  
  if (!inherits(model_spec, "ModelSpecification")) {
    stop("model_spec must be a ModelSpecification object", call. = FALSE)
  }
  
  # Get existing priors if any
  priors <- attr(model_spec, "priors")
  if (is.null(priors)) {
    priors <- PriorSpecification()
  }
  
  # Treatment effect parameter name
  trt_param <- paste0("b_", treatment_var)
  
  # Update or add treatment effect prior based on skepticism
  if (skepticism == "high") {
    # Highly skeptical prior: centered at zero with narrow scale
    priors <- add_prior(priors, trt_param, "normal", 0, 0.5)
  } else if (skepticism == "medium") {
    # Moderately skeptical prior: centered at zero with moderate scale
    priors <- add_prior(priors, trt_param, "normal", 0, 1.0)
  } else if (skepticism == "low") {
    # Weakly skeptical prior: centered at zero with wider scale
    priors <- add_prior(priors, trt_param, "normal", 0, 2.5)
  }
  
  # Attach updated priors to model spec
  attr(model_spec, "priors") <- priors
  
  return(model_spec)
}

#' Visualize prior distributions
#' 
#' @param prior_specs PriorSpecification object
#' @param parameters Vector of parameter names to visualize (NULL for all)
#' @param n_samples Number of samples to draw from each prior
#' @return A ggplot object
#' @export
plot_priors <- function(prior_specs, parameters = NULL, n_samples = 1000) {
  if (!inherits(prior_specs, "PriorSpecification")) {
    stop("prior_specs must be a PriorSpecification object", call. = FALSE)
  }
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # If no parameters specified, use all
  if (is.null(parameters)) {
    parameters <- prior_specs$priors$parameter
  } else {
    # Check if specified parameters exist
    missing_params <- setdiff(parameters, prior_specs$priors$parameter)
    if (length(missing_params) > 0) {
      warning(paste("Some parameters not found in prior_specs:", 
                   paste(missing_params, collapse = ", ")), call. = FALSE)
      parameters <- intersect(parameters, prior_specs$priors$parameter)
    }
  }
  
  # If no valid parameters, return NULL
  if (length(parameters) == 0) {
    warning("No valid parameters to plot", call. = FALSE)
    return(NULL)
  }
  
  # Draw samples from each prior
  samples_list <- list()
  
  for (param in parameters) {
    # Get prior specification
    prior_row <- prior_specs$priors[prior_specs$priors$parameter == param, ]
    
    # Draw samples based on distribution
    samples <- draw_prior_samples(prior_row, n_samples)
    
    # Store samples
    samples_list[[param]] <- data.frame(
      parameter = param,
      value = samples,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all samples
  all_samples <- do.call(rbind, samples_list)
  
  # Create plot
  plot <- ggplot2::ggplot(all_samples, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density(fill = "lightblue", alpha = 0.5) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Prior Distributions",
      x = "Parameter Value",
      y = "Density"
    )
  
  return(plot)
}

#' Draw samples from a prior distribution
#' 
#' @param prior_row Row from priors data frame
#' @param n Number of samples to draw
#' @return Vector of samples
#' @keywords internal
draw_prior_samples <- function(prior_row, n = 1000) {
  dist <- prior_row$distribution
  
  # Extract parameters
  location <- prior_row$location
  scale <- prior_row$scale
  lb <- prior_row$lb
  ub <- prior_row$ub
  
  # Handle NA parameters
  if (is.na(location)) location <- 0
  if (is.na(scale)) scale <- 1
  if (is.na(lb)) lb <- -Inf
  if (is.na(ub)) ub <- Inf
  
  # Draw samples based on distribution
  if (dist == "normal") {
    samples <- stats::rnorm(n, mean = location, sd = scale)
  } else if (dist == "student_t") {
    # Assuming df=3 if not specified
    df <- if (is.na(prior_row$df)) 3 else prior_row$df
    samples <- location + scale * stats::rt(n, df = df)
  } else if (dist == "cauchy") {
    samples <- location + scale * stats::rcauchy(n)
  } else if (dist == "gamma") {
    # For gamma, interpret location as shape and scale as scale
    samples <- stats::rgamma(n, shape = location, scale = scale)
  } else if (dist == "exponential") {
    # For exponential, interpret scale as rate
    samples <- stats::rexp(n, rate = 1/scale)
  } else if (dist == "beta") {
    # For beta, interpret location as alpha and scale as beta
    samples <- stats::rbeta(n, shape1 = location, shape2 = scale)
  } else if (dist == "uniform") {
    samples <- stats::runif(n, min = lb, max = ub)
  } else if (dist == "half_normal") {
    samples <- abs(stats::rnorm(n, mean = 0, sd = scale))
  } else if (dist == "half_cauchy") {
    samples <- abs(stats::rcauchy(n, location = 0, scale = scale))
  } else if (dist == "half_student_t") {
    # Assuming df=3 if not specified
    df <- if (is.na(prior_row$df)) 3 else prior_row$df
    samples <- abs(stats::rt(n, df = df) * scale)
  } else if (dist == "lognormal") {
    samples <- stats::rlnorm(n, meanlog = location, sdlog = scale)
  } else {
    # Default to normal
    warning(paste("Unrecognized distribution:", dist, "- using normal"), call. = FALSE)
    samples <- stats::rnorm(n, mean = location, sd = scale)
  }
  
  # Apply bounds
  samples <- pmax(samples, lb)
  samples <- pmin(samples, ub)
  
  return(samples)
}