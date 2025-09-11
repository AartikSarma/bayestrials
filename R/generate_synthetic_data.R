#' Generate Synthetic Clinical Trial Data
#'
#' Creates synthetic clinical trial data with known relationships for testing
#' and demonstration purposes.
#'
#' @param n_observations Number of observations to generate
#' @param n_variables Number of variables to include
#' @param treatment_effect Size of treatment effect
#' @param error_sd Standard deviation of error term
#' @param seed Random seed for reproducibility
#' @param outcome_type Type of outcome ("continuous", "binary", "count")
#' @param missing_rate Proportion of missing data
#'
#' @return A tibble with synthetic clinical trial data
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate continuous outcome data
#' data <- generate_synthetic_data(n_observations = 1000, seed = 123)
#' 
#' # Generate binary outcome data
#' data_binary <- generate_synthetic_data(
#'   n_observations = 500, 
#'   outcome_type = "binary", 
#'   seed = 456
#' )
#' }
generate_synthetic_data <- function(n_observations = 1000,
                                   n_variables = 15,
                                   treatment_effect = 0.5,
                                   error_sd = 1.0,
                                   seed = NULL,
                                   outcome_type = c("continuous", "binary", "count"),
                                   missing_rate = 0.05) {
  
  outcome_type <- match.arg(outcome_type)
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate baseline covariates
  age <- stats::rnorm(n_observations, mean = 65, sd = 10)
  age_centered <- scale(age, center = TRUE, scale = FALSE)[, 1]
  
  baseline_severity <- stats::rnorm(n_observations, mean = 50, sd = 15)
  baseline_severity_centered <- scale(baseline_severity, center = TRUE, scale = FALSE)[, 1]
  
  # Generate treatment assignment (randomized)
  treatment <- sample(c(0, 1), n_observations, replace = TRUE)
  
  # Additional covariates
  sex <- sample(c("Male", "Female"), n_observations, replace = TRUE, prob = c(0.6, 0.4))
  sex_male <- as.numeric(sex == "Male")
  
  # Biomarker with some correlation to baseline severity
  biomarker <- baseline_severity * 0.3 + stats::rnorm(n_observations, sd = 5)
  biomarker_centered <- scale(biomarker, center = TRUE, scale = FALSE)[, 1]
  
  # Generate additional noise variables
  additional_vars <- list()
  for (i in 1:(n_variables - 7)) {
    var_name <- paste0("covariate_", i)
    additional_vars[[var_name]] <- stats::rnorm(n_observations)
  }
  
  # Define true effects
  true_effects <- list(
    intercept = ifelse(outcome_type == "continuous", 10, -1),
    treatment_main = treatment_effect,
    age_effect = ifelse(outcome_type == "continuous", 0.02, 0.01),
    baseline_severity_effect = ifelse(outcome_type == "continuous", 0.8, 0.02),
    sex_effect = ifelse(outcome_type == "continuous", 0.3, 0.15),
    biomarker_effect = ifelse(outcome_type == "continuous", 0.1, 0.005),
    treatment_age_interaction = ifelse(outcome_type == "continuous", 0.01, 0.005)
  )
  
  # Generate linear predictor
  linear_predictor <- true_effects$intercept +
    true_effects$treatment_main * treatment +
    true_effects$age_effect * age_centered +
    true_effects$baseline_severity_effect * baseline_severity_centered +
    true_effects$sex_effect * sex_male +
    true_effects$biomarker_effect * biomarker_centered +
    true_effects$treatment_age_interaction * treatment * age_centered
  
  # Generate outcome based on type
  if (outcome_type == "continuous") {
    outcome_continuous <- linear_predictor + stats::rnorm(n_observations, sd = error_sd)
    outcome <- outcome_continuous
    
  } else if (outcome_type == "binary") {
    # Convert to probabilities using logistic link
    prob <- stats::plogis(linear_predictor)
    outcome_binary <- stats::rbinom(n_observations, 1, prob)
    outcome <- outcome_binary
    
  } else if (outcome_type == "count") {
    # Convert to rates using log link (ensure positive)
    rate <- exp(linear_predictor - max(linear_predictor) + 2)  # Shift to avoid overflow
    outcome_count <- stats::rpois(n_observations, lambda = rate)
    outcome <- outcome_count
  }
  
  # Create data frame
  data <- tibble::tibble(
    subject_id = 1:n_observations,
    treatment = treatment,
    age = age,
    age_centered = age_centered,
    baseline_severity = baseline_severity,
    baseline_severity_centered = baseline_severity_centered,
    sex = sex,
    sex_male = sex_male,
    biomarker = biomarker,
    biomarker_centered = biomarker_centered
  )
  
  # Add outcome variable(s)
  if (outcome_type == "continuous") {
    data$outcome_continuous <- outcome
  } else if (outcome_type == "binary") {
    data$outcome_binary <- outcome
  } else if (outcome_type == "count") {
    data$outcome_count <- outcome
  }
  
  # Add additional variables
  for (var_name in names(additional_vars)) {
    data[[var_name]] <- additional_vars[[var_name]]
  }
  
  # Add some missing data if requested
  if (missing_rate > 0) {
    missing_vars <- c("age", "baseline_severity", "biomarker")
    
    for (var in missing_vars) {
      if (var %in% names(data)) {
        n_missing <- floor(n_observations * missing_rate)
        missing_indices <- sample(n_observations, n_missing)
        data[[var]][missing_indices] <- NA
        
        # Update centered versions
        centered_var <- paste0(var, "_centered")
        if (centered_var %in% names(data)) {
          data[[centered_var]] <- scale(data[[var]], center = TRUE, scale = FALSE)[, 1]
        }
      }
    }
  }
  
  # Add attributes for true effects
  attr(data, "true_effects") <- true_effects
  attr(data, "outcome_type") <- outcome_type
  attr(data, "seed") <- seed
  attr(data, "n_observations") <- n_observations
  
  return(data)
}

#' Extract True Effects from Synthetic Data
#'
#' Extracts the true parameter values used to generate synthetic data.
#'
#' @param data Synthetic data generated by generate_synthetic_data()
#'
#' @return A list with true effect sizes
#' @export
get_true_effects <- function(data) {
  true_effects <- attr(data, "true_effects")
  
  if (is.null(true_effects)) {
    cli::cli_warn("No true effects found in data attributes")
    return(NULL)
  }
  
  return(true_effects)
}

#' Get Variable Descriptions for Synthetic Data
#'
#' Provides descriptions of variables in synthetic clinical trial data.
#'
#' @param data Synthetic data generated by generate_synthetic_data()
#'
#' @return A tibble with variable descriptions
#' @export
get_variable_descriptions <- function(data) {
  
  # Base descriptions
  descriptions <- tibble::tibble(
    variable = c("subject_id", "treatment", "age", "age_centered", 
                 "baseline_severity", "baseline_severity_centered",
                 "sex", "sex_male", "biomarker", "biomarker_centered"),
    description = c(
      "Subject identifier",
      "Treatment assignment (0 = control, 1 = treatment)",
      "Age in years",
      "Age centered at sample mean",
      "Baseline severity score",
      "Baseline severity centered at sample mean",
      "Sex (Male/Female)",
      "Sex indicator (1 = Male, 0 = Female)",
      "Biomarker value",
      "Biomarker centered at sample mean"
    ),
    type = c("identifier", "binary", "continuous", "continuous",
             "continuous", "continuous", "categorical", "binary",
             "continuous", "continuous"),
    role = c("id", "treatment", "covariate", "covariate",
             "covariate", "covariate", "covariate", "covariate", 
             "covariate", "covariate")
  )
  
  # Add outcome descriptions based on what's in the data
  if ("outcome_continuous" %in% names(data)) {
    descriptions <- dplyr::bind_rows(
      descriptions,
      tibble::tibble(
        variable = "outcome_continuous",
        description = "Continuous outcome measure",
        type = "continuous",
        role = "outcome"
      )
    )
  }
  
  if ("outcome_binary" %in% names(data)) {
    descriptions <- dplyr::bind_rows(
      descriptions,
      tibble::tibble(
        variable = "outcome_binary", 
        description = "Binary outcome (0/1)",
        type = "binary",
        role = "outcome"
      )
    )
  }
  
  if ("outcome_count" %in% names(data)) {
    descriptions <- dplyr::bind_rows(
      descriptions,
      tibble::tibble(
        variable = "outcome_count",
        description = "Count outcome",
        type = "count", 
        role = "outcome"
      )
    )
  }
  
  # Add additional covariate descriptions
  additional_vars <- names(data)[grepl("^covariate_", names(data))]
  if (length(additional_vars) > 0) {
    additional_descriptions <- tibble::tibble(
      variable = additional_vars,
      description = paste("Additional covariate", gsub("covariate_", "", additional_vars)),
      type = "continuous",
      role = "covariate"
    )
    
    descriptions <- dplyr::bind_rows(descriptions, additional_descriptions)
  }
  
  # Filter to only variables actually in the data
  descriptions <- descriptions[descriptions$variable %in% names(data), ]
  
  return(descriptions)
}