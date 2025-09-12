#' ModelSpecification class
#' 
#' @description S3 class for storing and validating model specifications
#' 
#' @export
ModelSpecification <- function(
  model_name = NULL,
  outcome = NULL,
  predictors = NULL,
  random_effects = NULL,
  family = "gaussian",
  link = "identity",
  additional_args = list()
) {
  structure(
    list(
      model_name = model_name,
      outcome = outcome,
      predictors = predictors,
      random_effects = random_effects,
      family = family,
      link = link,
      additional_args = additional_args
    ),
    class = "ModelSpecification"
  )
}

#' Print method for ModelSpecification
#' 
#' @param x A ModelSpecification object
#' @param ... Additional arguments passed to print
#' @return x invisibly
#' @export
print.ModelSpecification <- function(x, ...) {
  cat("Model Specification:", ifelse(is.null(x$model_name), "[Unnamed]", x$model_name), "\n")
  cat("  Formula: ", build_formula_string(x), "\n", sep = "")
  cat("  Family: ", x$family, "\n", sep = "")
  cat("  Link: ", x$link, "\n", sep = "")
  
  if (length(x$additional_args) > 0) {
    cat("  Additional arguments:\n")
    for (arg_name in names(x$additional_args)) {
      cat("    ", arg_name, ": ", format_value(x$additional_args[[arg_name]]), "\n", sep = "")
    }
  }
  
  invisible(x)
}

#' Build a formula string from ModelSpecification components
#' 
#' @param spec A ModelSpecification object
#' @return Character string representing the formula
#' @keywords internal
build_formula_string <- function(spec) {
  formula_str <- paste0(spec$outcome, " ~ ", spec$predictors)
  
  if (!is.null(spec$random_effects) && nchar(spec$random_effects) > 0) {
    formula_str <- paste0(formula_str, " + ", spec$random_effects)
  }
  
  return(formula_str)
}

#' Build a brms formula from ModelSpecification
#' 
#' @param spec A ModelSpecification object
#' @return A brms formula object
#' @export
build_formula <- function(spec) {
  if (!inherits(spec, "ModelSpecification")) {
    stop("Input must be a ModelSpecification object.", call. = FALSE)
  }
  
  formula_str <- build_formula_string(spec)
  
  # Check if brms is available
  if (requireNamespace("brms", quietly = TRUE)) {
    formula <- brms::brmsformula(formula_str)
    return(formula)
  } else {
    warning("brms package not available, returning formula as string", call. = FALSE)
    return(formula_str)
  }
}

#' Validate a ModelSpecification object
#' 
#' @param spec A ModelSpecification object
#' @return Logical indicating if the specification is valid
#' @export
validate_model_spec <- function(spec) {
  if (!inherits(spec, "ModelSpecification")) {
    stop("Input must be a ModelSpecification object.", call. = FALSE)
  }
  
  # Check required fields
  if (is.null(spec$outcome) || is.null(spec$predictors)) {
    stop("Model specification must include outcome and predictors.", call. = FALSE)
  }
  
  # Check family and link
  valid_families <- c("gaussian", "bernoulli", "binomial", "poisson", "negbinomial", 
                     "student", "beta", "gamma", "exponential", "weibull", "lognormal")
  
  if (!spec$family %in% valid_families) {
    stop(paste0("Invalid family: ", spec$family, ". Valid families are: ", 
               paste(valid_families, collapse = ", ")), call. = FALSE)
  }
  
  # Different families have different valid links
  valid_links <- switch(spec$family,
                        gaussian = c("identity", "log", "inverse"),
                        bernoulli = c("logit", "probit", "cauchit", "cloglog"),
                        binomial = c("logit", "probit", "cauchit", "cloglog"),
                        poisson = c("log", "identity"),
                        negbinomial = c("log", "identity"),
                        student = c("identity", "log", "inverse"),
                        gamma = c("inverse", "identity", "log"),
                        exponential = c("identity", "log"),
                        weibull = c("identity", "log"),
                        lognormal = c("identity", "log"),
                        beta = c("logit", "probit", "cloglog", "identity"),
                        c("identity", "log", "inverse", "logit", "probit")) # default set
  
  if (!spec$link %in% valid_links) {
    stop(paste0("Invalid link function '", spec$link, "' for family '", spec$family, 
               "'. Valid link functions are: ", paste(valid_links, collapse = ", ")), 
         call. = FALSE)
  }
  
  return(TRUE)
}

#' Read model specification from CSV
#' 
#' @param model_file Path to model specification CSV file
#' @param prior_file Path to prior specification CSV file (optional)
#' @param covariate_file Path to covariate specification CSV file (optional)
#' @param validate Logical indicating whether to validate specifications
#' @return List of ModelSpecification objects
#' @export
read_model_spec <- function(model_file, prior_file = NULL, covariate_file = NULL, validate = TRUE) {
  # Check if files exist
  check_file(model_file, "csv")
  
  if (!is.null(prior_file)) {
    check_file(prior_file, "csv")
  }
  
  if (!is.null(covariate_file)) {
    check_file(covariate_file, "csv")
  }
  
  # Read model specs
  model_specs <- utils::read.csv(model_file, stringsAsFactors = FALSE)
  
  # Read prior specs if provided
  prior_specs <- NULL
  if (!is.null(prior_file)) {
    prior_specs <- utils::read.csv(prior_file, stringsAsFactors = FALSE)
  }
  
  # Read covariate specs if provided
  covariate_specs <- NULL
  if (!is.null(covariate_file)) {
    covariate_specs <- utils::read.csv(covariate_file, stringsAsFactors = FALSE)
  }
  
  # Create model specification objects
  model_list <- list()
  
  for (i in 1:nrow(model_specs)) {
    row <- model_specs[i, ]
    
    # Create basic model spec
    spec <- ModelSpecification(
      model_name = row$model_name,
      outcome = row$outcome,
      predictors = row$predictors,
      random_effects = if ("random_effects" %in% names(row)) row$random_effects else NULL,
      family = row$family,
      link = row$link,
      additional_args = list()
    )
    
    # Set any additional columns as additional args
    extra_cols <- setdiff(names(row), c("model_name", "outcome", "predictors", "random_effects", "family", "link"))
    for (col in extra_cols) {
      if (!is.na(row[[col]]) && nchar(row[[col]]) > 0) {
        spec$additional_args[[col]] <- row[[col]]
      }
    }
    
    # Validate if requested
    if (validate) {
      validate_model_spec(spec)
    }
    
    model_list[[row$model_name]] <- spec
  }
  
  # Add the covariate and prior specs as attributes
  attr(model_list, "prior_specs") <- prior_specs
  attr(model_list, "covariate_specs") <- covariate_specs
  
  return(model_list)
}

#' Create a comprehensive model specification for Bayesian clinical trial analysis
#'
#' This function creates a `ModelSpecification` object that defines the structure
#' of your Bayesian model. The specification includes the outcome variable,
#' predictor variables, statistical family, and optional components like random
#' effects and priors. This is typically the first step in any bayestrials analysis.
#'
#' @param outcome Character string specifying the name of the outcome (dependent) 
#'   variable in your dataset. This should match exactly with a column name in 
#'   your data frame.
#'   
#' @param predictors Character string specifying the predictor (independent) 
#'   variables using standard R formula syntax. Examples:
#'   \itemize{
#'     \item \code{"treatment"} - Single predictor
#'     \item \code{"treatment + age + sex"} - Multiple predictors  
#'     \item \code{"treatment * age + sex"} - Interaction between treatment and age
#'     \item \code{"treatment + I(age^2)"} - Non-linear age effect
#'   }
#'   Variable names should match columns in your dataset.
#'   
#' @param random_effects Character string specifying random effects using lme4 
#'   syntax (optional). Common patterns:
#'   \itemize{
#'     \item \code{"(1|site)"} - Random intercepts by site
#'     \item \code{"(treatment|site)"} - Random slopes for treatment by site
#'     \item \code{"(1|site) + (1|patient)"} - Nested random effects
#'   }
#'   Leave as \code{NULL} for models without random effects.
#'   
#' @param family Character string specifying the statistical distribution family.
#'   Options include:
#'   \itemize{
#'     \item \code{"gaussian"} - Normal distribution (default, for continuous outcomes)
#'     \item \code{"binomial"} - Binomial distribution (for binary outcomes like success/failure)
#'     \item \code{"poisson"} - Poisson distribution (for count data)
#'     \item \code{"negbinomial"} - Negative binomial (for over-dispersed count data)
#'     \item \code{"gamma"} - Gamma distribution (for positive continuous data)
#'     \item \code{"student_t"} - Student-t distribution (robust to outliers)
#'     \item \code{"beta"} - Beta distribution (for proportions between 0 and 1)
#'   }
#'   
#' @param link Character string specifying the link function that connects the 
#'   linear predictor to the outcome. Common combinations:
#'   \itemize{
#'     \item \code{"identity"} - Direct linear relationship (default for gaussian)
#'     \item \code{"logit"} - Logistic link (default for binomial)  
#'     \item \code{"log"} - Log link (default for poisson, also used for gamma)
#'     \item \code{"probit"} - Probit link (alternative to logit for binary data)
#'     \item \code{"inverse"} - Inverse link (sometimes used with gamma)
#'   }
#'   
#' @param priors A `PriorSpecification` object containing prior distributions 
#'   for model parameters (optional). If provided, these priors will be attached
#'   to the model specification. Can also be added later using \code{attr()}.
#'   
#' @param model_name Character string providing a descriptive name for your model
#'   (optional). This name will appear in output summaries and reports. Examples:
#'   \code{"primary_efficacy"}, \code{"safety_analysis"}, \code{"sensitivity_1"}.
#'
#' @return A `ModelSpecification` object that can be used with [fit_model()].
#'   The object contains all the information needed to fit the Bayesian model
#'   and includes validation to ensure the specification is coherent.
#'
#' @details
#' ## Model Specification Process
#' 
#' The function performs several validation steps:
#' 1. Checks that the outcome variable name is provided
#' 2. Validates that predictors follow proper formula syntax
#' 3. Ensures family and link function combinations are valid
#' 4. Verifies random effects syntax (if provided)
#' 
#' ## Family-Link Combinations
#' 
#' Valid family-link combinations include:
#' - **Gaussian**: identity, log, inverse
#' - **Binomial**: logit, probit, cloglog, identity  
#' - **Poisson**: log, identity
#' - **Gamma**: inverse, identity, log
#' - **Student-t**: identity, log, inverse
#' 
#' ## Clinical Trial Context
#' 
#' In clinical trials, common patterns include:
#' - **Continuous outcomes** (e.g., blood pressure, weight): Use `family = "gaussian"`
#' - **Binary outcomes** (e.g., response/non-response): Use `family = "binomial"`  
#' - **Time-to-event** (with discrete time): Use `family = "binomial"` with appropriate predictors
#' - **Count outcomes** (e.g., number of adverse events): Use `family = "poisson"`
#' - **Biomarker data** (positive continuous): Consider `family = "gamma"` or `family = "lognormal"`
#'
#' @examples
#' # Basic continuous outcome model
#' continuous_spec <- create_model_spec(
#'   outcome = "systolic_bp",
#'   predictors = "treatment + age + baseline_bp",
#'   family = "gaussian",
#'   model_name = "primary_efficacy"
#' )
#' 
#' # Binary outcome with interaction
#' binary_spec <- create_model_spec(
#'   outcome = "response", 
#'   predictors = "treatment * biomarker_high + age + sex",
#'   family = "binomial",
#'   link = "logit",
#'   model_name = "response_analysis"
#' )
#' 
#' # Multi-site trial with random effects
#' hierarchical_spec <- create_model_spec(
#'   outcome = "pain_score",
#'   predictors = "treatment + baseline_pain + age",
#'   random_effects = "(1|site)",
#'   family = "gaussian", 
#'   model_name = "multi_site_analysis"
#' )
#' 
#' # Count data (adverse events)
#' count_spec <- create_model_spec(
#'   outcome = "ae_count",
#'   predictors = "treatment + age + comorbidity_score",
#'   family = "poisson",
#'   link = "log",
#'   model_name = "safety_analysis"
#' )
#' 
#' # With prior specifications
#' library(bayestrials)
#' priors <- PriorSpecification()
#' priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
#' priors <- add_prior(priors, "Intercept", "normal", 120, 20)  # For BP data
#' 
#' spec_with_priors <- create_model_spec(
#'   outcome = "systolic_bp",
#'   predictors = "treatment + age",
#'   family = "gaussian",
#'   priors = priors,
#'   model_name = "bp_analysis_with_priors"
#' )
#'
#' @seealso 
#' - [fit_model()] to fit the model using the specification
#' - [PriorSpecification()] to create prior specifications
#' - [validate_model_spec()] to check specification validity
#' - [update_model_spec()] to modify existing specifications
#' - [read_model_spec()] to read specifications from CSV files
#'
#' @export
create_model_spec <- function(outcome, 
                             predictors, 
                             random_effects = NULL,
                             family = "gaussian", 
                             link = "identity",
                             priors = NULL,
                             model_name = NULL) {
  
  spec <- ModelSpecification(
    model_name = model_name,
    outcome = outcome,
    predictors = predictors,
    random_effects = random_effects,
    family = family,
    link = link
  )
  
  # Validate the spec
  validate_model_spec(spec)
  
  # Attach priors if provided
  if (!is.null(priors)) {
    attr(spec, "priors") <- priors
  }
  
  return(spec)
}

#' Update existing model specification
#' 
#' @param model_spec A ModelSpecification object
#' @param ... Additional parameters to update
#' @param priors New prior specifications (optional)
#' @return Updated ModelSpecification object
#' @export
update_model_spec <- function(model_spec, ..., priors = NULL) {
  if (!inherits(model_spec, "ModelSpecification")) {
    stop("Input must be a ModelSpecification object.", call. = FALSE)
  }
  
  # Process additional arguments
  args <- list(...)
  
  # Update the model spec fields
  for (field in names(args)) {
    if (field %in% names(model_spec)) {
      model_spec[[field]] <- args[[field]]
    } else {
      model_spec$additional_args[[field]] <- args[[field]]
    }
  }
  
  # Update priors if provided
  if (!is.null(priors)) {
    attr(model_spec, "priors") <- priors
  }
  
  # Validate the updated spec
  validate_model_spec(model_spec)
  
  return(model_spec)
}

#' Save model specifications to CSV
#' 
#' @param model_specs List of ModelSpecification objects
#' @param model_file Path to save model specifications
#' @param prior_file Path to save prior specifications (optional)
#' @param covariate_file Path to save covariate specifications (optional)
#' @return Invisibly returns TRUE if successful
#' @export
save_model_specs <- function(model_specs, model_file, prior_file = NULL, covariate_file = NULL) {
  if (!is.list(model_specs) || !all(sapply(model_specs, inherits, "ModelSpecification"))) {
    stop("model_specs must be a list of ModelSpecification objects.", call. = FALSE)
  }
  
  # Convert model specs to data frame
  model_df <- data.frame(
    model_name = sapply(model_specs, function(x) x$model_name),
    outcome = sapply(model_specs, function(x) x$outcome),
    predictors = sapply(model_specs, function(x) x$predictors),
    random_effects = sapply(model_specs, function(x) x$random_effects),
    family = sapply(model_specs, function(x) x$family),
    link = sapply(model_specs, function(x) x$link),
    stringsAsFactors = FALSE
  )
  
  # Add any additional args that are present in all models
  all_add_args <- unique(unlist(lapply(model_specs, function(x) names(x$additional_args))))
  
  for (arg in all_add_args) {
    model_df[[arg]] <- sapply(model_specs, function(x) {
      if (arg %in% names(x$additional_args)) {
        return(x$additional_args[[arg]])
      } else {
        return(NA)
      }
    })
  }
  
  # Write model specifications to file
  utils::write.csv(model_df, model_file, row.names = FALSE)
  
  # Write prior specs if provided
  if (!is.null(prior_file) && !is.null(attr(model_specs, "prior_specs"))) {
    utils::write.csv(attr(model_specs, "prior_specs"), prior_file, row.names = FALSE)
  }
  
  # Write covariate specs if provided
  if (!is.null(covariate_file) && !is.null(attr(model_specs, "covariate_specs"))) {
    utils::write.csv(attr(model_specs, "covariate_specs"), covariate_file, row.names = FALSE)
  }
  
  invisible(TRUE)
}