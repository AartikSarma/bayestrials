#!/usr/bin/env Rscript

# Comprehensive Synthetic Data Workflow Test for bayestrials Package
# Tests the complete workflow across all outcome types with parameter recovery validation

cat("===========================================\n")
cat("bayestrials Synthetic Data Workflow Test\n")
cat("===========================================\n")

# Load required libraries
suppressPackageStartupMessages({
  library(bayestrials)
  library(tibble)
  library(dplyr)
})

# Test configuration
CONFIG <- list(
  n_observations = 300,  # Sample size for synthetic data
  n_iter = 1000,         # MCMC iterations (reduced for testing)
  n_chains = 2,          # Number of chains
  n_warmup = 500,        # Warmup iterations
  seed = 42,             # Random seed for reproducibility
  treatment_effect = 0.5 # True treatment effect
)

# Initialize test results tracker
test_results <- list(
  passed = 0,
  failed = 0,
  tests = list()
)

# Helper function to log test results
log_test <- function(test_name, passed, message = "") {
  status <- if (passed) "PASS" else "FAIL"
  icon <- if (passed) "\u2713" else "\u2717"
  cat(sprintf("  [%s] %s: %s\n", status, test_name, message))

  test_results$tests[[test_name]] <<- list(passed = passed, message = message)
  if (passed) {
    test_results$passed <<- test_results$passed + 1
  } else {
    test_results$failed <<- test_results$failed + 1
  }
}

# Helper function for checking parameter recovery
check_parameter_recovery <- function(estimates, true_effects, param_name, true_name, tolerance = 0.3) {
  param_row <- estimates[estimates$parameter == param_name, ]
  if (nrow(param_row) == 0) {
    return(list(recovered = FALSE, message = paste("Parameter", param_name, "not found")))
  }

  estimate <- param_row$mean[1]
  ci_lower <- param_row$`2.5%`[1]
  ci_upper <- param_row$`97.5%`[1]
  true_value <- true_effects[[true_name]]

  # Check if true value is within credible interval
  in_ci <- true_value >= ci_lower && true_value <= ci_upper

  # Check absolute bias
  bias <- abs(estimate - true_value)

  message <- sprintf(
    "True=%.3f, Est=%.3f, 95%%CI=[%.3f, %.3f], Bias=%.3f",
    true_value, estimate, ci_lower, ci_upper, bias
  )

  recovered <- in_ci || bias < tolerance

  return(list(recovered = recovered, message = message, bias = bias, in_ci = in_ci))
}

#==============================================================================
# SECTION 1: Continuous Outcome Workflow
#==============================================================================
cat("\n--- SECTION 1: Continuous Outcome Workflow ---\n")

tryCatch({
  # Step 1.1: Generate continuous outcome data
  cat("\nStep 1.1: Generating synthetic data (continuous outcome)...\n")
  data_continuous <- generate_synthetic_data(
    n_observations = CONFIG$n_observations,
    n_variables = 8,
    outcome_type = "continuous",
    treatment_effect = CONFIG$treatment_effect,
    seed = CONFIG$seed
  )

  log_test(
    "continuous_data_generation",
    nrow(data_continuous) == CONFIG$n_observations && "outcome_continuous" %in% names(data_continuous),
    sprintf("Generated %d observations with %d variables", nrow(data_continuous), ncol(data_continuous))
  )

  # Get true effects for validation
  true_effects_cont <- get_true_effects(data_continuous)
  log_test(
    "continuous_true_effects",
    !is.null(true_effects_cont) && "treatment_main" %in% names(true_effects_cont),
    sprintf("Treatment effect = %.3f", true_effects_cont$treatment_main)
  )

  # Step 1.2: Create model specification
  cat("\nStep 1.2: Creating model specification...\n")
  model_spec_cont <- create_model_spec(
    outcome = "outcome_continuous",
    predictors = "treatment + age_centered + baseline_severity_centered",
    family = "gaussian",
    model_name = "continuous_model"
  )

  log_test(
    "continuous_model_spec",
    inherits(model_spec_cont, "ModelSpecification") && validate_model_spec(model_spec_cont),
    sprintf("Model: %s ~ %s", model_spec_cont$outcome, model_spec_cont$predictors)
  )

  # Step 1.3: Create prior specification
  cat("\nStep 1.3: Creating prior specification...\n")
  priors_cont <- PriorSpecification()
  priors_cont <- add_prior(priors_cont, "b_treatment", "normal", 0, 2.5)
  priors_cont <- add_prior(priors_cont, "b_age_centered", "normal", 0, 1.0)
  priors_cont <- add_prior(priors_cont, "b_baseline_severity_centered", "normal", 0, 2.0)
  priors_cont <- add_prior(priors_cont, "Intercept", "normal", 10, 10)
  priors_cont <- add_prior(priors_cont, "sigma", "half_cauchy", NA, 5, 0)

  brms_priors_cont <- to_brms_prior(priors_cont)
  log_test(
    "continuous_priors",
    inherits(brms_priors_cont, "brmsprior") && nrow(brms_priors_cont) > 0,
    sprintf("Created %d priors", nrow(brms_priors_cont))
  )

  # Step 1.4: Fit model
  cat("\nStep 1.4: Fitting Bayesian model (this may take a few minutes)...\n")
  attr(model_spec_cont, "priors") <- priors_cont

  result_cont <- fit_model(
    model_spec = model_spec_cont,
    data = data_continuous,
    chains = CONFIG$n_chains,
    iter = CONFIG$n_iter,
    warmup = CONFIG$n_warmup,
    cores = 1,
    seed = CONFIG$seed,
    silent = TRUE
  )

  if (is.null(result_cont)) {
    log_test("continuous_model_fit", FALSE, "Model fitting failed - check C++ compiler")
    stop("Model fitting failed")
  }

  log_test(
    "continuous_model_fit",
    inherits(result_cont, "BayesianModel") && !is.null(result_cont$model),
    "Model fitted successfully"
  )

  # Step 1.5: Check diagnostics
  cat("\nStep 1.5: Checking model diagnostics...\n")
  diagnostics_cont <- check_diagnostics(result_cont)

  converged <- diagnostics_cont$converged %||% FALSE
  log_test(
    "continuous_convergence",
    converged,
    sprintf("Rhat_max=%.3f, ESS_bulk_min=%d, Divergences=%d",
            diagnostics_cont$max_rhat %||% NA,
            diagnostics_cont$min_ess_bulk %||% NA,
            diagnostics_cont$n_divergent %||% NA)
  )

  # Step 1.6: Extract estimates and check parameter recovery
  cat("\nStep 1.6: Extracting estimates and checking parameter recovery...\n")
  estimates_cont <- extract_estimates(result_cont)

  # Check treatment effect recovery
  treatment_check <- check_parameter_recovery(
    estimates_cont, true_effects_cont, "b_treatment", "treatment_main"
  )
  log_test(
    "continuous_treatment_recovery",
    treatment_check$recovered,
    treatment_check$message
  )

  # Check age effect recovery
  age_check <- check_parameter_recovery(
    estimates_cont, true_effects_cont, "b_age_centered", "age_effect"
  )
  log_test(
    "continuous_age_recovery",
    age_check$recovered,
    age_check$message
  )

}, error = function(e) {
  cat("Error in continuous workflow:", conditionMessage(e), "\n")
  log_test("continuous_workflow", FALSE, conditionMessage(e))
})

#==============================================================================
# SECTION 2: Binary Outcome Workflow
#==============================================================================
cat("\n--- SECTION 2: Binary Outcome Workflow ---\n")

tryCatch({
  # Step 2.1: Generate binary outcome data
  cat("\nStep 2.1: Generating synthetic data (binary outcome)...\n")
  data_binary <- generate_synthetic_data(
    n_observations = CONFIG$n_observations,
    n_variables = 8,
    outcome_type = "binary",
    treatment_effect = CONFIG$treatment_effect,
    seed = CONFIG$seed + 1
  )

  log_test(
    "binary_data_generation",
    nrow(data_binary) == CONFIG$n_observations && "outcome_binary" %in% names(data_binary),
    sprintf("Generated %d observations, outcome mean = %.3f",
            nrow(data_binary), mean(data_binary$outcome_binary, na.rm = TRUE))
  )

  true_effects_bin <- get_true_effects(data_binary)

  # Step 2.2: Create model specification
  cat("\nStep 2.2: Creating binary model specification...\n")
  model_spec_bin <- create_model_spec(
    outcome = "outcome_binary",
    predictors = "treatment + age_centered + baseline_severity_centered",
    family = "binomial",
    link = "logit",
    model_name = "binary_model"
  )

  log_test(
    "binary_model_spec",
    inherits(model_spec_bin, "ModelSpecification") && model_spec_bin$family == "binomial",
    sprintf("Family: %s, Link: %s", model_spec_bin$family, model_spec_bin$link)
  )

  # Step 2.3: Create priors for binary model
  cat("\nStep 2.3: Creating prior specification for binary model...\n")
  priors_bin <- PriorSpecification()
  priors_bin <- add_prior(priors_bin, "b_treatment", "normal", 0, 2.0)
  priors_bin <- add_prior(priors_bin, "b_age_centered", "normal", 0, 1.0)
  priors_bin <- add_prior(priors_bin, "b_baseline_severity_centered", "normal", 0, 1.0)
  priors_bin <- add_prior(priors_bin, "Intercept", "normal", 0, 5)

  attr(model_spec_bin, "priors") <- priors_bin

  # Step 2.4: Fit binary model
  cat("\nStep 2.4: Fitting binary model...\n")
  result_bin <- fit_model(
    model_spec = model_spec_bin,
    data = data_binary,
    chains = CONFIG$n_chains,
    iter = CONFIG$n_iter,
    warmup = CONFIG$n_warmup,
    cores = 1,
    seed = CONFIG$seed,
    silent = TRUE
  )

  if (!is.null(result_bin)) {
    log_test(
      "binary_model_fit",
      inherits(result_bin, "BayesianModel"),
      "Binary model fitted successfully"
    )

    # Check diagnostics
    diagnostics_bin <- check_diagnostics(result_bin)
    log_test(
      "binary_convergence",
      diagnostics_bin$converged %||% FALSE,
      sprintf("Rhat_max=%.3f", diagnostics_bin$max_rhat %||% NA)
    )

    # Extract and check estimates
    estimates_bin <- extract_estimates(result_bin)
    treatment_check_bin <- check_parameter_recovery(
      estimates_bin, true_effects_bin, "b_treatment", "treatment_main", tolerance = 0.5
    )
    log_test(
      "binary_treatment_recovery",
      treatment_check_bin$recovered,
      treatment_check_bin$message
    )
  } else {
    log_test("binary_model_fit", FALSE, "Model fitting failed")
  }

}, error = function(e) {
  cat("Error in binary workflow:", conditionMessage(e), "\n")
  log_test("binary_workflow", FALSE, conditionMessage(e))
})

#==============================================================================
# SECTION 3: Count Outcome Workflow
#==============================================================================
cat("\n--- SECTION 3: Count Outcome Workflow ---\n")

tryCatch({
  # Step 3.1: Generate count outcome data
  cat("\nStep 3.1: Generating synthetic data (count outcome)...\n")
  data_count <- generate_synthetic_data(
    n_observations = CONFIG$n_observations,
    n_variables = 8,
    outcome_type = "count",
    treatment_effect = CONFIG$treatment_effect,
    seed = CONFIG$seed + 2
  )

  log_test(
    "count_data_generation",
    nrow(data_count) == CONFIG$n_observations && "outcome_count" %in% names(data_count),
    sprintf("Generated %d observations, outcome mean = %.2f",
            nrow(data_count), mean(data_count$outcome_count, na.rm = TRUE))
  )

  # Step 3.2: Create model specification
  cat("\nStep 3.2: Creating count model specification...\n")
  model_spec_count <- create_model_spec(
    outcome = "outcome_count",
    predictors = "treatment + age_centered",
    family = "poisson",
    model_name = "count_model"
  )

  log_test(
    "count_model_spec",
    inherits(model_spec_count, "ModelSpecification") && model_spec_count$family == "poisson",
    sprintf("Family: %s", model_spec_count$family)
  )

  # Step 3.3: Create priors for count model
  cat("\nStep 3.3: Creating prior specification for count model...\n")
  priors_count <- PriorSpecification()
  priors_count <- add_prior(priors_count, "b_treatment", "normal", 0, 1.0)
  priors_count <- add_prior(priors_count, "b_age_centered", "normal", 0, 0.5)
  priors_count <- add_prior(priors_count, "Intercept", "normal", 0, 5)

  attr(model_spec_count, "priors") <- priors_count

  # Step 3.4: Fit count model
  cat("\nStep 3.4: Fitting count model...\n")
  result_count <- fit_model(
    model_spec = model_spec_count,
    data = data_count,
    chains = CONFIG$n_chains,
    iter = CONFIG$n_iter,
    warmup = CONFIG$n_warmup,
    cores = 1,
    seed = CONFIG$seed,
    silent = TRUE
  )

  if (!is.null(result_count)) {
    log_test(
      "count_model_fit",
      inherits(result_count, "BayesianModel"),
      "Count model fitted successfully"
    )

    diagnostics_count <- check_diagnostics(result_count)
    log_test(
      "count_convergence",
      diagnostics_count$converged %||% FALSE,
      sprintf("Rhat_max=%.3f", diagnostics_count$max_rhat %||% NA)
    )
  } else {
    log_test("count_model_fit", FALSE, "Model fitting failed")
  }

}, error = function(e) {
  cat("Error in count workflow:", conditionMessage(e), "\n")
  log_test("count_workflow", FALSE, conditionMessage(e))
})

#==============================================================================
# SECTION 4: Sensitivity Analysis
#==============================================================================
cat("\n--- SECTION 4: Sensitivity Analysis ---\n")

tryCatch({
  # Only run if continuous model was fit successfully
  if (exists("result_cont") && !is.null(result_cont)) {

    cat("\nStep 4.1: Creating skeptical prior specification...\n")
    priors_skeptical <- PriorSpecification()
    priors_skeptical <- add_prior(priors_skeptical, "b_treatment", "normal", 0, 0.5)  # More skeptical
    priors_skeptical <- add_prior(priors_skeptical, "b_age_centered", "normal", 0, 0.5)
    priors_skeptical <- add_prior(priors_skeptical, "b_baseline_severity_centered", "normal", 0, 1.0)
    priors_skeptical <- add_prior(priors_skeptical, "Intercept", "normal", 10, 10)
    priors_skeptical <- add_prior(priors_skeptical, "sigma", "half_cauchy", NA, 5, 0)

    model_spec_skeptical <- create_model_spec(
      outcome = "outcome_continuous",
      predictors = "treatment + age_centered + baseline_severity_centered",
      family = "gaussian",
      model_name = "skeptical_model"
    )
    attr(model_spec_skeptical, "priors") <- priors_skeptical

    cat("\nStep 4.2: Fitting model with skeptical priors...\n")
    result_skeptical <- fit_model(
      model_spec = model_spec_skeptical,
      data = data_continuous,
      chains = CONFIG$n_chains,
      iter = CONFIG$n_iter,
      warmup = CONFIG$n_warmup,
      cores = 1,
      seed = CONFIG$seed,
      silent = TRUE
    )

    if (!is.null(result_skeptical)) {
      log_test(
        "sensitivity_skeptical_fit",
        inherits(result_skeptical, "BayesianModel"),
        "Skeptical model fitted successfully"
      )

      # Compare treatment effect estimates
      est_neutral <- extract_estimates(result_cont)
      est_skeptical <- extract_estimates(result_skeptical)

      treat_neutral <- est_neutral$mean[est_neutral$parameter == "b_treatment"]
      treat_skeptical <- est_skeptical$mean[est_skeptical$parameter == "b_treatment"]

      log_test(
        "sensitivity_comparison",
        TRUE,
        sprintf("Neutral prior: %.3f, Skeptical prior: %.3f, Diff: %.3f",
                treat_neutral, treat_skeptical, abs(treat_neutral - treat_skeptical))
      )
    } else {
      log_test("sensitivity_skeptical_fit", FALSE, "Failed to fit skeptical model")
    }
  } else {
    log_test("sensitivity_analysis", FALSE, "Skipped - no base model available")
  }

}, error = function(e) {
  cat("Error in sensitivity analysis:", conditionMessage(e), "\n")
  log_test("sensitivity_analysis", FALSE, conditionMessage(e))
})

#==============================================================================
# SECTION 5: Visualization Tests
#==============================================================================
cat("\n--- SECTION 5: Visualization Tests ---\n")

tryCatch({
  if (exists("result_cont") && !is.null(result_cont)) {

    cat("\nStep 5.1: Creating posterior distribution plot...\n")
    post_plot <- plot_posterior(result_cont)
    log_test(
      "visualization_posterior",
      inherits(post_plot, "gg") || inherits(post_plot, "ggplot"),
      "Posterior plot created successfully"
    )

    cat("\nStep 5.2: Creating forest plot...\n")
    estimates_for_forest <- extract_estimates(result_cont)
    forest_data <- tibble(
      trial = "Test Model",
      estimate = estimates_for_forest$mean[estimates_for_forest$parameter == "b_treatment"],
      ci_lower = estimates_for_forest$`2.5%`[estimates_for_forest$parameter == "b_treatment"],
      ci_upper = estimates_for_forest$`97.5%`[estimates_for_forest$parameter == "b_treatment"]
    )

    forest_plot <- plot_forest(forest_data, parameter = "treatment")
    log_test(
      "visualization_forest",
      inherits(forest_plot, "gg") || inherits(forest_plot, "ggplot"),
      "Forest plot created successfully"
    )

  } else {
    log_test("visualization_tests", FALSE, "Skipped - no model results available")
  }

}, error = function(e) {
  cat("Error in visualization tests:", conditionMessage(e), "\n")
  log_test("visualization_tests", FALSE, conditionMessage(e))
})

#==============================================================================
# SECTION 6: Export and Reproducibility
#==============================================================================
cat("\n--- SECTION 6: Export and Reproducibility ---\n")

tryCatch({
  # Test prior export/import
  cat("\nStep 6.1: Testing prior export/import cycle...\n")
  temp_file <- tempfile(fileext = ".csv")

  if (exists("priors_cont") && !is.null(priors_cont)) {
    export_priors(priors_cont, temp_file)
    priors_reimported <- import_priors(temp_file)

    log_test(
      "export_import_priors",
      inherits(priors_reimported, "PriorSpecification") &&
        nrow(priors_reimported$priors) == nrow(priors_cont$priors),
      sprintf("Exported and reimported %d priors", nrow(priors_reimported$priors))
    )

    unlink(temp_file)
  }

  # Test results export
  cat("\nStep 6.2: Testing results export...\n")
  if (exists("estimates_cont") && !is.null(estimates_cont)) {
    temp_results <- tempfile(fileext = ".csv")
    write.csv(estimates_cont, temp_results, row.names = FALSE)

    log_test(
      "export_estimates",
      file.exists(temp_results) && file.size(temp_results) > 0,
      sprintf("Exported estimates to CSV (%d bytes)", file.size(temp_results))
    )

    unlink(temp_results)
  }

  # Test session info
  cat("\nStep 6.3: Testing session info capture...\n")
  session_info <- get_session_info()
  log_test(
    "session_info",
    is.list(session_info) && "bayestrials_version" %in% names(session_info),
    sprintf("Version: %s", session_info$bayestrials_version)
  )

}, error = function(e) {
  cat("Error in export tests:", conditionMessage(e), "\n")
  log_test("export_tests", FALSE, conditionMessage(e))
})

#==============================================================================
# TEST SUMMARY
#==============================================================================
cat("\n===========================================\n")
cat("TEST SUMMARY\n")
cat("===========================================\n")

total_tests <- test_results$passed + test_results$failed
pass_rate <- if (total_tests > 0) test_results$passed / total_tests * 100 else 0

cat(sprintf("\nTotal Tests: %d\n", total_tests))
cat(sprintf("Passed: %d (%s)\n", test_results$passed, "\u2713"))
cat(sprintf("Failed: %d (%s)\n", test_results$failed, "\u2717"))
cat(sprintf("Pass Rate: %.1f%%\n", pass_rate))

if (test_results$failed > 0) {
  cat("\nFailed Tests:\n")
  for (test_name in names(test_results$tests)) {
    test_info <- test_results$tests[[test_name]]
    if (!test_info$passed) {
      cat(sprintf("  - %s: %s\n", test_name, test_info$message))
    }
  }
}

cat("\n===========================================\n")

if (test_results$failed == 0) {
  cat("\nAll synthetic data workflow tests PASSED!\n")
  cat("The bayestrials package is working correctly.\n")
} else {
  cat("\nSome tests FAILED. Please review the errors above.\n")
}

cat("\nTest completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Return exit code for CI/CD integration
quit(save = "no", status = if (test_results$failed > 0) 1 else 0)
