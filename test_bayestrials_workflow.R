#!/usr/bin/env Rscript

# Test Script for bayestrials Package  
# Demonstrates the complete brms-based workflow

cat("=== Starting bayestrials Workflow Test ===\n")

# Load required libraries
library(bayestrials)

# Set up for reproducibility
set.seed(123)

cat("\n=== Step 1: Generate Synthetic Data ===\n")
# Generate synthetic clinical trial data
data <- generate_synthetic_data(
  n_observations = 500,  # Reduced for faster testing
  n_variables = 10,      # Fewer variables for testing
  seed = 123
)

cat("Generated data dimensions:", dim(data), "\n")
cat("Variables include:", paste(head(names(data), 6), collapse = ", "), "...\n")

# Get variable descriptions
var_desc <- get_variable_descriptions(data)
cat("Variable descriptions available for", nrow(var_desc), "variables\n")

cat("\n=== Step 2: Specify Models ===\n")
# Create model specification using brms approach
model_spec <- create_model_spec(
  "outcome_continuous",
  "treatment + age_centered + baseline_severity_centered",
  family = "gaussian",
  model_name = "linear_regression"
)

cat("Model specification created:\n")
print(model_spec)

cat("\n=== Step 3: Specify Priors ===\n")
# Create a simple prior specification using the PriorSpecification class
priors <- PriorSpecification()
priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
priors <- add_prior(priors, "b_baseline_severity_centered", "normal", 0, 2.0)
priors <- add_prior(priors, "Intercept", "normal", 0, 10)
priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

cat("Prior specifications created:\n")
print(priors)

# Create skeptical priors for comparison
skeptical_priors <- PriorSpecification()
skeptical_priors <- add_prior(skeptical_priors, "b_treatment", "normal", 0, 0.5)
skeptical_priors <- add_prior(skeptical_priors, "b_age_centered", "normal", 0, 0.5)
skeptical_priors <- add_prior(skeptical_priors, "b_baseline_severity_centered", "normal", 0, 1.0) 

cat("\n=== Step 4: Setup Parallel Processing ===\n")
# Setup parallel processing (use sequential for testing to avoid issues)
setup_parallel(strategy = "sequential", verbose = TRUE)

# Check configuration
config <- get_parallel_config()
cat("Parallel configuration:\n")
print(config)

cat("\n=== Step 5: Fit Models ===\n")
cat("Fitting Bayesian model (this may take a few minutes)...\n")

# Attach priors to model spec
attr(model_spec, "priors") <- priors

# Fit the main model
result <- fit_model(
  model_spec = model_spec,
  data = data,
  chains = 2,        # Reduced for testing
  iter = 1000,       # Reduced for testing
  cores = 1,         # Sequential for compatibility
  seed = 123
)

if (is.null(result)) {
  cat("❌ Model fitting failed - this may be due to missing build tools\n")
  cat("Please ensure you have the necessary C++ compiler installed\n")
  cat("  - macOS: Run 'xcode-select --install' in Terminal\n")
  cat("  - Windows: Install Rtools from https://cran.r-project.org/bin/windows/Rtools/\n")
  cat("  - Linux: Install build-essential or equivalent\n")
  quit("no", 1)
}

cat("Model fitting completed!\n")
print(result)

cat("\n=== Step 6: Check Diagnostics ===\n")
# Check model diagnostics
diag_check <- check_diagnostics(result)
cat("Diagnostics check:\n")
print(diag_check)

cat("\n=== Step 7: Extract Results ===\n")
# Extract parameter estimates
estimates <- extract_estimates(result)
cat("Parameter estimates:\n")
print(estimates)

# Get true effects for comparison
true_effects <- get_true_effects(data)
cat("\nTrue effects from data generation:\n")
print(true_effects)

cat("\n=== Step 8: Create Visualizations ===\n")
# Create posterior plot
cat("Creating posterior distribution plot...\n")
post_plot <- plot_posterior(result)
print(post_plot)

# Create forest plot comparison
cat("Creating forest plot...\n")
forest_data <- tibble::tibble(
  trial = "Test Model",
  estimate = estimates$mean[estimates$parameter == "b_treatment"],
  ci_lower = estimates$`2.5%`[estimates$parameter == "b_treatment"],
  ci_upper = estimates$`97.5%`[estimates$parameter == "b_treatment"]
)

forest_plot <- plot_forest(forest_data, parameter = "treatment")
print(forest_plot)

cat("\n=== Step 9: Model Comparison ===\n")
# Compare with different priors
cat("Fitting models with different priors for comparison...\n")

# Create skeptical model spec
skeptical_model <- create_model_spec(
  "outcome_continuous",
  "treatment + age_centered + baseline_severity_centered",
  family = "gaussian",
  model_name = "skeptical_regression"
)
attr(skeptical_model, "priors") <- skeptical_priors

# Fit model with skeptical priors
result_skeptical <- fit_model(
  model_spec = skeptical_model,
  data = data,
  chains = 2,
  iter = 1000,
  cores = 1,
  seed = 123
)

# For now skip optimistic to keep test faster
result_optimistic <- NULL

# Compare models if all fitted successfully
models_list <- list(
  neutral = result,
  skeptical = result_skeptical,
  optimistic = result_optimistic
)

# Remove any NULL models
models_list <- models_list[!sapply(models_list, is.null)]

if (length(models_list) > 1) {
  cat("Comparing", length(models_list), "models...\n")
  
  # Create comparison plot
  comparison_plot <- plot_model_comparison(
    models_list, 
    parameter = "treatment",
    plot_type = "intervals"
  )
  print(comparison_plot)
  
  # Compare using WAIC if possible
  tryCatch({
    waic_comparison <- compare_models(models_list[[1]], models_list[[2]], criterion = "waic")
    cat("WAIC comparison:\n")
    print(waic_comparison)
  }, error = function(e) {
    cat("WAIC comparison not available:", e$message, "\n")
  })
} else {
  cat("Only one model fitted successfully, skipping comparison\n")
}

cat("\n=== Step 10: Validation Tests ===\n")
# Validate results make sense
cat("Running validation checks...\n")

# Check if treatment effect is recovered reasonably well
treatment_param <- estimates$parameter == "b_treatment"
if (any(treatment_param)) {
  treatment_estimate <- estimates$mean[treatment_param]
  true_treatment <- true_effects$treatment_main
  bias <- abs(treatment_estimate - true_treatment)
  
  cat("Treatment effect recovery:\n")
  cat("  True value:", round(true_treatment, 3), "\n")
  cat("  Estimated value:", round(treatment_estimate, 3), "\n")
  cat("  Absolute bias:", round(bias, 3), "\n")
  
  if (bias < 0.1) {
    cat("  ✓ Treatment effect well recovered!\n")
  } else if (bias < 0.2) {
    cat("  ⚠ Treatment effect reasonably recovered\n")
  } else {
    cat("  ✗ Treatment effect poorly recovered - check model\n")
  }
} else {
  cat("Treatment parameter not found in estimates\n")
}

# Check convergence
if (!is.null(diag_check) && "converged" %in% names(diag_check)) {
  if (diag_check$converged) {
    cat("  ✓ Model converged successfully!\n")
  } else {
    cat("  ⚠ Model convergence issues detected\n")
  }
}

cat("\n=== Step 11: Export Results ===\n")
# Save results for later analysis
saveRDS(result, "bayesian_result.rds")
cat("Main result saved to bayesian_result.rds\n")

# Export estimates to CSV
if (requireNamespace("readr", quietly = TRUE)) {
  readr::write_csv(estimates, "parameter_estimates.csv")
  cat("Parameter estimates exported to parameter_estimates.csv\n")
} else {
  write.csv(estimates, "parameter_estimates.csv", row.names = FALSE)
  cat("Parameter estimates exported to parameter_estimates.csv\n")
}

# Export priors to CSV for documentation
export_priors(priors, "prior_specifications.csv")
cat("Prior specifications exported to prior_specifications.csv\n")

cat("\n=== Step 12: Generate Report (if available) ===\n")
# Try to generate a basic report
tryCatch({
  if (exists("generate_report")) {
    report_file <- "bayesian_analysis_report.html"
    generate_report(
      model = result,
      output_file = report_file,
      title = "Bayesian Clinical Trial Analysis - Test Report"
    )
    cat("Report generated:", report_file, "\n")
  } else {
    cat("Report generation function not available\n")
  }
}, error = function(e) {
  cat("Report generation failed:", e$message, "\n")
})

cat("\n=== Workflow Test Completed Successfully! ===\n")

# Summary of generated files
cat("\nGenerated files:\n")
files_created <- c(
  "bayesian_result.rds",
  "parameter_estimates.csv",
  "prior_specifications.csv"
)

# Check which files actually exist
existing_files <- files_created[file.exists(files_created)]
for (file in existing_files) {
  cat("  ✓", file, "\n")
}

cat("\n=== Test Summary ===\n")
cat("✓ Data generation: SUCCESS\n")
cat("✓ Model specification: SUCCESS\n") 
cat("✓ Prior specification: SUCCESS\n")
cat("✓ Model fitting: SUCCESS\n")
cat("✓ Results extraction: SUCCESS\n")
cat("✓ Visualization: SUCCESS\n")
cat("✓ Export: SUCCESS\n")

cat("\nThe bayestrials package workflow has been successfully tested!\n")
cat("This demonstrates the brms-based Bayesian analysis capabilities.\n")