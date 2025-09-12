#!/usr/bin/env Rscript

# Verbose version of workflow test to debug the caching issue

devtools::load_all()

cat("=== Step 1: Generate Data ===\n")
set.seed(123)
data <- generate_synthetic_data(
  n_observations = 50,  # Smaller for faster testing
  n_variables = 5,
  outcome_type = "continuous",
  seed = 123
)

cat("=== Step 2: Create Model Spec ===\n")
model_spec <- create_model_spec(
  outcome = "outcome_continuous",
  predictors = "treatment + age_centered + baseline_severity_centered",
  family = "gaussian",
  model_name = "workflow_test"
)

cat("=== Step 3: Create Priors ===\n")
priors <- PriorSpecification()
priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
priors <- add_prior(priors, "b_baseline_severity_centered", "normal", 0, 2.0)
priors <- add_prior(priors, "Intercept", "normal", 0, 10)
priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

# Attach priors
attr(model_spec, "priors") <- priors

cat("=== Step 4: Test Prior Conversion ===\n")
converted_priors <- to_brms_prior(priors, model_spec)
cat("Converted priors:\n")
print(converted_priors)

cat("=== Step 5: Test fit_model ===\n")
tryCatch({
  result <- fit_model(
    model_spec = model_spec,
    data = data,
    chains = 1,
    iter = 50,
    cores = 1,
    seed = 456  # Different seed
  )
  
  cat("SUCCESS: Workflow completed!\n")
  print(result)
  
}, error = function(e) {
  cat("ERROR in workflow:\n")
  cat(paste("Error message:", e$message, "\n"))
  
  # Try to get the actual Stan code if possible
  if (grepl("half_cauchy", e$message, ignore.case = TRUE)) {
    cat("\nThis suggests the Stan model still has half_cauchy_lpdf\n")
    cat("Let's check if this is a caching issue...\n")
  }
})