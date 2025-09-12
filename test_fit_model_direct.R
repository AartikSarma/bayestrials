#!/usr/bin/env Rscript

# Test fit_model function directly

devtools::load_all()

# Generate data
set.seed(123) 
data <- generate_synthetic_data(
  n_observations = 100,
  n_variables = 5,
  outcome_type = "continuous",
  seed = 123
)

# Create model specification  
model_spec <- create_model_spec(
  outcome = "outcome_continuous",
  predictors = "treatment + age_centered + baseline_severity_centered",
  family = "gaussian",
  model_name = "test_fit_model"
)

# Create priors
priors <- PriorSpecification()
priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
priors <- add_prior(priors, "b_baseline_severity_centered", "normal", 0, 2.0)
priors <- add_prior(priors, "Intercept", "normal", 0, 10)
priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

# Attach priors to model spec
attr(model_spec, "priors") <- priors

cat("=== Testing fit_model function ===\n")
cat("Model spec created with priors attached\n")

tryCatch({
  result <- fit_model(
    model_spec = model_spec,
    data = data,
    chains = 1,
    iter = 50,
    cores = 1,
    seed = 123
  )
  
  cat("SUCCESS: fit_model function worked!\n")
  print(result)
  
}, error = function(e) {
  cat("ERROR in fit_model function:\n")
  cat(paste("Error:", e$message, "\n"))
})