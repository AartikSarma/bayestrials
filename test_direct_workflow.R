#!/usr/bin/env Rscript

# Test workflow with direct brms call to isolate the issue

devtools::load_all()
library(brms)

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
  model_name = "test_direct"
)

# Create priors
priors <- PriorSpecification()
priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
priors <- add_prior(priors, "b_baseline_severity_centered", "normal", 0, 2.0)
priors <- add_prior(priors, "Intercept", "normal", 0, 10)
priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

# Convert priors
brms_priors <- to_brms_prior(priors)
print("Converted priors:")
print(brms_priors)

# Build formula  
formula <- build_formula(model_spec)
print("Formula:")
print(formula)

# Test direct brms call
cat("\n=== Testing direct brms call ===\n")

tryCatch({
  model <- brm(
    formula = formula,
    data = data,
    prior = brms_priors,
    family = gaussian(),
    chains = 1,
    iter = 50,
    cores = 1,
    refresh = 0
  )
  
  cat("SUCCESS: Direct brms call worked!\n")
  print(summary(model))
  
}, error = function(e) {
  cat("ERROR in direct brms call:\n")
  cat(paste("Error:", e$message, "\n"))
  
  # Try to get more details
  if(grepl("half_cauchy", e$message)) {
    cat("This confirms the half_cauchy issue still exists\n")
  }
})