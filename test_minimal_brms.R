#!/usr/bin/env Rscript

# Minimal test to debug brms prior issue
library(brms)
library(rstan)

# Generate minimal data
set.seed(123)
data <- data.frame(
  y = rnorm(20),
  x = rnorm(20)
)

# Create priors including coefficients and sigma
my_prior <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "Intercept"), 
  prior(student_t(3, 0, 5), class = "sigma")
)
print("Priors specified:")
print(my_prior)

cat("\n=== Testing brms model with student_t prior ===\n")

# Try to fit a simple model
tryCatch({
  model <- brm(
    formula = y ~ x,
    data = data,
    prior = my_prior,
    chains = 1,
    iter = 100,
    cores = 1,
    refresh = 0,  # Suppress iteration output
    sample_prior = "only",  # Sample from prior only for speed
    backend = "rstan"
  )
  
  cat("SUCCESS: Model fitted successfully!\n")
  print(model)
  
}, error = function(e) {
  cat("ERROR in model fitting:\n")
  cat(paste("Error:", e$message, "\n"))
})