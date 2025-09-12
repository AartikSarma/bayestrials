# BayesTrials: A Comprehensive Framework for Reproducible Bayesian Clinical Trial Analysis

[![R-CMD-check](https://github.com/AartikSarma/bayestrials/workflows/R-CMD-check/badge.svg)](https://github.com/AartikSarma/bayestrials/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> A comprehensive R package for reproducible Bayesian analyses of clinical datasets using brms

BayesTrials provides a standardized framework for conducting Bayesian clinical trial analyses with systematic prior specification, automated model fitting, and comprehensive reporting capabilities.

## Features

- **Standardized Prior Specification**: Tibble-based format for defining neutral, optimistic, and skeptical priors
- **Multiple Model Support**: Gaussian, Binomial, Poisson, and other distributions via brms
- **Parallel Processing**: Automated fitting of all prior × model permutations using `future`/`furrr`
- **Convergence Diagnostics**: Built-in Rhat, ESS, and divergence monitoring
- **Rich Visualizations**: Prior vs posterior plots, forest plots, and diagnostic dashboards
- **Automated Reporting**: Parameterized RMarkdown templates with HTML/PDF export
- **Model Comparison**: WAIC and LOO-IC for Bayesian model selection
- **CSV-based Workflow**: Define models and priors through simple CSV files
- **Data Harmonization**: Work with multiple trials with different variable names
- **Effect Modifier Identification**: Discover and validate treatment effect heterogeneity

## Installation

### System Requirements

Before installing bayestrials, ensure you have:

- **R** (≥ 4.0.0)
- **C++ compiler** for Stan compilation:
  - **macOS**: Install Xcode Command Line Tools: `xcode-select --install`
  - **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  - **Linux**: Install `g++` and related development tools

### Package Installation

```r
# Install development version from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("AartikSarma/bayestrials")

# Load the package
library(bayestrials)
```

### Key Dependencies

bayestrials uses [brms](https://paul-buerkner.github.io/brms/) (Bayesian Regression Models using Stan) as its modeling backend. This provides:
- Robust Stan model compilation and caching
- Comprehensive family support (Gaussian, Binomial, Poisson, etc.)  
- Advanced model diagnostics and convergence checking

## Use Cases

bayestrials is designed for various clinical trial analysis scenarios:

### 🏥 **Primary Analysis of RCTs**
- Analyze treatment effects with appropriate uncertainty quantification
- Handle missing data with principled Bayesian approaches
- Incorporate historical or expert prior information
- Generate regulatory-ready analysis reports

### 📊 **Secondary and Exploratory Analyses**  
- Subgroup analyses with proper multiplicity adjustment
- Time-to-event analyses with survival models
- Longitudinal analyses with mixed-effects models
- Biomarker and pharmacokinetic modeling

### 🔍 **Meta-Analysis and Evidence Synthesis**
- Bayesian meta-analysis of multiple trials
- Network meta-analysis for indirect comparisons
- Integration of real-world evidence with trial data
- Historical borrowing with dynamic priors

### 🎯 **Regulatory Submissions**
- FDA/EMA compliant Bayesian analysis plans
- Sensitivity analyses for regulatory review  
- Prior justification and robustness assessment
- Comprehensive documentation and reproducibility

## Quick Start

```r
library(bayestrials)

# 1. Generate synthetic clinical trial data
set.seed(123)
data <- generate_synthetic_data(
  n_observations = 200,    # Start smaller for quicker testing
  outcome_type = "continuous",
  seed = 123
)

# 2. Create model specification
model_spec <- create_model_spec(
  outcome = "outcome_continuous",
  predictors = "treatment + age_centered + baseline_severity_centered",
  family = "gaussian",
  model_name = "linear_regression"
)

# 3. Specify priors (using the improved PriorSpecification system)
priors <- PriorSpecification()
priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
priors <- add_prior(priors, "b_baseline_severity_centered", "normal", 0, 2.0)
priors <- add_prior(priors, "Intercept", "normal", 0, 10)
priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)  # Automatically converts to student_t

# Attach priors to model
attr(model_spec, "priors") <- priors

# 4. Set up parallel processing (optional)
setup_parallel(strategy = "sequential")  # Use sequential for initial testing

# 5. Fit Bayesian model
result <- fit_model(
  model_spec = model_spec,
  data = data,
  chains = 2,      # Fewer chains for quicker testing
  iter = 1000      # Fewer iterations for quicker testing
)

# 6. Check model and diagnostics
print(result)
diagnostics <- check_diagnostics(result)
print(diagnostics)

# 7. Extract and visualize results
estimates <- extract_estimates(result)
print(estimates)

# Create visualizations
plot_posterior(result)
forest_plot <- plot_forest(estimates)
```

### Expected Output

When the above code runs successfully, you should see:

```r
# Model summary output:
Bayesian Clinical Trial Model
----------------------------
Model Specification: linear_regression 
  Formula: outcome_continuous ~ treatment + age_centered + baseline_severity_centered
  Family: gaussian
  Link: identity

Model Summary:
  Algorithm: sampling
  Samples: 1000
  Chains: 2

Treatment Effect:
  Estimate: 1.256
  95% CI: [0.892, 1.634]

Diagnostics:
  Warnings: 0
  Rhat issues: 0 parameters
```

**Key indicators of success:**
- **Rhat values ≈ 1.00**: All parameters converged properly
- **ESS > 400**: Adequate effective sample size for reliable estimates  
- **No divergences**: MCMC sampling was efficient
- **Treatment effect**: Credible interval and point estimate

## CSV-Based Workflow

For larger projects and reproducible analyses, bayestrials supports CSV-based model and prior specification:

### 📁 **File Structure**
```
analysis_project/
├── data/
│   ├── trial_data.csv           # Your clinical trial data
│   └── covariate_specs.csv      # Variable definitions and transformations
├── specifications/
│   ├── model_specs.csv          # Model specifications  
│   └── prior_specs.csv          # Prior specifications
└── analysis/
    └── run_analysis.R           # Main analysis script
```

### 📋 **Model Specification CSV Format**
```csv
model_name,outcome,predictors,family,link
primary_analysis,response,treatment + age + sex + baseline_score,gaussian,identity
sensitivity_1,response,treatment + age + sex,gaussian,identity  
subgroup_analysis,response,treatment * age + sex + baseline_score,gaussian,identity
```

### 📊 **Prior Specification CSV Format**
```csv
parameter,distribution,location,scale,lb,ub
b_treatment,normal,0,2.5,NA,NA
b_age,normal,0,1.0,NA,NA
b_sex,normal,0,1.5,NA,NA
Intercept,normal,50,20,NA,NA
sigma,half_cauchy,NA,5,0,NA
```

### 🚀 **Complete Workflow Example**

```r
# Load and prepare data with comprehensive preprocessing
data <- load_clinical_trial("data/trial_data.csv")
data <- preprocess_data(data, specs = "data/covariate_specs.csv")

# Print data summary
cat("Dataset loaded successfully:\n")
cat("- Observations:", nrow(data), "\n")
cat("- Variables:", ncol(data), "\n")
cat("- Treatment groups:", table(data$treatment), "\n")

# Read all model specifications from CSV
model_specs <- read_model_spec(
  model_file = "specifications/model_specs.csv",
  prior_file = "specifications/prior_specs.csv",
  validate = TRUE  # Validate all specifications
)

# Display available models
cat("\nAvailable models:\n")
for(i in seq_along(model_specs)) {
  cat("-", names(model_specs)[i], "\n")
}

# Fit primary model
primary_model <- fit_model(
  model_spec = model_specs$primary_analysis,
  data = data,
  chains = 4,
  iter = 2000,
  seed = 123
)

# Generate comprehensive analysis report
report_file <- generate_report(
  model = primary_model,
  template = "basic_report",
  output_file = "primary_analysis_report.html",
  params = list(
    title = "Primary Efficacy Analysis",
    subtitle = "Bayesian Analysis of Treatment Effect"
  )
)

cat("Analysis complete! Report saved to:", report_file, "\n")
```

### 📈 **Batch Processing Multiple Models**

```r
# Fit all specified models in parallel
all_models <- fit_models(
  model_specs = model_specs,
  data = data,
  parallel = TRUE,
  cores = 4
)

# Compare models using information criteria
model_comparison <- compare_models(all_models)
print(model_comparison)

# Generate comparison report
comparison_report <- generate_report(
  models = all_models,
  template = "model_comparison",
  output_file = "model_comparison_report.html"
)
```

## Multi-Trial Workflow

```r
# Load multiple trials
trials <- load_trials(c("trial1.csv", "trial2.csv", "trial3.csv"))

# Harmonize variables across trials
harmonized_trials <- harmonize_trials(trials, mapping = "variable_mapping.csv")

# Prepare data
prep_data <- preprocess_data(harmonized_trials, specs = "covariate_specs.csv")

# Specify and fit hierarchical model
model_spec <- read_model_spec("multi_trial_model.csv", 
                            "multi_trial_priors.csv")
model <- fit_model(model_spec, prep_data)

# Analyze treatment effects across trials
trial_effects <- extract_trial_effects(model, "treatment")
plot_forest(trial_effects)

# Test for effect modifiers
modifiers <- identify_effect_modifiers(model, prep_data, 
                                    treatment_var = "treatment",
                                    candidate_modifiers = c("age", "sex", "biomarker"))
```

## Sensitivity Analysis Workflow

```r
# Specify base model
base_spec <- read_model_spec("base_model.csv", "base_priors.csv")
base_model <- fit_model(base_spec, data)

# Define sensitivity analyses
sensitivity_specs <- list(
  prior_skeptical = update_priors(base_spec, skepticism = "high"),
  prior_optimistic = update_priors(base_spec, skepticism = "low"),
  missing_mar = update_model(base_spec, missing_mechanism = "mar"),
  missing_mnar = update_model(base_spec, missing_mechanism = "mnar")
)

# Run sensitivity analyses
sensitivity_results <- run_sensitivity_analyses(base_model, 
                                              data, 
                                              sensitivity_specs)

# Visualize sensitivity
plot_sensitivity(sensitivity_results, 
                parameter = "treatment_effect")

# Generate summary report
sensitivity_report(sensitivity_results, 
                  template = "sensitivity_report")
```

## Package Structure

```
bayestrials/
├── R/
│   ├── model_spec.R             # Model specification framework
│   ├── prior_spec.R             # Prior specification system  
│   ├── model_fitting.R          # brms-based model fitting
│   ├── visualization.R          # Plotting functions
│   ├── diagnostics.R            # Convergence diagnostics
│   ├── generate_synthetic_data.R # Test data generation
│   ├── standard_priors.R        # Standard prior sets
│   ├── enhanced_diagnostics.R   # Enhanced diagnostic tools
│   └── parallel.R               # Parallel processing utilities
├── inst/
│   ├── extdata/                 # Example CSV files
│   └── rmarkdown/               # Report templates
├── tests/
│   └── testthat/                # Unit tests
└── vignettes/                   # Documentation
```

## Dependencies

### Core Statistical
- [brms](https://paul-buerkner.github.io/brms/) - Bayesian regression models using Stan
- [posterior](https://mc-stan.org/posterior/) - Tools for working with posterior distributions
- [bridgesampling](https://github.com/quentingronau/bridgesampling) - Model comparison via bridge sampling

### Data Management
- [tidyverse](https://www.tidyverse.org/) - Data manipulation and visualization
- [tibble](https://tibble.tidyverse.org/) - Modern data frames
- [dplyr](https://dplyr.tidyverse.org/) - Data manipulation
- [tidyr](https://tidyr.tidyverse.org/) - Data tidying

### Parallel Processing  
- [future](https://future.futureverse.org/) - Parallel and distributed processing
- [furrr](https://furrr.futureverse.org/) - Apply mapping functions in parallel
- [progressr](https://progressr.futureverse.org/) - Progress reporting

### Visualization & Reporting
- [ggplot2](https://ggplot2.tidyverse.org/) - Grammar of graphics
- [plotly](https://plotly.com/r/) - Interactive visualizations
- [rmarkdown](https://rmarkdown.rstudio.com/) - Dynamic documents
- [gt](https://gt.rstudio.com/) - Publication-ready tables

## Development

### Setup
```r
# Install development dependencies
devtools::install_dev_deps()

# Load package for development
devtools::load_all()

# Generate documentation
devtools::document()
```

### Testing
```r
# Run all tests
devtools::test()

# Test the complete workflow
source("test_bayestrials_workflow.R")
```

### Package Checks
```r
# Run R CMD check
devtools::check()
```

## Troubleshooting

### 🔧 **Installation Issues**

#### **C++ Compiler Problems**

**macOS - Xcode Command Line Tools**
```bash
# Check if tools are installed
xcode-select -p

# Install if missing
xcode-select --install

# Verify installation
gcc --version
make --version
```

**Windows - Rtools Installation**
```r
# Check current R version
R.version.string

# Download appropriate Rtools from https://cran.r-project.org/bin/windows/Rtools/
# For R 4.3+: Rtools43
# For R 4.2: Rtools42

# Verify Rtools installation
Sys.which("make")
Sys.which("gcc")

# If empty, add to PATH or reinstall Rtools with "Add to PATH" checked
```

**Linux - Build Tools**
```bash
# Ubuntu/Debian
sudo apt update
sudo apt install build-essential r-base-dev

# CentOS/RHEL/Fedora  
sudo yum groupinstall "Development Tools"
sudo yum install R-devel

# Verify installation
gcc --version
make --version
```

#### **R Package Dependencies**

```r
# Check for missing dependencies
required_packages <- c("brms", "rstan", "rstantools", "bridgesampling", 
                      "tidyverse", "future", "furrr", "progressr")

missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if(length(missing) > 0) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing)
}

# For brms specifically (if installation fails)
install.packages("brms", dependencies = TRUE)

# For development version of brms (if needed)
remotes::install_github("paul-buerkner/brms")
```

#### **Stan Configuration Issues**

```r
# Configure Stan compilation options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Test Stan installation
library(rstan)
example(stan_model, package = "rstan", run.dontrun = TRUE)

# If Stan fails, try rebuilding
remove.packages("rstan")
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

### 🎯 **Model Fitting Issues**

#### **Convergence Diagnostics**

**High Rhat Values (> 1.01)**
```r
# Increase iterations and warmup
result <- fit_model(
  model_spec, data,
  iter = 4000,        # Double iterations
  warmup = 2000,      # Increase warmup
  chains = 4
)

# Increase adaptation
result <- fit_model(
  model_spec, data,
  control = list(
    adapt_delta = 0.95,    # Increase from default 0.8
    max_treedepth = 12     # Increase from default 10
  )
)

# Check specific parameters with issues
diagnostics <- check_diagnostics(result)
print(diagnostics$rhat_summary)
```

**Divergent Transitions**
```r
# Increase adapt_delta (most common solution)
result <- fit_model(
  model_spec, data,
  control = list(adapt_delta = 0.99)  # Up to 0.999 if needed
)

# Reparameterize model if persistent divergences
# Check for highly correlated parameters
estimates <- extract_estimates(result)
cor_matrix <- cor(extract_posterior(result))
print(cor_matrix)

# Consider different prior specifications
# Hierarchical centering vs non-centering
```

**Low Effective Sample Size (ESS)**
```r
# Increase total iterations
result <- fit_model(
  model_spec, data,
  iter = 8000,        # More iterations
  chains = 4
)

# Check for poor mixing - examine trace plots
plot_diagnostics(result, type = "trace")

# Consider thinning (last resort)
result <- fit_model(
  model_spec, data,
  thin = 2           # Keep every 2nd sample
)
```

#### **Prior Specification Issues**

**Understanding Prior Conversion**
```r
# Check how your priors are converted
priors <- PriorSpecification()
priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

# Preview conversion
brms_priors <- to_brms_prior(priors)
print(brms_priors)  # Shows: student_t(3, 0, 5)

# Visualize prior implications
plot_priors(priors)
```

**Prior-Data Conflicts**
```r
# Check for unreasonable priors
summary(your_data$outcome)  # Check outcome range
# If outcome ranges 0-100, don't use Normal(0, 1000) for intercept

# Use prior predictive checks
prior_predictive <- brm(
  formula, data,
  prior = your_priors,
  sample_prior = "only",  # Sample from prior only
  chains = 2, iter = 500
)

# Check if prior predictions make sense
plot(prior_predictive)
```

#### **Data-Related Issues**

**Missing Data**
```r
# Check missing data patterns
library(VIM)
VIM::aggr(data, col = c('navyblue','red'), numbers = TRUE, sortVars = TRUE)

# Handle missing data explicitly
complete_data <- na.omit(data)  # Complete case analysis
# Or use multiple imputation approaches

# Model can handle some missingness
result <- fit_model(model_spec, data)  # brms handles NAs in predictors
```

**Scaling and Centering Issues**
```r
# Check variable scales
summary(data[, c("age", "baseline_score", "biomarker")])

# Standardize continuous variables
data$age_scaled <- scale(data$age)[, 1]
data$baseline_scaled <- scale(data$baseline_score)[, 1]

# Update model specification
model_spec <- create_model_spec(
  outcome = "response",
  predictors = "treatment + age_scaled + baseline_scaled"
)
```

**Extreme Outliers**
```r
# Identify outliers
boxplot(data$outcome)
outliers <- which(abs(scale(data$outcome)) > 3)

# Options:
# 1. Remove outliers (with justification)
clean_data <- data[-outliers, ]

# 2. Use robust models (Student-t errors instead of normal)
robust_spec <- create_model_spec(
  outcome = "response",
  predictors = "treatment + covariates",
  family = "student_t"  # Heavy-tailed errors
)

# 3. Transform outcome
data$log_outcome <- log(data$outcome + 1)  # Log transform
```

### 🚀 **Performance Optimization**

**Slow Model Fitting**
```r
# Use more cores
result <- fit_model(
  model_spec, data,
  cores = parallel::detectCores() - 1  # Leave one core free
)

# Reduce complexity for initial fitting
small_data <- data[sample(nrow(data), 100), ]  # Subset for testing
simple_spec <- create_model_spec(
  outcome = "response", 
  predictors = "treatment"  # Minimal model first
)

# Optimize Stan settings
result <- fit_model(
  model_spec, data,
  algorithm = "meanfield"  # Faster variational inference (less accurate)
)
```

**Memory Issues**
```r
# Monitor memory usage
gc()  # Garbage collection
object.size(result)  # Check model size

# Reduce memory footprint
result <- fit_model(
  model_spec, data,
  save_warmup = FALSE,     # Don't save warmup samples
  save_all_pars = FALSE    # Don't save all parameters
)

# Clear workspace regularly
rm(large_objects)
gc()
```

### 📋 **Diagnostic Checklist**

Before reporting issues, check:

```r
# 1. System information
sessionInfo()
cat("bayestrials version:", packageVersion("bayestrials"), "\n")

# 2. Data integrity  
str(data)
summary(data)
sum(is.na(data))

# 3. Model specification
print(model_spec)
validate_model_spec(model_spec)

# 4. Prior specifications
print(priors)
brms_version <- to_brms_prior(priors)
print(brms_version)

# 5. Fitting diagnostics
diagnostics <- check_diagnostics(result)
print(diagnostics)
```

## Getting Help

- **Vignettes**: Run `browseVignettes("bayestrials")` for detailed tutorials
- **Function Help**: Use `?function_name` for specific function documentation  
- **Issues**: Report bugs and feature requests on [GitHub](https://github.com/AartikSarma/bayestrials/issues)
- **Stan Resources**: [Stan User's Guide](https://mc-stan.org/users/documentation/) for underlying methodology

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use `bayestrials` in your research, please cite:

```bibtex
@software{bayestrials,
  author = {Sarma, Aartik},
  title = {BayesTrials: A Comprehensive Framework for Reproducible Bayesian Clinical Trial Analysis},
  url = {https://github.com/AartikSarma/bayestrials},
  version = {0.1.0},
  year = {2024}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.