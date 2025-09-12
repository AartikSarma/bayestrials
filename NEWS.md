# bayestrials 0.2.0 (Development Version)

## Major Changes

### Migration from rstan to brms Backend

* **BREAKING CHANGE**: bayestrials now uses [brms](https://paul-buerkner.github.io/brms/) as its Bayesian modeling backend instead of direct rstan interface
* This change provides:
  - More robust Stan model compilation and caching
  - Better convergence diagnostics and warnings
  - Expanded family support and link functions
  - Improved error handling and user feedback

### Prior Specification Improvements

* **Automatic Prior Conversion**: `half_cauchy` priors are now automatically converted to `student_t(3, 0, scale)` format for brms compatibility
* **Enhanced Prior Validation**: Improved error checking and validation for prior specifications
* **Backward Compatibility**: Existing prior specifications continue to work with automatic conversion

### New Features

* **Enhanced Diagnostics**: Improved convergence checking with clearer warning messages
* **Better Error Messages**: More informative error messages for common issues
* **Comprehensive Testing**: Added extensive test suite with 69+ tests covering core functionality
* **Improved Documentation**: Updated README and vignettes with working examples

## Bug Fixes

* Fixed non-standard evaluation (NSE) issues in `to_brms_prior()` function (#789b2b1)
* Resolved Stan compilation errors with half-Cauchy priors
* Fixed duplicated prior specifications warnings
* Corrected function name inconsistencies in documentation

## Documentation Updates

* **README**: Updated with working examples, installation guide, and troubleshooting section
* **Vignettes**: Intro vignette now has executable code with real output
* **Function Help**: Improved documentation with practical examples
* **Installation Guide**: Added system requirements and dependency information

## Under the Hood

* Migrated from `rstan::stan()` to `brms::brm()` for model fitting
* Updated prior specification system to handle brms format requirements
* Improved parallel processing configuration
* Enhanced error handling and logging system

## Breaking Changes

While we've maintained backward compatibility for most user-facing functions, some internal changes may affect advanced users:

* Prior conversion now uses brms format internally
* Some model diagnostics output formats have changed
* Error message formats have been updated

## Migration Guide

For users upgrading from previous versions:

1. **Installation**: Ensure you have the required C++ compiler tools
2. **Priors**: Your existing prior specifications will work unchanged
3. **Models**: Model fitting syntax remains the same
4. **Outputs**: BayesianModel objects maintain the same structure

## Thanks

This release includes contributions and testing from the development team, with special recognition for resolving the challenging prior specification system migration.

---

# bayestrials 0.1.0

## Initial Release

* Initial package structure and core functionality
* Basic Bayesian model fitting with rstan backend
* Prior specification system
* Parallel processing support
* Visualization and diagnostics functions