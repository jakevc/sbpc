# Changelog

All notable changes to SBPC will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2025-05-28
### Changed
- **BREAKING**: Replaced p-values with posterior probabilities in the Bayesian model
- Added rust-bio crate dependency for Bayesian statistical modeling
- Changed CLI parameter from `--pval` to `--posterior-threshold` with default of 0.95
- Updated peak merging logic to multiply posterior probabilities for uncertainty-aware detection
- Modified BED output format to display posterior probabilities
- Enabled local FDR filtering with direct posterior probability thresholds
- Added support for comparing experiments by multiplying posteriors across samples

## [0.2.0] - 2025-05-21
### Changed
- **BREAKING**: Removed `-o/--prefix` option; peaks are now printed to stdout
- **BREAKING**: Changed default step size from 100 to 200
- Added `--metrics-file` flag to optionally write metrics to a file
- Metrics are now printed as JSON to stdout by default
- Simplified BayesianModel constructor by removing the step parameter
- Improved code organization and removed unused code
- Removed utils.rs module
- Updated documentation to reflect new command-line usage

## [0.1.3] - 2025-05-19
### Changed
- Updated installation instructions in README

## [0.1.1] - 2025-05-13
### Changed
- Updated release workflow to only trigger on semver tags pushed to main branch

## [0.1.0] - 2025-05-12

### Added
- Initial release of SBPC (Simple Bayesian Peak Caller)
- Bayesian statistical framework for peak detection
- P-value calculation for each detected peak region
- Memory-efficient implementation
- Support for CUT&TAG/CUT&RUN data
- Command-line interface with various configuration options
- Output in BED format with p-values in the score column
- Performance metrics in JSON format
- Comprehensive test suite
- Synthetic data generation for testing
