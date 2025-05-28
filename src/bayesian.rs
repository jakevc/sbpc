use anyhow::Result;
use bio::stats::{LogProb, Prob};
use log::{debug, info};
use statrs::distribution::Beta;
use statrs::statistics::{Distribution, Statistics};

use crate::bam::GenomicRange;

pub struct BayesianModel {
    significance_threshold: f64,
    min_reads: u32,
}

impl BayesianModel {
    pub fn new(significance_threshold: f64, min_reads: u32) -> Self {
        Self {
            significance_threshold,
            min_reads,
        }
    }

    pub fn identify_significant_bins(
        &self,
        bin_counts: &[(GenomicRange, usize)],
        total_reads: usize,
    ) -> Result<Vec<GenomicRange>> {
        info!("Applying Bayesian model to identify significant bins using log-space probabilities");

        let (alpha, beta) = self.calculate_prior_parameters(bin_counts, total_reads);

        debug!("Prior parameters: alpha={}, beta={}", alpha, beta);

        let mut significant_bins = Vec::new();
        let mut posterior_probs = Vec::new();

        for (bin, count) in bin_counts {
            if *count < self.min_reads as usize {
                continue;
            }

            let posterior_alpha = alpha + *count as f64;
            let posterior_beta = beta + (total_reads - *count) as f64;

            let posterior_prob =
                self.calculate_posterior_probability(posterior_alpha, posterior_beta)?;

            posterior_probs.push((bin.clone(), posterior_prob));
        }

        let significant =
            self.apply_posterior_threshold(posterior_probs, self.significance_threshold)?;

        for (bin, _) in significant {
            significant_bins.push(bin);
        }

        info!("Found {} significant bins", significant_bins.len());

        Ok(significant_bins)
    }

    fn calculate_prior_parameters(
        &self,
        bin_counts: &[(GenomicRange, usize)],
        total_reads: usize,
    ) -> (f64, f64) {
        let non_zero_counts: Vec<f64> = bin_counts
            .iter()
            .filter_map(|(_, count)| {
                if *count > 0 {
                    Some(*count as f64)
                } else {
                    None
                }
            })
            .collect();

        let mean = non_zero_counts.clone().mean();
        let variance = non_zero_counts.variance();

        let total_mean = mean / total_reads as f64;
        let total_variance = variance / (total_reads * total_reads) as f64;

        if total_variance >= total_mean * (1.0 - total_mean) {
            return (1.0, 99.0);
        }

        let alpha_beta_sum = total_mean * (1.0 - total_mean) / total_variance - 1.0;
        let alpha = total_mean * alpha_beta_sum;
        let beta = (1.0 - total_mean) * alpha_beta_sum;

        (alpha.max(0.01), beta.max(0.01))
    }

    fn calculate_posterior_probability(
        &self,
        posterior_alpha: f64,
        posterior_beta: f64,
    ) -> Result<f64> {
        let posterior_distribution = Beta::new(posterior_alpha, posterior_beta)?;

        let signal_prob = posterior_distribution.mean().unwrap_or(0.5);

        let log_prob = LogProb::from(Prob(signal_prob));

        Ok((*Prob::from(log_prob)).clamp(0.0, 1.0))
    }

    fn apply_posterior_threshold(
        &self,
        posterior_probs: Vec<(GenomicRange, f64)>,
        threshold: f64,
    ) -> Result<Vec<(GenomicRange, f64)>> {
        let mut significant = Vec::new();

        for (bin, posterior_prob) in posterior_probs {
            if posterior_prob >= threshold {
                let mut bin_clone = bin.clone();
                bin_clone.posterior_prob = posterior_prob;
                significant.push((bin_clone, posterior_prob));
            }
        }

        Ok(significant)
    }
}
