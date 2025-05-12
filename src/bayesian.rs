use anyhow::Result;
use log::{debug, info};
use statrs::distribution::{Beta, ContinuousCDF};
use statrs::statistics::{Distribution, Statistics};
use std::collections::HashMap;

use crate::bam::GenomicRange;

pub struct BayesianModel {
    significance_threshold: f64,
    min_reads: u32,
    step: u32,
    slide: u32,
}

impl BayesianModel {
    pub fn new(significance_threshold: f64, min_reads: u32, step: u32, slide: u32) -> Self {
        Self {
            significance_threshold,
            min_reads,
            step,
            slide,
        }
    }
    
    pub fn identify_significant_bins(
        &self,
        bin_counts: &[(GenomicRange, usize)],
        total_reads: usize,
    ) -> Result<Vec<GenomicRange>> {
        info!("Applying Bayesian model to identify significant bins");
        
        let (alpha, beta) = self.calculate_prior_parameters(bin_counts, total_reads);
        
        debug!("Prior parameters: alpha={}, beta={}", alpha, beta);
        
        let prior_distribution = Beta::new(alpha, beta)?;
        
        let mut significant_bins = Vec::new();
        let mut p_values = Vec::new();
        
        for (bin, count) in bin_counts {
            if *count < self.min_reads as usize {
                continue;
            }
            
            let posterior_alpha = alpha + *count as f64;
            let posterior_beta = beta + (total_reads - *count) as f64;
            
            let p_value = self.calculate_bayesian_p_value(
                *count as f64 / total_reads as f64,
                &prior_distribution,
                posterior_alpha,
                posterior_beta,
            )?;
            
            p_values.push((bin.clone(), p_value));
        }
        
        let significant = self.apply_fdr_correction(p_values, self.significance_threshold)?;
        
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
            .filter_map(|(_, count)| if *count > 0 { Some(*count as f64) } else { None })
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
    
    fn calculate_bayesian_p_value(
        &self,
        observed_proportion: f64,
        prior_distribution: &Beta,
        posterior_alpha: f64,
        posterior_beta: f64,
    ) -> Result<f64> {
        let p_value = 1.0 - prior_distribution.cdf(observed_proportion);
        
        let posterior_distribution = Beta::new(posterior_alpha, posterior_beta)?;
        let posterior_mean = posterior_distribution.mean();
        
        let posterior_mean_val = posterior_mean.unwrap_or(0.5);
        let prior_mean_val = prior_distribution.mean().unwrap_or(0.5);
        
        let adjusted_p_value = p_value * (posterior_mean_val / prior_mean_val);
        
        Ok(adjusted_p_value.min(1.0).max(0.0))
    }
    
    fn apply_fdr_correction(
        &self,
        p_values: Vec<(GenomicRange, f64)>,
        threshold: f64,
    ) -> Result<Vec<(GenomicRange, f64)>> {
        let mut sorted_p_values = p_values;
        sorted_p_values.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        
        let n = sorted_p_values.len();
        let mut significant = Vec::new();
        
        for (i, (bin, p_value)) in sorted_p_values.iter().enumerate() {
            let adjusted_p_value = p_value * n as f64 / (i + 1) as f64;
            
            if adjusted_p_value <= threshold {
                let mut bin_clone = bin.clone();
                bin_clone.p_value = *p_value;
                significant.push((bin_clone, *p_value));
            }
        }
        
        Ok(significant)
    }
}
