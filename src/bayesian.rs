use anyhow::Result;
use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::{LogProb, Prob};
use log::info;
use statrs::statistics::Statistics;

use crate::bam::GenomicRange;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GenomicEvent {
    pub bin_count: usize,
    pub total_reads: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ReadCountData {
    pub observed_count: usize,
}

#[allow(dead_code)]
pub struct GenomicPrior {
    pub alpha: f64,
    pub beta: f64,
}

impl Prior for GenomicPrior {
    type Event = GenomicEvent;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        let prior_signal = 0.01; // P(signal): Prior probability of true signal
        LogProb::from(Prob(prior_signal))
    }
}

impl GenomicPrior {
    pub fn from_bin_counts(bin_counts: &[(GenomicRange, usize)], total_reads: usize) -> Self {
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
            return Self {
                alpha: 1.0,
                beta: 99.0,
            };
        }

        let alpha_beta_sum = total_mean * (1.0 - total_mean) / total_variance - 1.0;
        let alpha = total_mean * alpha_beta_sum;
        let beta = (1.0 - total_mean) * alpha_beta_sum;

        Self {
            alpha: alpha.max(0.01),
            beta: beta.max(0.01),
        }
    }
}

pub struct GenomicLikelihood;

impl Likelihood for GenomicLikelihood {
    type Event = GenomicEvent;
    type Data = ReadCountData;

    fn compute(&self, event: &Self::Event, data: &Self::Data, _payload: &mut ()) -> LogProb {
        let observed_proportion = data.observed_count as f64 / event.total_reads as f64;

        let likelihood_given_signal = observed_proportion.max(0.5);

        LogProb::from(Prob(likelihood_given_signal))
    }
}

pub struct GenomicPosterior;

impl Posterior for GenomicPosterior {
    type Event = GenomicEvent;
    type BaseEvent = GenomicEvent;
    type Data = ReadCountData;

    fn compute<F>(&self, event: &Self::Event, data: &Self::Data, joint_prob: &mut F) -> LogProb
    where
        F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb,
    {
        // Calculate joint probability for signal hypothesis
        let joint_prob_signal = joint_prob(event, data);

        let noise_event = GenomicEvent {
            bin_count: (data.observed_count as f64 * 0.1) as usize,
            total_reads: event.total_reads,
        };

        // Calculate joint probability for noise hypothesis
        let joint_prob_noise = joint_prob(&noise_event, data);

        // Calculate evidence (marginal likelihood)
        let evidence = joint_prob_signal.ln_add_exp(joint_prob_noise);

        joint_prob_signal - evidence
    }
}

pub struct BayesianModel {
    significance_threshold: f64,
    min_reads: u32,
    model: Model<GenomicLikelihood, GenomicPrior, GenomicPosterior>,
}

impl BayesianModel {
    pub fn new(significance_threshold: f64, min_reads: u32) -> Self {
        let likelihood = GenomicLikelihood;
        let prior = GenomicPrior {
            alpha: 1.0,
            beta: 99.0,
        };
        let posterior = GenomicPosterior;
        let model = Model::new(likelihood, prior, posterior);

        Self {
            significance_threshold,
            min_reads,
            model,
        }
    }

    pub fn significance_threshold(&self) -> f64 {
        self.significance_threshold
    }

    pub fn min_reads(&self) -> u32 {
        self.min_reads
    }

    pub fn identify_significant_bins(
        &mut self, // Note: needs to be mutable to update the model
        bin_counts: &[(GenomicRange, usize)],
        total_reads: usize,
    ) -> Result<Vec<GenomicRange>> {
        info!("Applying rust-bio Bayesian model to identify significant bins");

        let prior = GenomicPrior::from_bin_counts(bin_counts, total_reads);
        *self.model.prior_mut() = prior;

        let mut significant_bins = Vec::new();
        let mut posterior_probs = Vec::new();

        for (bin, count) in bin_counts {
            if *count < self.min_reads as usize {
                continue;
            }

            let event = GenomicEvent {
                bin_count: *count,
                total_reads,
            };

            let data = ReadCountData {
                observed_count: *count,
            };

            let universe = vec![event.clone()];
            let model_instance = self.model.compute(universe, &data);

            let posterior_log = model_instance
                .posterior(&event)
                .unwrap_or(LogProb::ln_zero());
            let posterior_prob = (*Prob::from(posterior_log)).clamp(0.0, 1.0);

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
