use anyhow::Result;
use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::{LogProb, Prob};
use log::info;
use statrs::distribution::{Discrete, NegativeBinomial};
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
#[derive(Clone)]
pub struct GenomicPrior {
    pub r: f64, // number of successes parameter
    pub p: f64, // success probability parameter
}

impl Prior for GenomicPrior {
    type Event = GenomicEvent;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        let prior_signal = 0.01; // P(signal): Prior probability of true signal
        LogProb::from(Prob(prior_signal))
    }
}

impl GenomicPrior {
    pub fn from_bin_counts(bin_counts: &[(GenomicRange, usize)]) -> Self {
        let counts: Vec<f64> = bin_counts.iter().map(|(_, count)| *count as f64).collect();

        if counts.is_empty() || counts.len() < 2 {
            return Self { r: 2.0, p: 0.3 };
        }

        let mean = counts.clone().mean();
        let variance = counts.clone().variance();

        info!(
            "Count statistics: mean={}, variance={}, n={}",
            mean,
            variance,
            counts.len()
        );

        if variance > mean && variance.is_finite() && mean > 0.0 {
            let r = (mean * mean) / (variance - mean);
            let p = mean / variance;

            if r.is_finite() && r > 0.0 && p.is_finite() && p > 0.0 && p < 1.0 {
                let final_r = r.clamp(0.1, 100.0); // More reasonable bounds
                let final_p = p.clamp(0.05, 0.95); // More reasonable bounds
                info!("Using method of moments: r={}, p={}", final_r, final_p);
                return Self {
                    r: final_r,
                    p: final_p,
                };
            }
        }

        info!("Using conservative fallback parameters");
        Self { r: 2.0, p: 0.3 }
    }
}

pub struct GenomicLikelihood {
    pub r: f64,
    pub p: f64,
}

impl Likelihood for GenomicLikelihood {
    type Event = GenomicEvent;
    type Data = ReadCountData;

    fn compute(&self, _event: &Self::Event, data: &Self::Data, _payload: &mut ()) -> LogProb {
        match NegativeBinomial::new(self.r, self.p) {
            Ok(nb_dist) => {
                let log_likelihood = nb_dist.ln_pmf(data.observed_count as u64);
                info!(
                    "NB likelihood: count={}, r={}, p={}, ln_pmf={}",
                    data.observed_count, self.r, self.p, log_likelihood
                );
                LogProb::from(log_likelihood)
            }
            Err(e) => {
                info!(
                    "NegativeBinomial::new failed: r={}, p={}, error={:?}",
                    self.r, self.p, e
                );
                LogProb::ln_zero()
            }
        }
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
        // Calculate signal likelihood using current parameters
        let signal_likelihood = joint_prob(event, data);

        // Calculate noise likelihood using background parameters (conservative, high p)
        let noise_likelihood = match NegativeBinomial::new(1.0, 0.8) {
            Ok(nb_dist) => {
                let log_likelihood = nb_dist.ln_pmf(data.observed_count as u64);
                LogProb::from(log_likelihood)
            }
            Err(_) => LogProb::ln_zero(),
        };

        let prior_signal = LogProb::from(Prob(0.5));
        let prior_noise = LogProb::from(Prob(0.5));

        let joint_prob_signal = signal_likelihood + prior_signal;
        let joint_prob_noise = noise_likelihood + prior_noise;

        // Calculate evidence (marginal likelihood)
        let evidence = joint_prob_signal.ln_add_exp(joint_prob_noise);

        let posterior = joint_prob_signal - evidence;
        info!("Posterior calculation: signal_lik={:?}, noise_lik={:?}, signal_joint={:?}, noise_joint={:?}, evidence={:?}, posterior={:?}", 
              signal_likelihood, noise_likelihood, joint_prob_signal, joint_prob_noise, evidence, posterior);

        posterior
    }
}

pub struct BayesianModel {
    significance_threshold: f64,
    min_reads: u32,
    model: Model<GenomicLikelihood, GenomicPrior, GenomicPosterior>,
}

impl BayesianModel {
    pub fn new(significance_threshold: f64, min_reads: u32) -> Self {
        let likelihood = GenomicLikelihood { r: 1.0, p: 0.5 };
        let prior = GenomicPrior { r: 1.0, p: 0.5 };
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

        let prior = GenomicPrior::from_bin_counts(bin_counts);
        info!(
            "Estimated negative binomial parameters: r={}, p={}",
            prior.r, prior.p
        );
        *self.model.prior_mut() = prior.clone();

        self.model.likelihood_mut().r = prior.r;
        self.model.likelihood_mut().p = prior.p;

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

            let mut joint_prob_fn = |event: &GenomicEvent, data: &ReadCountData| -> LogProb {
                let likelihood = self.model.likelihood().compute(event, data, &mut ());
                let prior = self.model.prior().compute(event);
                likelihood + prior
            };

            let posterior_log = self
                .model
                .posterior()
                .compute(&event, &data, &mut joint_prob_fn);
            let posterior_prob = (*Prob::from(posterior_log)).clamp(0.0, 1.0);

            info!(
                "Bin {}:{}-{} count={} posterior_log={:?} posterior_prob={}",
                bin.chrom, bin.start, bin.end, count, posterior_log, posterior_prob
            );
            posterior_probs.push((bin.clone(), posterior_prob));
        }

        let significant =
            self.apply_posterior_threshold(posterior_probs, self.significance_threshold)?;

        let significant_bins: Vec<GenomicRange> =
            significant.into_iter().map(|(bin, _)| bin).collect();

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
                info!(
                    "Setting bin {}:{}-{} posterior_prob to: {}",
                    bin_clone.chrom, bin_clone.start, bin_clone.end, bin_clone.posterior_prob
                );
                significant.push((bin_clone, posterior_prob));
            }
        }

        Ok(significant)
    }
}
