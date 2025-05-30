use crate::bam::{BamProcessor, GenomicRange};
use crate::bayesian::BayesianModel;
use crate::cli::Cli;
use crate::genome::Genome;
use anyhow::Result;
use bio::io::bed::{Record, Writer};

use log::info;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io;

pub struct PeakCaller {
    cli: Cli,
    genome: Genome,
    bam_processor: BamProcessor,
    bayesian_model: BayesianModel,
}

impl PeakCaller {
    pub fn new(cli: &Cli) -> Result<Self> {
        let mut cli_copy = cli.clone();
        cli_copy.adjust_for_broad_peaks();

        let genome = if let Some(chromsize_path) = &cli.chromsize {
            Genome::from_file(chromsize_path)?
        } else {
            Genome::from_bam(&cli.bam)?
        };

        let bam_processor = BamProcessor::new(&cli.bam, cli.control.as_deref())?;

        let bayesian_model = BayesianModel::new(cli.posterior_threshold, cli.minreads);

        Ok(Self {
            cli: cli_copy,
            genome,
            bam_processor,
            bayesian_model,
        })
    }

    pub fn call_peaks(&mut self) -> Result<Peaks> {
        info!("Starting peak calling with Bayesian framework");

        // Parallelize per chromosome
        let chroms: Vec<String> = self.genome.seqnames.clone();
        let step = self.cli.step;
        let mdist = self.cli.mdist;
        let minwidth = self.cli.minwidth;
        let bam_processor = &self.bam_processor;
        let total_reads = bam_processor.total_reads();

        // Create a thread-local clone of the bayesian model for parallel processing
        let model_params = (
            self.bayesian_model.significance_threshold(),
            self.bayesian_model.min_reads(),
        );

        let all_peaks: Vec<GenomicRange> = chroms
            .par_iter()
            .flat_map(|chrom| {
                // Create bins for this chromosome only
                let bins: Vec<GenomicRange> = self.genome.create_bins(step, Some(chrom)).unwrap();

                if bins.is_empty() {
                    return Vec::new();
                }

                // Count reads in bins for this chromosome
                let bin_counts = bam_processor.count_reads_in_bins(&bins).unwrap();

                // Create a thread-local model instance
                let mut thread_local_model = BayesianModel::new(model_params.0, model_params.1);

                // Identify significant bins using the thread-local model
                let significant_bins = thread_local_model
                    .identify_significant_bins(&bin_counts, total_reads)
                    .unwrap();

                // Merge bins into peaks for this chromosome
                let merged_peaks = self.merge_bins_into_peaks(significant_bins, mdist).unwrap();

                // Filter peaks by width and return directly
                self.filter_peaks_by_width(merged_peaks, minwidth).unwrap()
            })
            .collect();

        info!("Peak calling completed, found {} peaks", all_peaks.len());

        Ok(Peaks { ranges: all_peaks })
    }

    fn merge_bins_into_peaks(
        &self,
        bins: Vec<GenomicRange>,
        max_distance: u32,
    ) -> Result<Vec<GenomicRange>> {
        let mut chrom_bins: HashMap<String, Vec<GenomicRange>> = HashMap::new();

        for bin in bins {
            chrom_bins.entry(bin.chrom.clone()).or_default().push(bin);
        }

        let merged: Vec<GenomicRange> = chrom_bins
            .par_iter()
            .flat_map(|(_chrom, bins)| {
                let mut sorted_bins = bins.clone();
                sorted_bins.sort_by_key(|b| b.start);

                let mut result = Vec::new();
                if sorted_bins.is_empty() {
                    return result;
                }

                let mut current_peak = sorted_bins[0].clone();

                for bin in sorted_bins.iter().skip(1) {
                    if bin.start <= current_peak.end + max_distance {
                        current_peak.end = bin.end.max(current_peak.end);

                        current_peak.posterior_prob = current_peak.posterior_prob.max(bin.posterior_prob);
                    } else {
                        result.push(current_peak.clone());
                        current_peak = bin.clone();
                    }
                }

                result.push(current_peak);
                result
            })
            .collect();

        Ok(merged)
    }

    fn filter_peaks_by_width(
        &self,
        peaks: Vec<GenomicRange>,
        min_width: u32,
    ) -> Result<Vec<GenomicRange>> {
        let filtered = peaks
            .into_iter()
            .filter(|peak| peak.end - peak.start >= min_width)
            .collect();

        Ok(filtered)
    }
}

pub struct Peaks {
    pub ranges: Vec<GenomicRange>,
}

impl Peaks {
    pub fn write_to_stdout_bed(&self) -> usize {
        let stdout = io::stdout();
        let handle = stdout.lock();
        let mut writer = Writer::new(handle);

        for (i, range) in self.ranges.iter().enumerate() {
            let mut record = Record::new();
            record.set_chrom(&range.chrom);
            record.set_start(range.start as u64);
            record.set_end(range.end as u64);
            record.set_name(&format!("peak{}", i + 1));
            record.set_score(&format!("{:.6}", range.posterior_prob));
            record.push_aux("."); // Add strand field as "." (unknown)

            writer.write(&record).unwrap();
        }
        self.ranges.len()
    }
}
