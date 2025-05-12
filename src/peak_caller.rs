use crate::bam::{BamProcessor, GenomicRange};
use crate::bayesian::BayesianModel;
use crate::cli::Cli;
use crate::genome::Genome;
use anyhow::{Context, Result};
use log::{debug, info};
use rayon::prelude::*;
use std::collections::HashMap;

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
        
        let bayesian_model = BayesianModel::new(
            cli.pval,
            cli.minreads,
            cli.step,
            cli.slide,
        );
        
        Ok(Self {
            cli: cli_copy,
            genome,
            bam_processor,
            bayesian_model,
        })
    }
    
    pub fn call_peaks(&mut self) -> Result<Peaks> {
        info!("Starting peak calling with Bayesian framework");
        
        let bins = self.genome.create_bins(self.cli.step, self.cli.slide)?;
        
        let bin_counts = self.bam_processor.count_reads_in_bins(&bins)?;
        
        let significant_bins = self.bayesian_model.identify_significant_bins(
            &bin_counts,
            self.bam_processor.total_reads(),
        )?;
        
        let merged_peaks = self.merge_bins_into_peaks(significant_bins, self.cli.mdist)?;
        
        let filtered_peaks = self.filter_peaks_by_width(merged_peaks, self.cli.minwidth)?;
        
        info!("Peak calling completed, found {} peaks", filtered_peaks.len());
        
        Ok(Peaks { ranges: filtered_peaks })
    }
    
    fn merge_bins_into_peaks(&self, bins: Vec<GenomicRange>, max_distance: u32) -> Result<Vec<GenomicRange>> {
        let mut chrom_bins: HashMap<String, Vec<GenomicRange>> = HashMap::new();
        
        for bin in bins {
            chrom_bins.entry(bin.chrom.clone()).or_default().push(bin);
        }
        
        let merged: Vec<GenomicRange> = chrom_bins
            .par_iter()
            .flat_map(|(chrom, bins)| {
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
                        current_peak.p_value *= bin.p_value;
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
    
    fn filter_peaks_by_width(&self, peaks: Vec<GenomicRange>, min_width: u32) -> Result<Vec<GenomicRange>> {
        let filtered = peaks
            .into_iter()
            .filter(|peak| peak.end - peak.start >= min_width)
            .collect();
        
        Ok(filtered)
    }
}

pub struct Peaks {
    ranges: Vec<GenomicRange>,
}

impl Peaks {
    pub fn write_to_bed(&self, path: &str) -> Result<usize> {
        use std::fs::File;
        use std::io::{BufWriter, Write};
        
        let file = File::create(path).context("Failed to create BED file")?;
        let mut writer = BufWriter::new(file);
        
        for range in &self.ranges {
            writeln!(
                writer,
                "{}\t{}\t{}\t.\t{:.6e}",
                range.chrom, range.start, range.end, range.p_value
            )?;
        }
        
        Ok(self.ranges.len())
    }
    
    pub fn len(&self) -> usize {
        self.ranges.len()
    }
    
    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }
}
