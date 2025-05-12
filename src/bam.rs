use anyhow::{Context, Result};
use log::{debug, info};
use rust_htslib::bam::{self, Read, Record};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct GenomicRange {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub p_value: f64,
}

pub struct BamProcessor {
    bam_path: String,
    control_path: Option<String>,
    total_reads: usize,
    control_reads: Option<usize>,
}

impl BamProcessor {
    pub fn new(bam_path: &str, control_path: Option<&str>) -> Result<Self> {
        let total_reads = Self::count_total_reads(bam_path)?;
        
        let control_reads = if let Some(control) = control_path {
            Some(Self::count_total_reads(control)?)
        } else {
            None
        };
        
        Ok(Self {
            bam_path: bam_path.to_string(),
            control_path: control_path.map(String::from),
            total_reads,
            control_reads,
        })
    }
    
    fn count_total_reads(path: &str) -> Result<usize> {
        let mut bam = bam::Reader::from_path(path)
            .context(format!("Failed to open BAM file: {}", path))?;
        
        let mut count = 0;
        for record in bam.records() {
            let _rec = record?;
            count += 1;
        }
        
        Ok(count)
    }
    
    pub fn total_reads(&self) -> usize {
        self.total_reads
    }
    
    pub fn count_reads_in_bins(&self, bins: &[GenomicRange]) -> Result<Vec<(GenomicRange, usize)>> {
        info!("Counting reads in {} bins", bins.len());
        
        let mut chrom_bins: HashMap<String, Vec<(usize, &GenomicRange)>> = HashMap::new();
        for (i, bin) in bins.iter().enumerate() {
            chrom_bins
                .entry(bin.chrom.clone())
                .or_default()
                .push((i, bin));
        }
        
        let mut bin_counts = vec![0; bins.len()];
        
        let mut bam = bam::Reader::from_path(&self.bam_path)
            .context(format!("Failed to open BAM file: {}", self.bam_path))?;
        
        let target_names = bam.header().target_names().to_owned();
        
        for record in bam.records() {
            let rec = record?;
            
            if rec.is_unmapped() {
                continue;
            }
            
            let tid = rec.tid();
            let chrom = if tid >= 0 && (tid as usize) < target_names.len() {
                std::str::from_utf8(target_names[tid as usize])?.to_string()
            } else {
                continue;
            };
            
            let pos = rec.pos() as u32;
            
            if let Some(chrom_bin_list) = chrom_bins.get(&chrom) {
                for (bin_idx, bin) in chrom_bin_list {
                    if pos >= bin.start && pos < bin.end {
                        bin_counts[*bin_idx] += 1;
                    }
                }
            }
        }
        
        if let Some(control_path) = &self.control_path {
            self.normalize_to_control(bins, &mut bin_counts, control_path)?;
        }
        
        let result = bins
            .iter()
            .zip(bin_counts.iter())
            .map(|(bin, &count)| {
                let mut bin_clone = bin.clone();
                bin_clone.p_value = 1.0; // Initialize p-value to 1.0, will be calculated by Bayesian model
                (bin_clone, count)
            })
            .collect();
        
        Ok(result)
    }
    
    fn normalize_to_control(
        &self,
        bins: &[GenomicRange],
        bin_counts: &mut [usize],
        control_path: &str,
    ) -> Result<()> {
        info!("Normalizing counts using control BAM: {}", control_path);
        
        let mut chrom_bins: HashMap<String, Vec<(usize, &GenomicRange)>> = HashMap::new();
        for (i, bin) in bins.iter().enumerate() {
            chrom_bins
                .entry(bin.chrom.clone())
                .or_default()
                .push((i, bin));
        }
        
        let mut control_counts = vec![0; bins.len()];
        
        let mut control_bam = bam::Reader::from_path(control_path)
            .context(format!("Failed to open control BAM file: {}", control_path))?;
        
        let target_names = control_bam.header().target_names().to_owned();
        
        for record in control_bam.records() {
            let rec = record?;
            
            if rec.is_unmapped() {
                continue;
            }
            
            let tid = rec.tid();
            let chrom = if tid >= 0 && (tid as usize) < target_names.len() {
                std::str::from_utf8(target_names[tid as usize])?.to_string()
            } else {
                continue;
            };
            
            let pos = rec.pos() as u32;
            
            if let Some(chrom_bin_list) = chrom_bins.get(&chrom) {
                for (bin_idx, bin) in chrom_bin_list {
                    if pos >= bin.start && pos < bin.end {
                        control_counts[*bin_idx] += 1;
                    }
                }
            }
        }
        
        let scale_factor = self.total_reads as f64 / self.control_reads.unwrap_or(1) as f64;
        
        for i in 0..bin_counts.len() {
            let treatment_count = bin_counts[i] as f64;
            let control_count = control_counts[i] as f64 * scale_factor;
            
            let normalized_count = (treatment_count - control_count).max(0.0);
            bin_counts[i] = normalized_count.round() as usize;
        }
        
        Ok(())
    }
}
