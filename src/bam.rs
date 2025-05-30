use anyhow::{Context, Result};
use log::info;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct GenomicRange {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub posterior_prob: f64,
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
        let mut bam =
            bam::Reader::from_path(path).context(format!("Failed to open BAM file: {}", path))?;

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
        use rayon::prelude::*;
        use rust_htslib::bam::IndexedReader;
        info!(
            "Counting reads in {} bins using fast arithmetic bin assignment",
            bins.len()
        );

        // Group bins by chromosome and build bin lookup tables
        let mut chrom_bins: HashMap<String, Vec<(usize, &GenomicRange)>> = HashMap::new();
        for (i, bin) in bins.iter().enumerate() {
            chrom_bins
                .entry(bin.chrom.clone())
                .or_default()
                .push((i, bin));
        }

        // Parallelize by chromosome, collect local results
        let per_chrom_results: Vec<Vec<(usize, usize)>> = chrom_bins
            .par_iter()
            .map(|(chrom, bin_list)| {
                let mut bam =
                    IndexedReader::from_path(&self.bam_path).expect("Failed to open BAM file");
                let header = bam.header().to_owned();
                let tid = header
                    .target_names()
                    .iter()
                    .position(|name| std::str::from_utf8(name).unwrap() == chrom);
                let tid = match tid {
                    Some(tid) => tid as u32,
                    None => return Vec::new(),
                };
                let min_start = bin_list[0].1.start;
                let bin_size = bin_list[0].1.end - bin_list[0].1.start;
                let num_bins = bin_list.len();
                let max_end = bin_list[num_bins - 1].1.end;
                if bam.fetch((tid, min_start, max_end)).is_err() {
                    return Vec::new();
                }
                let mut local_counts = vec![0usize; num_bins];
                for rec in bam.records() {
                    let rec = match rec {
                        Ok(r) => r,
                        Err(_) => continue,
                    };
                    if rec.is_unmapped() {
                        continue;
                    }
                    let pos = rec.pos() as u32;
                    if pos < min_start {
                        continue;
                    }
                    let bin_idx = ((pos - min_start) / bin_size) as usize;
                    if bin_idx < num_bins {
                        local_counts[bin_idx] += 1;
                    }
                }
                // Return (global_idx, count) for each bin in this chromosome
                bin_list
                    .iter()
                    .enumerate()
                    .map(|(i, (global_idx, _))| (*global_idx, local_counts[i]))
                    .collect()
            })
            .collect();

        // Merge local results into global bin_counts
        let mut bin_counts = vec![0; bins.len()];
        for chrom_result in per_chrom_results {
            for (global_idx, count) in chrom_result {
                bin_counts[global_idx] = count;
            }
        }

        // If control_path is set, normalize to control
        if let (Some(control_path), Some(_)) = (&self.control_path, &self.control_reads) {
            self.normalize_to_control(bins, &mut bin_counts, control_path)?;
        }

        let result = bins
            .iter()
            .cloned()
            .zip(bin_counts.iter().cloned())
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

        let target_names: Vec<Vec<u8>> = control_bam
            .header()
            .target_names()
            .iter()
            .map(|&name| name.to_vec())
            .collect();

        for record in control_bam.records() {
            let rec = record?;

            if rec.is_unmapped() {
                continue;
            }

            let tid = rec.tid();
            let chrom = if tid >= 0 && (tid as usize) < target_names.len() {
                std::str::from_utf8(&target_names[tid as usize])?.to_string()
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
