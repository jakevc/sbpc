use anyhow::{Context, Result};
use log::info;
use rust_htslib::bam::{self, Read};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::bam::GenomicRange;

pub struct Genome {
    pub seqnames: Vec<String>,
    pub lengths: Vec<u32>,
}

impl Genome {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        info!("Loading genome from file: {:?}", path.as_ref());
        
        let file = File::open(path).context("Failed to open chromosome sizes file")?;
        let reader = BufReader::new(file);
        
        let mut seqnames = Vec::new();
        let mut lengths = Vec::new();
        
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            
            if parts.len() >= 2 {
                let chrom = parts[0].to_string();
                let length = parts[1].parse::<u32>()
                    .context(format!("Failed to parse chromosome length: {}", parts[1]))?;
                
                seqnames.push(chrom);
                lengths.push(length);
            }
        }
        
        if seqnames.is_empty() {
            anyhow::bail!("No chromosomes found in the file");
        }
        
        Ok(Self { seqnames, lengths })
    }
    
    pub fn from_bam<P: AsRef<Path>>(path: P) -> Result<Self> {
        info!("Loading genome from BAM header: {:?}", path.as_ref());
        
        let bam = bam::Reader::from_path(path)
            .context("Failed to open BAM file")?;
        
        let header = bam.header();
        let target_names = header.target_names();
        
        let mut seqnames = Vec::new();
        let mut lengths = Vec::new();
        
        for (i, name) in target_names.iter().enumerate() {
            let chrom = std::str::from_utf8(name)
                .context("Failed to parse chromosome name")?
                .to_string();
            
            let length = header.target_len(i as u32)
                .context(format!("Failed to get length for chromosome {}", chrom))?
                as u32;
            
            seqnames.push(chrom);
            lengths.push(length);
        }
        
        if seqnames.is_empty() {
            anyhow::bail!("No chromosomes found in the BAM header");
        }
        
        let filtered = Self::filter_known_chromosomes(seqnames, lengths);
        
        Ok(filtered)
    }
    
    pub fn filter_known_chromosomes(seqnames: Vec<String>, lengths: Vec<u32>) -> Self {
        let filter_patterns = ["Un", "_", "EBV", "N", "M"];
        
        let mut filtered_seqnames = Vec::new();
        let mut filtered_lengths = Vec::new();
        
        for (i, chrom) in seqnames.iter().enumerate() {
            if !filter_patterns.iter().any(|&pattern| chrom.contains(pattern)) {
                filtered_seqnames.push(chrom.clone());
                filtered_lengths.push(lengths[i]);
            }
        }
        
        Self {
            seqnames: filtered_seqnames,
            lengths: filtered_lengths,
        }
    }
    
    pub fn create_bins(&self, step: u32, slide: u32) -> Result<Vec<GenomicRange>> {
        info!("Creating genomic bins with step={}, slide={}", step, slide);
        
        let mut bins = Vec::new();
        
        for (i, chrom) in self.seqnames.iter().enumerate() {
            let length = self.lengths[i];
            let mut start = 0;
            
            while start <= length - step {
                let end = start + step;
                
                bins.push(GenomicRange {
                    chrom: chrom.clone(),
                    start,
                    end,
                    p_value: 1.0, // Initialize p-value to 1.0
                });
                
                start += slide;
            }
        }
        
        info!("Created {} genomic bins", bins.len());
        
        Ok(bins)
    }
}
