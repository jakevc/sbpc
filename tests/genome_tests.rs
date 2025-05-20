use anyhow::Result;
use sbpc::genome::Genome;
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn test_create_bins() -> Result<()> {
    let seqnames = vec!["chr1".to_string(), "chr2".to_string()];
    let lengths = vec![1000, 500];
    let genome = Genome { seqnames, lengths };

    let step = 100;
    let slide = 50;

    // With non-overlapping bins, slide is ignored and only step is used
    let bins = genome.create_bins(step, slide)?;

    // chr1: 1000/100 = 10 bins, chr2: 500/100 = 5 bins, total = 15
    assert_eq!(bins.len(), 15, "Expected 15 bins for non-overlapping bins");

    assert_eq!(bins[0].chrom, "chr1");
    assert_eq!(bins[0].start, 0);
    assert_eq!(bins[0].end, 100);

    assert_eq!(bins[9].chrom, "chr1");
    assert_eq!(bins[9].start, 900);
    assert_eq!(bins[9].end, 1000);

    assert_eq!(bins[10].chrom, "chr2");
    assert_eq!(bins[10].start, 0);
    assert_eq!(bins[10].end, 100);

    Ok(())
}

#[test]
fn test_from_file() -> Result<()> {
    let dir = tempdir()?;
    let file_path = dir.path().join("test_chromsizes.txt");

    let mut file = File::create(&file_path)?;
    writeln!(file, "chr1\t1000")?;
    writeln!(file, "chr2\t500")?;
    writeln!(file, "chr3\t2000")?;

    let genome = Genome::from_file(&file_path)?;

    assert_eq!(genome.seqnames.len(), 3);
    assert_eq!(genome.lengths.len(), 3);

    assert_eq!(genome.seqnames[0], "chr1");
    assert_eq!(genome.lengths[0], 1000);

    assert_eq!(genome.seqnames[1], "chr2");
    assert_eq!(genome.lengths[1], 500);

    assert_eq!(genome.seqnames[2], "chr3");
    assert_eq!(genome.lengths[2], 2000);

    Ok(())
}

#[test]
fn test_filter_known_chromosomes() {
    let seqnames = vec![
        "chr1".to_string(),
        "chr2".to_string(),
        "chrUn_random".to_string(),
        "chr_EBV".to_string(),
        "chrM".to_string(),
    ];
    let lengths = vec![1000, 500, 300, 200, 100];

    let filtered = Genome::filter_known_chromosomes(seqnames, lengths);

    assert_eq!(filtered.seqnames.len(), 2);
    assert_eq!(filtered.lengths.len(), 2);

    assert_eq!(filtered.seqnames[0], "chr1");
    assert_eq!(filtered.lengths[0], 1000);

    assert_eq!(filtered.seqnames[1], "chr2");
    assert_eq!(filtered.lengths[1], 500);
}
