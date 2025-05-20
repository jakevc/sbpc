use anyhow::Result;
use sbpc::bam::{BamProcessor, GenomicRange};

#[test]
fn test_bam_processor_total_reads() -> Result<()> {
    // Create a small BAM file for testing (simulate with a dummy file for now)
    // In a real test, use a real BAM file or a fixture
    let bam_path = "tests/data/test_sample.bam";
    let processor = BamProcessor::new(bam_path, None)?;
    let total = processor.total_reads();
    assert!(total > 0, "Total reads should be positive");
    Ok(())
}

#[test]
fn test_count_reads_in_bins_empty() -> Result<()> {
    let bam_path = "tests/data/test_sample.bam";
    let processor = BamProcessor::new(bam_path, None)?;
    let bins: Vec<GenomicRange> = vec![];
    let counts = processor.count_reads_in_bins(&bins)?;
    assert!(counts.is_empty(), "Counts should be empty for empty bins");
    Ok(())
}

#[test]
fn test_count_reads_in_bins_basic() -> Result<()> {
    let bam_path = "tests/data/test_sample.bam";
    let processor = BamProcessor::new(bam_path, None)?;
    let bins = vec![
        GenomicRange {
            chrom: "chr1".to_string(),
            start: 0,
            end: 100,
            p_value: 1.0,
        },
        GenomicRange {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            p_value: 1.0,
        },
    ];
    let counts = processor.count_reads_in_bins(&bins)?;
    assert_eq!(counts.len(), 2);
    for (bin, _count) in counts {
        assert_eq!(bin.chrom, "chr1");
        assert!(bin.start < bin.end);
    }
    Ok(())
}
