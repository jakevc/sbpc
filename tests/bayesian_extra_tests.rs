use anyhow::Result;
use sbpc::bam::GenomicRange;
use sbpc::bayesian::BayesianModel;

#[test]
fn test_bayesian_model_min_reads_filter() -> Result<()> {
    let mut model = BayesianModel::new(0.05, 10);
    let bins = vec![
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 0,
                end: 100,
                posterior_prob: 1.0,
            },
            5,
        ), // below min_reads
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                posterior_prob: 1.0,
            },
            15,
        ), // above min_reads
    ];
    let total_reads = 100;
    let significant = model.identify_significant_bins(&bins, total_reads)?;
    assert!(significant.iter().all(|bin| bin.start != 0));
    Ok(())
}

#[test]
fn test_bayesian_model_posterior_prob_range() -> Result<()> {
    let mut model = BayesianModel::new(0.05, 1);
    let bins = vec![
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 0,
                end: 100,
                posterior_prob: 1.0,
            },
            50,
        ),
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                posterior_prob: 1.0,
            },
            25,
        ),
    ];
    let total_reads = 100;
    let significant = model.identify_significant_bins(&bins, total_reads)?;
    for bin in significant {
        assert!(bin.posterior_prob >= 0.0 && bin.posterior_prob <= 1.0);
    }
    Ok(())
}
