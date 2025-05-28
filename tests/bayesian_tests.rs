use anyhow::Result;
use sbpc::bam::GenomicRange;
use sbpc::bayesian::BayesianModel;

#[test]
fn test_identify_significant_bins() -> Result<()> {
    let mut model = BayesianModel::new(0.05, 5);

    let bins = vec![
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 0,
                end: 100,
                posterior_prob: 1.0,
            },
            20, // High count, should be significant
        ),
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                posterior_prob: 1.0,
            },
            2, // Low count, should not be significant
        ),
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 200,
                end: 300,
                posterior_prob: 1.0,
            },
            15, // Medium count, might be significant
        ),
    ];

    let total_reads = 100;

    let significant_bins = model.identify_significant_bins(&bins, total_reads)?;

    assert!(
        !significant_bins.is_empty(),
        "Should have at least one significant bin"
    );

    let low_count_bins: Vec<_> = significant_bins
        .iter()
        .filter(|bin| bin.start == 100 && bin.end == 200)
        .collect();

    assert!(
        low_count_bins.is_empty(),
        "Low count bin should not be significant"
    );

    Ok(())
}

#[test]
fn test_bayesian_posterior_probability_calculation() -> Result<()> {
    let mut model = BayesianModel::new(0.05, 5);

    let bins = vec![
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 0,
                end: 100,
                posterior_prob: 1.0,
            },
            50, // Very high count
        ),
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                posterior_prob: 1.0,
            },
            25, // High count
        ),
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 200,
                end: 300,
                posterior_prob: 1.0,
            },
            10, // Medium count
        ),
        (
            GenomicRange {
                chrom: "chr1".to_string(),
                start: 300,
                end: 400,
                posterior_prob: 1.0,
            },
            5, // Low count
        ),
    ];

    let total_reads = 200;

    let significant_bins = model.identify_significant_bins(&bins, total_reads)?;

    for bin in &significant_bins {
        assert!(
            bin.posterior_prob >= 0.0 && bin.posterior_prob <= 1.0,
            "Posterior probability should be in range [0, 1], got: {}",
            bin.posterior_prob
        );
    }

    if significant_bins.len() >= 2 {
        let mut sorted_bins = significant_bins.clone();
        sorted_bins.sort_by(|a, b| a.start.cmp(&b.start));

        for i in 1..sorted_bins.len() {
            let prev_bin = &sorted_bins[i - 1];
            let curr_bin = &sorted_bins[i];

            let prev_count = bins
                .iter()
                .find(|(bin, _)| bin.start == prev_bin.start && bin.end == prev_bin.end)
                .map(|(_, count)| *count)
                .unwrap_or(0);

            let curr_count = bins
                .iter()
                .find(|(bin, _)| bin.start == curr_bin.start && bin.end == curr_bin.end)
                .map(|(_, count)| *count)
                .unwrap_or(0);

            if prev_count > curr_count {
                assert!(
                    prev_bin.posterior_prob >= curr_bin.posterior_prob,
                    "Bin with higher count should have higher posterior probability"
                );
            }
        }
    }

    Ok(())
}
