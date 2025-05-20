use anyhow::Result;
use sbpc::genome::Genome;

#[test]
fn test_create_bins_non_overlapping() -> Result<()> {
    let seqnames = vec!["chr1".to_string()];
    let lengths = vec![1000];
    let genome = Genome { seqnames, lengths };
    let bins = genome.create_bins(100, 100)?;
    assert_eq!(bins.len(), 10);
    for (i, bin) in bins.iter().enumerate() {
        assert_eq!(bin.start, i as u32 * 100);
        assert_eq!(bin.end, (i as u32 + 1) * 100);
    }
    Ok(())
}

#[test]
fn test_filter_known_chromosomes_excludes() {
    let seqnames = vec!["chr1".to_string(), "chrUn_random".to_string(), "chrM".to_string()];
    let lengths = vec![1000, 500, 100];
    let filtered = Genome::filter_known_chromosomes(seqnames, lengths);
    assert_eq!(filtered.seqnames, vec!["chr1"]);
    assert_eq!(filtered.lengths, vec![1000]);
}
