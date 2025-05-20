use anyhow::Result;
use sbpc::metrics::Metrics;
use std::fs;
use std::time::Duration;
use tempfile::tempdir;

#[test]
fn test_metrics_creation() {
    let version = "0.1.3";
    let prefix = "test_sample";
    let command = "sbpc -b test.bam -o test_sample";
    let peaks = 1000;
    let elapsed = Duration::from_secs(60);

    let metrics = Metrics::new(version, prefix, command, peaks, elapsed);

    assert_eq!(metrics.version, version);
    assert_eq!(metrics.prefix, prefix);
    assert_eq!(metrics.command, command);
    assert_eq!(metrics.peaks, peaks);

    assert!(metrics.elapsed.contains("60s"));
}

#[test]
fn test_write_to_file() -> Result<()> {
    let dir = tempdir()?;
    let file_path = dir.path().join("test_metrics.json");

    let version = "0.1.3";
    let prefix = "test_sample";
    let command = "sbpc -b test.bam -o test_sample";
    let peaks = 1000;
    let elapsed = Duration::from_secs(60);

    let metrics = Metrics::new(version, prefix, command, peaks, elapsed);

    metrics.write_to_file(file_path.to_str().unwrap())?;

    assert!(file_path.exists());

    let contents = fs::read_to_string(file_path)?;

    assert!(contents.contains("\"sbpc_version\""));
    assert!(contents.contains("\"date\""));
    assert!(contents.contains("\"elapsed\""));
    assert!(contents.contains("\"prefix\""));
    assert!(contents.contains("\"command\""));
    assert!(contents.contains("\"peak_counts\""));

    assert!(contents.contains(&format!("\"sbpc_version\": \"{}\"", version)));
    assert!(contents.contains(&format!("\"prefix\": \"{}\"", prefix)));
    assert!(contents.contains(&format!("\"command\": \"{}\"", command)));
    assert!(contents.contains(&format!("\"peak_counts\": {}", peaks)));

    Ok(())
}
