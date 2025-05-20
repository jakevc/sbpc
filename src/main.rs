use clap::Parser;
use log::info;
use std::time::Instant;

mod bam;
mod bayesian;
mod cli;
mod genome;
mod metrics;
mod peak_caller;

fn main() -> anyhow::Result<()> {
    // Set RUST_LOG if --verbose is passed, before env_logger::init()
    let cli = cli::Cli::parse();
    if cli.verbose {
        // Only set if not already set by user
        if std::env::var("RUST_LOG").is_err() {
            std::env::set_var("RUST_LOG", "info");
        }
    }
    env_logger::init();

    let start_time = Instant::now();

    let mut peak_caller = peak_caller::PeakCaller::new(&cli)?;

    let peaks = peak_caller.call_peaks()?;

    let peak_count = peaks.write_to_stdout_bed();

    let metrics = metrics::Metrics::new(
        env!("CARGO_PKG_VERSION"),
        &std::env::args().collect::<Vec<String>>().join(" "),
        peak_count,
        start_time.elapsed(),
    );

    info!("{}", serde_json::to_string_pretty(&metrics)?);
    // Write metrics to file if requested
    if cli.metrics_file {
        use std::path::Path;
        let bam_stem = Path::new(&cli.bam)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("metrics");
        let metrics_path = format!("{}_sbpc.json", bam_stem);
        metrics.write_to_file(&metrics_path)?;
    }

    info!("Peak calling completed successfully");
    info!("Found {} peaks", peak_count);
    info!(
        "Results printed to stdout; metrics {}",
        if cli.metrics_file { "also written to file" } else { "printed to stdout only" }
    );

    Ok(())
}
