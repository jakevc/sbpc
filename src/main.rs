use clap::Parser;
use log::info;
use std::time::Instant;

mod bam;
mod bayesian;
mod cli;
mod genome;
mod metrics;
mod peak_caller;
mod utils;

fn main() -> anyhow::Result<()> {
    env_logger::init();

    let cli = cli::Cli::parse();

    if cli.version {
        println!("SBPC version: {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }

    let start_time = Instant::now();

    let mut peak_caller = peak_caller::PeakCaller::new(&cli)?;

    let peaks = peak_caller.call_peaks()?;

    let peak_count = peaks.write_to_bed(&format!("{}_peaks.bed", cli.prefix))?;

    let metrics = metrics::Metrics::new(
        env!("CARGO_PKG_VERSION"),
        &cli.prefix,
        &std::env::args().collect::<Vec<String>>().join(" "),
        peak_count,
        start_time.elapsed(),
    );

    metrics.write_to_file(&format!("{}_sbpc.json", cli.prefix))?;

    info!("Peak calling completed successfully");
    info!("Found {} peaks", peak_count);
    info!(
        "Results written to {}_peaks.bed and {}_sbpc.json",
        cli.prefix, cli.prefix
    );

    Ok(())
}
