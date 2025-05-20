# SBPC - Simple Bayesian Peak Caller

SBPC is a peak caller for genomic data that implements a Bayesian statistical framework to evaluate the probabilities of regions being peaks. It is designed for CUT&TAG/CUT&RUN sequencing data and provides p-values for each detected peak region.

## Features

- Bayesian statistical framework for peak calling
- P-value calculation for each detected peak region
- Memory-efficient implementation that doesn't scale linearly with BAM file size
- Same input/output format as [gopeaks](https://github.com/maxsonBraunLab/gopeaks)
- Optimized for CUT&TAG/CUT&RUN data

## Installation

### Pre-built Binaries

#### Linux

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/v0.1.0/sbpc-v0.1.0-x86_64-unknown-linux-gnu.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

#### macOS (Intel)

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/v0.1.0/sbpc-v0.1.0-x86_64-apple-darwin.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

#### macOS (Apple Silicon)

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/v0.1.0/sbpc-v0.1.0-aarch64-apple-darwin.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

### Using Cargo

```bash
cargo install sbpc
```

### From Source

```bash
git clone https://github.com/jakevc/sbpc.git
cd sbpc
cargo build --release
```

The binary will be available at `target/release/sbpc`.

### Using Docker

You can also run SBPC using Docker, which avoids the need to install system dependencies:

```bash
# Pull the latest image
docker pull ghcr.io/jakevc/sbpc:latest

# Run SBPC with Docker
docker run --rm -v /path/to/data:/data ghcr.io/jakevc/sbpc:latest -b /data/sample.bam -c /data/control.bam -o /data/output_prefix
```

This approach eliminates issues with glibc dependencies across different platforms, making it easier to use SBPC regardless of your operating system.

## Usage

```
sbpc -b <sample>.bam -c <control>.bam -o path/to/output/<sample>
```

### Command Line Arguments

```
Usage: sbpc [OPTIONS] --bam <BAM>

Options:
  -b, --bam <BAM>              Input BAM file (must be paired-end reads)
  -c, --control <CONTROL>      Input BAM file with control signal to be normalized (e.g. IgG, Input)
  -s, --chromsize <CHROMSIZE>  Chromosome sizes for the genome if not found in the bam header
  -m, --mdist <MDIST>          Merge peaks within <mdist> base pairs [default: 1000]
  -r, --minreads <MINREADS>    Test genome bins with at least <minreads> read pairs [default: 15]
  -p, --pval <PVAL>            Define significance threshold <pval> with multiple hypothesis correction via Benjamini-Hochberg [default: 0.05]
  -t, --step <STEP>            Bin size for coverage bins [default: 100]
  -w, --minwidth <MINWIDTH>    Minimum width (bp) of a peak [default: 150]
  -o, --prefix <PREFIX>        Output prefix to write peaks and metrics file [default: sample]
  -v, --version                Print the current SBPC version
      --broad                  Run SBPC on broad marks (--step 5000)
      --verbose                Run SBPC in verbose mode
  -h, --help                   Print help
```

## Output

Two output files are generated each with the output prefix ${prefix}, set to "sample" by default.

- sample_peaks.bed - BED file containing peak locations
- sample_sbpc.json - JSON file containing metrics and p-values

## How It Works

SBPC uses a Bayesian statistical framework to evaluate the probabilities of regions being peaks. Unlike sliding window approaches, SBPC's Bayesian model allows for p-value multiplication of independent statistical tests, providing more accurate peak detection.

### Bayesian Statistical Framework

The algorithm implements a Bayesian approach to peak calling that:

1. **Prior Distribution**: Uses a Beta distribution as the prior for read enrichment, with parameters derived from the global read distribution. This captures the background expectation of read counts across the genome.

2. **Posterior Calculation**: For each genomic bin, calculates a posterior probability by updating the prior with observed read counts. This uses Bayes' theorem to determine the probability that a region is enriched beyond what would be expected by chance.

3. **Adaptive Parameter Estimation**: Automatically estimates the alpha and beta parameters of the Beta distribution based on the mean and variance of non-zero read counts, making the model adaptable to different datasets.

4. **Independent Tests**: Treats each genomic bin as an independent statistical test, allowing for the multiplication of p-values when merging adjacent significant regions.

### P-value Calculation

P-values are calculated through a multi-step process:

1. **Bayesian P-value**: For each bin, calculates a p-value representing the probability of observing the given read enrichment or greater under the background model.

2. **Posterior Adjustment**: Adjusts raw p-values based on the ratio of posterior to prior means, which helps account for local biases in read distribution.

3. **FDR Correction**: Applies Benjamini-Hochberg false discovery rate (FDR) correction to control for multiple testing, ensuring that the reported significant regions maintain the specified significance threshold.

4. **Peak Merging**: When adjacent significant bins are merged into peaks, their p-values are multiplied together, resulting in a combined p-value that represents the joint significance of the entire peak region.

5. **BED File Output**: Each peak's p-value is written to the score column (column 5) of the output BED file in scientific notation with 6 decimal places, allowing for downstream filtering and analysis.

### Memory Optimization

SBPC implements several strategies to optimize memory usage and avoid linear scaling with BAM file size:

1. **Efficient Binning**: Uses a smart binning strategy that only stores information for bins with reads, rather than pre-allocating bins for the entire genome.

2. **Streaming Processing**: Processes BAM files in a streaming fashion, avoiding loading the entire file into memory.

3. **Chromosome Filtering**: Automatically filters out non-standard chromosomes (e.g., unplaced contigs, mitochondrial DNA) to reduce memory footprint.

4. **Parallel Processing**: Utilizes Rayon for parallel processing of independent chromosomes, improving performance without increasing memory usage.

5. **Memory-Efficient Data Structures**: Uses compact data structures to represent genomic ranges and read counts, minimizing memory overhead.

These optimizations allow SBPC to process large BAM files with a relatively constant memory footprint, unlike approaches that scale linearly with input size.

## License

MIT
