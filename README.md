# SBPC - Simple Bayesian Peak Caller

SBPC is a peak caller for genomic data that implements a Bayesian statistical framework to evaluate the probabilities of regions being peaks. It is designed for CUT&TAG/CUT&RUN sequencing data and provides p-values for each detected peak region.

## Usage

```
sbpc -b <sample>.bam -c <control>.bam > <sample>_peaks.bed
```

## Installation

### Pre-built Binaries

#### Linux

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/v0.1.3/sbpc-v0.1.3-x86_64-unknown-linux-gnu.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

#### macOS (Intel)

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/v0.1.3/sbpc-v0.1.3-x86_64-apple-darwin.tar.gz | tar xz
sudo mv sbpc /usr/local/bin/
```

#### macOS (Apple Silicon)

```bash
curl -L https://github.com/jakevc/sbpc/releases/download/v0.1.3/sbpc-v0.1.3-aarch64-apple-darwin.tar.gz | tar xz
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

## Requirements

SBPC requires input BAM files to be **coordinate-sorted** and **indexed**. You can prepare your BAM files using [samtools](http://www.htslib.org/doc/samtools.html):

```bash
samtools sort -o sorted.bam input.bam
samtools index sorted.bam
```

Both the sample and control BAM files must be sorted and have a corresponding `.bai` index file in the same directory.

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
  -t, --step <STEP>            Bin size for coverage bins [default: 200]
  -w, --minwidth <MINWIDTH>    Minimum width (bp) of a peak [default: 150]
  -v, --version                Print the current SBPC version
      --broad                  Run SBPC on broad marks (--step 5000)
      --verbose                Run SBPC in verbose mode
  -h, --help                   Print help
```

## Output

- Peaks are always printed in BED format to stdout.
- Metrics are printed as JSON to stdout by default. To also write metrics to a file (e.g., <prefix>_sbpc.json), use the `--metrics-file` flag.

## License

MIT
