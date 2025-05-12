# SBPC - Simple Bayesian Peak Caller

SBPC is a peak caller for genomic data that implements a Bayesian statistical framework to evaluate the probabilities of regions being peaks. It is designed for CUT&TAG/CUT&RUN sequencing data and provides p-values for each detected peak region.

## Features

- Bayesian statistical framework for peak calling
- P-value calculation for each detected peak region
- Memory-efficient implementation that doesn't scale linearly with BAM file size
- Same input/output format as [gopeaks](https://github.com/maxsonBraunLab/gopeaks)
- Optimized for CUT&TAG/CUT&RUN data

## Installation

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
  -l, --slide <SLIDE>          Slide size for coverage bins [default: 50]
  -w, --minwidth <MINWIDTH>    Minimum width (bp) of a peak [default: 150]
  -o, --prefix <PREFIX>        Output prefix to write peaks and metrics file [default: sample]
  -v, --version                Print the current SBPC version
      --broad                  Run SBPC on broad marks (--step 5000 & --slide 1000)
      --verbose                Run SBPC in verbose mode
  -h, --help                   Print help
```

## Output

Two output files are generated each with the output prefix ${prefix}, set to "sample" by default.

- sample_peaks.bed - BED file containing peak locations
- sample_sbpc.json - JSON file containing metrics and p-values

## How It Works

SBPC uses a Bayesian statistical framework to evaluate the probabilities of regions being peaks. Unlike sliding window approaches, SBPC's Bayesian model allows for p-value multiplication of independent statistical tests, providing more accurate peak detection.

## License

MIT
