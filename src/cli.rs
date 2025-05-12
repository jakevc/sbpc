use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[arg(short = 'b', long)]
    pub bam: String,

    #[arg(short = 'c', long)]
    pub control: Option<String>,

    #[arg(short = 's', long)]
    pub chromsize: Option<String>,

    #[arg(short = 'm', long, default_value_t = 1000)]
    pub mdist: u32,

    #[arg(short = 'r', long, default_value_t = 15)]
    pub minreads: u32,

    #[arg(short = 'p', long, default_value_t = 0.05)]
    pub pval: f64,

    #[arg(short = 't', long, default_value_t = 100)]
    pub step: u32,

    #[arg(short = 'l', long, default_value_t = 50)]
    pub slide: u32,

    #[arg(short = 'w', long, default_value_t = 150)]
    pub minwidth: u32,

    #[arg(short = 'o', long, default_value = "sample")]
    pub prefix: String,

    #[arg(short = 'v', long)]
    pub version: bool,

    #[arg(long)]
    pub broad: bool,

    #[arg(long)]
    pub verbose: bool,
}

impl Cli {
    pub fn adjust_for_broad_peaks(&mut self) {
        if self.broad {
            self.step = 5000;
            self.slide = 1000;
        }
    }
}
