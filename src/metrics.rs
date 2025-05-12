use anyhow::Result;
use chrono::Local;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::time::Duration;

#[derive(Serialize, Deserialize, Debug)]
pub struct Metrics {
    #[serde(rename = "sbpc_version")]
    pub version: String,
    pub date: String,
    pub elapsed: String,
    pub prefix: String,
    pub command: String,
    #[serde(rename = "peak_counts")]
    pub peaks: usize,
}

impl Metrics {
    pub fn new(
        version: &str,
        prefix: &str,
        command: &str,
        peaks: usize,
        elapsed: Duration,
    ) -> Self {
        Self {
            version: version.to_string(),
            date: Local::now().format("%Y-%m-%d %I:%M:%S %p").to_string(),
            elapsed: format!("{:?}", elapsed),
            prefix: prefix.to_string(),
            command: command.to_string(),
            peaks,
        }
    }

    pub fn write_to_file(&self, path: &str) -> Result<()> {
        let json = serde_json::to_string_pretty(self)?;
        let mut file = File::create(path)?;
        file.write_all(json.as_bytes())?;
        file.write_all(b"\n")?;
        Ok(())
    }
}
