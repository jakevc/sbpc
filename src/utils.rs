use anyhow::Result;
use log::debug;
use std::path::Path;

pub fn file_exists<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref().exists()
}

pub fn ensure_dir<P: AsRef<Path>>(path: P) -> Result<()> {
    let path = path.as_ref();
    if !path.exists() {
        std::fs::create_dir_all(path)?;
    }
    Ok(())
}

pub fn file_stem<P: AsRef<Path>>(path: P) -> Option<String> {
    path.as_ref()
        .file_stem()
        .and_then(|s| s.to_str())
        .map(String::from)
}

pub fn file_extension<P: AsRef<Path>>(path: P) -> Option<String> {
    path.as_ref()
        .extension()
        .and_then(|s| s.to_str())
        .map(String::from)
}

pub fn format_with_commas(num: usize) -> String {
    let mut s = String::new();
    let num_str = num.to_string();
    let len = num_str.len();
    
    for (i, c) in num_str.chars().enumerate() {
        s.push(c);
        if (len - i - 1) % 3 == 0 && i < len - 1 {
            s.push(',');
        }
    }
    
    s
}

pub fn calculate_memory_usage() -> Result<f64> {
    
    let mut memory_usage = 0.0;
    
    #[cfg(target_os = "linux")]
    {
        if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmRSS:") {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let Ok(kb) = parts[1].parse::<f64>() {
                            memory_usage = kb / 1024.0; // Convert KB to MB
                            break;
                        }
                    }
                }
            }
        }
    }
    
    debug!("Current memory usage: {:.2} MB", memory_usage);
    
    Ok(memory_usage)
}
