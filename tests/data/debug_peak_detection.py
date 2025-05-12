"""
Debug script to test SBPC peak detection with various parameters.
This script generates test data and runs SBPC with different parameter combinations.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.append(str(Path(__file__).parent))
from generate_test_data import main as generate_data

def run_sbpc_with_params(params):
    """Run SBPC with the given parameters and return the output."""
    cmd = [
        "../../target/release/sbpc",
        "-b", "test_sample.bam",
        "-c", "test_control.bam",
        "-o", "test_output"
    ] + params
    
    print(f"Running SBPC with parameters: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd, 
            cwd=str(Path(__file__).parent),
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error running SBPC: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return None, e.stderr

def count_peaks(output_file):
    """Count the number of peaks in the output file."""
    try:
        with open(output_file, 'r') as f:
            return len(f.readlines())
    except FileNotFoundError:
        return 0

def main():
    print("Building SBPC...")
    subprocess.run(
        ["cargo", "build", "--release"],
        cwd=str(Path(__file__).parent.parent.parent),
        check=True
    )
    
    print("Generating test data...")
    sys.argv = ["generate_test_data.py", "--reads", "10000", "--peaks", "5"]
    generate_data()
    
    parameter_sets = [
        ["--verbose", "-r", "1", "-p", "0.5", "-m", "500", "-w", "100"],
        ["--verbose", "-r", "0", "-p", "0.1"],
        ["--verbose", "-r", "5", "-p", "0.1", "-t", "50", "-l", "25"],
        ["--verbose", "-r", "5", "-p", "0.1", "--broad"],
        ["--verbose", "-r", "1", "-p", "0.9", "-m", "100", "-w", "50"],
    ]
    
    for i, params in enumerate(parameter_sets):
        print(f"\n=== Test {i+1}: {' '.join(params)} ===")
        stdout, stderr = run_sbpc_with_params(params)
        
        peak_file = Path(__file__).parent / "test_output_peaks.bed"
        metrics_file = Path(__file__).parent / "test_output_sbpc.json"
        
        num_peaks = count_peaks(peak_file)
        print(f"Detected {num_peaks} peaks")
        
        if num_peaks > 0:
            print("SUCCESS! Peaks detected.")
            print(f"Parameters that worked: {' '.join(params)}")
            
            print("\nDetected peaks:")
            with open(peak_file, 'r') as f:
                print(f.read())
                
            print("\nMetrics:")
            with open(metrics_file, 'r') as f:
                print(f.read())
                
            with open(Path(__file__).parent / "successful_params.txt", 'w') as f:
                f.write(f"./target/release/sbpc -b tests/data/test_sample.bam -c tests/data/test_control.bam -o tests/data/test_output {' '.join(params)}")
            
            return
    
    print("\nNo peaks detected with any parameter combination.")
    print("Consider modifying the test data generation to create more pronounced peaks.")

if __name__ == "__main__":
    main()
