"""
Generate a small synthetic BAM file for testing SBPC.
"""
import os
import pysam
import numpy as np
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Generate synthetic BAM file for testing')
    parser.add_argument('--output', default='test_sample.bam', help='Output BAM file')
    parser.add_argument('--control', default='test_control.bam', help='Output control BAM file')
    parser.add_argument('--reads', type=int, default=1000, help='Number of reads to generate')
    parser.add_argument('--peaks', type=int, default=5, help='Number of peaks to generate')
    return parser.parse_args()

def create_header():
    """Create a simple header with a few chromosomes."""
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'LN': 10000, 'SN': 'chr1'},
                     {'LN': 8000, 'SN': 'chr2'},
                     {'LN': 5000, 'SN': 'chr3'}],
              'PG': [{'ID': 'synthetic', 'PN': 'synthetic.py', 'VN': '1.0'}]}
    return header

def generate_peak_regions(num_peaks, header):
    """Generate random peak regions."""
    peaks = []
    for _ in range(num_peaks):
        chrom_idx = np.random.randint(0, len(header['SQ']))
        chrom = header['SQ'][chrom_idx]['SN']
        chrom_len = header['SQ'][chrom_idx]['LN']
        
        width = np.random.randint(150, 300)
        start = np.random.randint(0, chrom_len - width)
        end = start + width
        
        peaks.append((chrom, start, end))
    
    return peaks

def generate_reads(num_reads, peaks, header, is_control=False):
    """Generate random reads, with enrichment in peak regions for sample."""
    reads = []
    
    peak_prob = 0.9 if not is_control else 0.1
    
    if not is_control:
        for peak_idx, (chrom, peak_start, peak_end) in enumerate(peaks):
            for i in range(50):
                read_id = num_reads + peak_idx * 50 + i
                
                pos = np.random.randint(peak_start, peak_end - 50)
                tlen = np.random.randint(50, 100)  # Shorter fragments for more precise peaks
                
                read1 = pysam.AlignedSegment()
                read1.query_name = f"peak_read_{read_id}_1"
                read1.reference_id = next(idx for idx, sq in enumerate(header['SQ']) if sq['SN'] == chrom)
                read1.reference_start = pos
                read1.mapping_quality = 60
                read1.cigartuples = [(0, 75)]  # 75M (75bp match)
                read1.flag = 99  # Paired, mapped, first in pair
                read1.template_length = tlen
                read1.next_reference_id = read1.reference_id
                read1.next_reference_start = pos + tlen - 75
                
                read2 = pysam.AlignedSegment()
                read2.query_name = f"peak_read_{read_id}_1"
                read2.reference_id = read1.reference_id
                read2.reference_start = read1.next_reference_start
                read2.mapping_quality = 60
                read2.cigartuples = [(0, 75)]  # 75M (75bp match)
                read2.flag = 147  # Paired, mapped, second in pair, reverse strand
                read2.template_length = -tlen
                read2.next_reference_id = read1.reference_id
                read2.next_reference_start = read1.reference_start
                
                reads.append(read1)
                reads.append(read2)
    
    for i in range(num_reads):
        in_peak = np.random.random() < peak_prob
        
        if in_peak and peaks:
            peak_idx = np.random.randint(0, len(peaks))
            chrom, peak_start, peak_end = peaks[peak_idx]
            
            pos = np.random.randint(peak_start, peak_end)
            tlen = np.random.randint(50, 150)  # Fragment length
        else:
            chrom_idx = np.random.randint(0, len(header['SQ']))
            chrom = header['SQ'][chrom_idx]['SN']
            chrom_len = header['SQ'][chrom_idx]['LN']
            pos = np.random.randint(0, chrom_len)
            tlen = np.random.randint(50, 150)  # Fragment length
        
        read1 = pysam.AlignedSegment()
        read1.query_name = f"read_{i}_1"
        read1.reference_id = next(idx for idx, sq in enumerate(header['SQ']) if sq['SN'] == chrom)
        read1.reference_start = pos
        read1.mapping_quality = 60
        read1.cigartuples = [(0, 75)]  # 75M (75bp match)
        read1.flag = 99  # Paired, mapped, first in pair
        read1.template_length = tlen
        read1.next_reference_id = read1.reference_id
        read1.next_reference_start = pos + tlen - 75
        
        read2 = pysam.AlignedSegment()
        read2.query_name = f"read_{i}_1"
        read2.reference_id = read1.reference_id
        read2.reference_start = read1.next_reference_start
        read2.mapping_quality = 60
        read2.cigartuples = [(0, 75)]  # 75M (75bp match)
        read2.flag = 147  # Paired, mapped, second in pair, reverse strand
        read2.template_length = -tlen
        read2.next_reference_id = read1.reference_id
        read2.next_reference_start = read1.reference_start
        
        reads.append(read1)
        reads.append(read2)
    
    return reads

def write_bam(reads, header, filename):
    """Write reads to a BAM file."""
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    
    with pysam.AlignmentFile(filename, "wb", header=header) as outf:
        for read in reads:
            outf.write(read)
    
    pysam.index(filename)
    
    return filename

def main():
    args = parse_args()
    
    header = create_header()
    peaks = generate_peak_regions(args.peaks, header)
    
    sample_reads = generate_reads(args.reads, peaks, header, is_control=False)
    sample_bam = write_bam(sample_reads, header, args.output)
    
    control_reads = generate_reads(args.reads, peaks, header, is_control=True)
    control_bam = write_bam(control_reads, header, args.control)
    
    with open("test_true_peaks.bed", "w") as f:
        for i, (chrom, start, end) in enumerate(peaks):
            f.write(f"{chrom}\t{start}\t{end}\tpeak_{i}\t1000\t.\n")
    
    print(f"Generated {len(sample_reads)} reads in sample BAM: {args.output}")
    print(f"Generated {len(control_reads)} reads in control BAM: {args.control}")
    print(f"Generated {len(peaks)} true peaks in: test_true_peaks.bed")

if __name__ == "__main__":
    main()
