#!/usr/bin/env python3
import sys
import os

def compare_sam_files(sam_file1, sam_file2, output_file="converted_comparison_final.txt"):
    """Compare reads between two SAM files"""
    
    # Read first SAM file
    small_reads = {}
    with open(sam_file1, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 1:
                read_name = fields[0]
                if read_name not in small_reads:
                    small_reads[read_name] = []
                small_reads[read_name].append(line.strip())
    
    # Read second SAM file
    conv_reads = {}
    with open(sam_file2, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 1:
                read_name = fields[0]
                if read_name not in conv_reads:
                    conv_reads[read_name] = []
                conv_reads[read_name].append(line.strip())
    
    # Get all unique read names from both files
    all_reads = set(small_reads.keys()) | set(conv_reads.keys())
    
    # Write comparison file
    with open(output_file, 'w') as f:
        f.write(f"# Comparison of {sam_file1} vs {sam_file2}\n")
        f.write("# Format: [SOURCE] followed by alignment line(s)\n")
        f.write(f"# FILE1 = {sam_file1}, FILE2 = {sam_file2}\n\n")
        
        for read_name in sorted(all_reads):
            f.write(f"=== READ: {read_name} ===\n")
            
            # Write first SAM file lines
            if read_name in small_reads:
                for line in small_reads[read_name]:
                    f.write(f"[FILE1] {line}\n")
            else:
                f.write(f"[FILE1] (No alignment found)\n")
            
            # Write second SAM file lines
            if read_name in conv_reads:
                for line in conv_reads[read_name]:
                    f.write(f"[FILE2] {line}\n")
            else:
                f.write(f"[FILE2] (No alignment found)\n")
            
            f.write("\n")
    
    print(f"Comparison complete. Found {len(all_reads)} unique reads.")
    print(f"{sam_file1}: {len(small_reads)} reads")
    print(f"{sam_file2}: {len(conv_reads)} reads")
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 compare_conversions.py <sam_file1> <sam_file2> [output_file]")
        print("Example: python3 compare_conversions.py file1.sam file2.sam comparison.txt")
        sys.exit(1)
    
    sam_file1 = sys.argv[1]
    sam_file2 = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else "converted_comparison_final.txt"
    
    # Check if input files exist
    if not os.path.exists(sam_file1):
        print(f"Error: {sam_file1} does not exist")
        sys.exit(1)
    if not os.path.exists(sam_file2):
        print(f"Error: {sam_file2} does not exist")
        sys.exit(1)
    
    compare_sam_files(sam_file1, sam_file2, output_file)