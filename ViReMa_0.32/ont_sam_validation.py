#!/usr/bin/env python3

import sys
import os
import argparse
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Validate and compare SAM file statistics')
    parser.add_argument('sam_file', help='Primary SAM file to analyze (e.g., output_mm2_ont)')
    parser.add_argument('--output', '-o', help='Output file for results (default: stdout)')
    parser.add_argument('--comparison-files', '-c', nargs='*', 
                       default=['output.sam', 'multiRound'],
                       help='Additional files to compare against (default: output.sam multiRound)')
    return parser.parse_args()

def analyze_sam_file(file_path):
    """
    Analyze a SAM file and return statistics
    Returns dictionary with metrics
    """
    stats = {
        'file_path': file_path,
        'total_entries': 0,
        'header_lines': 0,
        'alignment_entries': 0,
        'unique_read_names': set(),
        'entries_with_n_gaps': 0,
        'mapped_entries': 0,
        'unmapped_entries': 0,
        'primary_alignments': 0,
        'supplementary_alignments': 0,
        'secondary_alignments': 0,
        'reads_with_n_gaps': set(),
        'cigar_operations': defaultdict(int),
        'file_exists': False,
        'error': None
    }
    
    # Check if file exists
    if not os.path.exists(file_path):
        stats['error'] = f"File not found: {file_path}"
        return stats
    
    stats['file_exists'] = True
    
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                stats['total_entries'] += 1
                
                if line.startswith('@'):
                    stats['header_lines'] += 1
                    continue
                
                # Parse SAM alignment line
                fields = line.split('\t')
                if len(fields) < 9:
                    continue  # Skip malformed lines
                
                stats['alignment_entries'] += 1
                
                read_name = fields[0]
                flag = int(fields[1]) if fields[1].isdigit() else 0
                cigar = fields[5] if len(fields) > 5 else '*'
                
                # Track unique read names
                stats['unique_read_names'].add(read_name)
                
                # Analyze flag
                if flag & 4:  # Unmapped
                    stats['unmapped_entries'] += 1
                else:
                    stats['mapped_entries'] += 1
                
                if flag & 256:  # Secondary alignment
                    stats['secondary_alignments'] += 1
                elif flag & 2048:  # Supplementary alignment
                    stats['supplementary_alignments'] += 1
                else:
                    stats['primary_alignments'] += 1
                
                # Analyze CIGAR string
                if cigar != '*':
                    if 'N' in cigar:
                        stats['entries_with_n_gaps'] += 1
                        stats['reads_with_n_gaps'].add(read_name)
                    
                    # Count CIGAR operations
                    import re
                    cigar_ops = re.findall(r'\d+([MIDNSHPX=])', cigar)
                    for op in cigar_ops:
                        stats['cigar_operations'][op] += 1
    
    except Exception as e:
        stats['error'] = f"Error reading file: {str(e)}"
    
    # Convert sets to counts for final output
    stats['unique_read_count'] = len(stats['unique_read_names'])
    stats['unique_reads_with_n_gaps'] = len(stats['reads_with_n_gaps'])
    
    return stats

def print_file_statistics(stats, output_file=None):
    """Print detailed statistics for a single file"""
    if output_file is None:
        output_file = sys.stdout
    
    file_name = os.path.basename(stats['file_path'])
    
    output_file.write(f"\n{'='*60}\n")
    output_file.write(f"STATISTICS FOR: {file_name}\n")
    output_file.write(f"{'='*60}\n")
    
    if stats['error']:
        output_file.write(f"ERROR: {stats['error']}\n")
        return
    
    output_file.write(f"File path: {stats['file_path']}\n")
    output_file.write(f"File exists: {stats['file_exists']}\n\n")
    
    output_file.write("ENTRY COUNTS:\n")
    output_file.write(f"  Total entries: {stats['total_entries']:,}\n")
    output_file.write(f"  Header lines: {stats['header_lines']:,}\n")
    output_file.write(f"  Alignment entries: {stats['alignment_entries']:,}\n\n")
    
    output_file.write("READ ANALYSIS:\n")
    output_file.write(f"  Unique read names: {stats['unique_read_count']:,}\n")
    output_file.write(f"  Entries with N gaps: {stats['entries_with_n_gaps']:,}\n")
    output_file.write(f"  Unique reads with N gaps: {stats['unique_reads_with_n_gaps']:,}\n")
    
    if stats['unique_read_count'] > 0:
        n_gap_percentage = (stats['unique_reads_with_n_gaps'] / stats['unique_read_count']) * 100
        output_file.write(f"  Percentage of reads with N gaps: {n_gap_percentage:.2f}%\n")
    
    output_file.write(f"\nMAPPING STATUS:\n")
    output_file.write(f"  Mapped entries: {stats['mapped_entries']:,}\n")
    output_file.write(f"  Unmapped entries: {stats['unmapped_entries']:,}\n")
    
    output_file.write(f"\nALIGNMENT TYPES:\n")
    output_file.write(f"  Primary alignments: {stats['primary_alignments']:,}\n")
    output_file.write(f"  Supplementary alignments: {stats['supplementary_alignments']:,}\n")
    output_file.write(f"  Secondary alignments: {stats['secondary_alignments']:,}\n")
    
    if stats['cigar_operations']:
        output_file.write(f"\nCIGAR OPERATIONS:\n")
        for op in sorted(stats['cigar_operations'].keys()):
            output_file.write(f"  {op}: {stats['cigar_operations'][op]:,}\n")

def print_comparison_table(all_stats, output_file=None):
    """Print comparison table across all files"""
    if output_file is None:
        output_file = sys.stdout
    
    output_file.write(f"\n{'='*80}\n")
    output_file.write(f"COMPARISON TABLE\n")
    output_file.write(f"{'='*80}\n")
    
    # Create table
    headers = ['Metric', 'File', 'Value']
    
    # Filter out files with errors
    valid_stats = [s for s in all_stats if not s['error']]
    
    if not valid_stats:
        output_file.write("No valid files to compare.\n")
        return
    
    # Key metrics to compare
    metrics = [
        ('total_entries', 'Total Entries'),
        ('alignment_entries', 'Alignment Entries'),
        ('unique_read_count', 'Unique Read Names'),
        ('entries_with_n_gaps', 'Entries with N Gaps'),
        ('unique_reads_with_n_gaps', 'Unique Reads with N Gaps'),
        ('mapped_entries', 'Mapped Entries'),
        ('unmapped_entries', 'Unmapped Entries'),
        ('primary_alignments', 'Primary Alignments'),
        ('supplementary_alignments', 'Supplementary Alignments')
    ]
    
    # Print table header
    output_file.write(f"{'Metric':<30} {'File':<20} {'Value':>15}\n")
    output_file.write(f"{'-'*30} {'-'*20} {'-'*15}\n")
    
    for metric_key, metric_name in metrics:
        for i, stats in enumerate(valid_stats):
            file_name = os.path.basename(stats['file_path'])
            value = stats.get(metric_key, 0)
            
            if i == 0:  # First file for this metric
                output_file.write(f"{metric_name:<30} {file_name:<20} {value:>15,}\n")
            else:
                output_file.write(f"{'':30} {file_name:<20} {value:>15,}\n")
        output_file.write("\n")

def print_summary_insights(all_stats, primary_file, output_file=None):
    """Print summary insights and potential issues"""
    if output_file is None:
        output_file = sys.stdout
    
    output_file.write(f"\n{'='*80}\n")
    output_file.write(f"SUMMARY INSIGHTS\n")
    output_file.write(f"{'='*80}\n")
    
    valid_stats = [s for s in all_stats if not s['error']]
    
    if len(valid_stats) < 2:
        output_file.write("Need at least 2 valid files for comparison insights.\n")
        return
    
    primary_stats = next((s for s in valid_stats if s['file_path'] == primary_file), None)
    if not primary_stats:
        output_file.write(f"Primary file {primary_file} not found in valid statistics.\n")
        return
    
    output_file.write(f"Primary file: {os.path.basename(primary_file)}\n\n")
    
    # Compare with other files
    for stats in valid_stats:
        if stats['file_path'] == primary_file:
            continue
        
        file_name = os.path.basename(stats['file_path'])
        output_file.write(f"Comparison with {file_name}:\n")
        
        # Unique reads comparison
        if stats['unique_read_count'] != primary_stats['unique_read_count']:
            diff = stats['unique_read_count'] - primary_stats['unique_read_count']
            output_file.write(f"  • Unique reads difference: {diff:+,}\n")
        
        # N-gap entries comparison
        if stats['entries_with_n_gaps'] != primary_stats['entries_with_n_gaps']:
            diff = stats['entries_with_n_gaps'] - primary_stats['entries_with_n_gaps']
            output_file.write(f"  • N-gap entries difference: {diff:+,}\n")
        
        # Unique reads with N-gaps comparison
        if stats['unique_reads_with_n_gaps'] != primary_stats['unique_reads_with_n_gaps']:
            diff = stats['unique_reads_with_n_gaps'] - primary_stats['unique_reads_with_n_gaps']
            output_file.write(f"  • Unique reads with N-gaps difference: {diff:+,}\n")
        
        output_file.write("\n")

def main():
    args = parse_arguments()
    
    # Collect all files to analyze
    all_files = [args.sam_file] + args.comparison_files
    all_stats = []
    
    print(f"Analyzing {len(all_files)} files...")
    
    # Analyze each file
    for file_path in all_files:
        print(f"Processing: {file_path}")
        stats = analyze_sam_file(file_path)
        all_stats.append(stats)
    
    # Open output file if specified
    output_file = open(args.output, 'w') if args.output else sys.stdout
    
    try:
        # Print individual file statistics
        for stats in all_stats:
            print_file_statistics(stats, output_file)
        
        # Print comparison table
        print_comparison_table(all_stats, output_file)
        
        # Print summary insights
        print_summary_insights(all_stats, args.sam_file, output_file)
        
    finally:
        if args.output:
            output_file.close()
            print(f"\nResults saved to: {args.output}")

if __name__ == "__main__":
    main()