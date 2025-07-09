#!/usr/bin/env python3

import argparse
import sys
import re

def filter_sam_by_tc(input_file, output_file, remove_tc_values):
    """
    Filter SAM file to remove entries with specified TC:i: values
    
    Args:
        input_file: Path to input SAM file
        output_file: Path to output SAM file  
        remove_tc_values: List of TC values to remove (e.g., [2, 3])
    """
    removed_count = 0
    kept_count = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Keep header lines
            if line.startswith('@'):
                outfile.write(line)
                continue
            
            # Check if line has TC:i: tag
            tc_match = re.search(r'TC:i:(\d+)', line)
            
            if tc_match:
                tc_value = int(tc_match.group(1))
                if tc_value in remove_tc_values:
                    # Remove this line
                    removed_count += 1
                    continue
            
            # Keep this line
            outfile.write(line)
            kept_count += 1
    
    return removed_count, kept_count

def analyze_sam_tc_distribution(input_file):
    """Analyze the distribution of TC:i: values in the SAM file"""
    tc_counts = {}
    total_records = 0
    
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('@'):
                continue
            
            total_records += 1
            tc_match = re.search(r'TC:i:(\d+)', line)
            
            if tc_match:
                tc_value = int(tc_match.group(1))
                tc_counts[tc_value] = tc_counts.get(tc_value, 0) + 1
            else:
                # Records without TC:i: (single records)
                tc_counts['none'] = tc_counts.get('none', 0) + 1
    
    return tc_counts, total_records

def main():
    parser = argparse.ArgumentParser(
        description='Filter SAM file to remove multi-record entries with specified TC:i: values',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Remove all records with TC:i:2
  python filter_sam_multirecords.py input.sam output.sam --remove-tc 2
  
  # Remove all records with TC:i:2 or TC:i:3
  python filter_sam_multirecords.py input.sam output.sam --remove-tc 2 3
  
  # Analyze TC distribution without filtering
  python filter_sam_multirecords.py input.sam --analyze-only
        """
    )
    
    parser.add_argument('input_sam', help='Input SAM file path')
    parser.add_argument('output_sam', nargs='?', help='Output SAM file path (required unless --analyze-only)')
    parser.add_argument('--remove-tc', type=int, nargs='+', 
                        help='TC:i: values to remove (e.g., --remove-tc 2 3)')
    parser.add_argument('--analyze-only', action='store_true',
                        help='Only analyze TC distribution, do not filter')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.analyze_only and not args.output_sam:
        parser.error('Output SAM file is required unless using --analyze-only')
    
    if not args.analyze_only and not args.remove_tc:
        parser.error('--remove-tc is required unless using --analyze-only')
    
    # Check if input file exists
    try:
        with open(args.input_sam, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_sam}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Analyze TC distribution
    print(f"Analyzing SAM file: {args.input_sam}")
    tc_counts, total_records = analyze_sam_tc_distribution(args.input_sam)
    
    print(f"\nTC:i: distribution:")
    print("-" * 30)
    for tc_value in sorted(tc_counts.keys(), key=lambda x: (x != 'none', x)):
        if tc_value == 'none':
            print(f"  No TC:i: (single records): {tc_counts[tc_value]:,}")
        else:
            print(f"  TC:i:{tc_value}: {tc_counts[tc_value]:,}")
    
    print(f"\nTotal records: {total_records:,}")
    
    # If analyze-only mode, exit here
    if args.analyze_only:
        return
    
    # Filter the file
    print(f"\nFiltering SAM file...")
    print(f"Removing records with TC:i: values: {args.remove_tc}")
    
    removed_count, kept_count = filter_sam_by_tc(args.input_sam, args.output_sam, args.remove_tc)
    
    print(f"\nFiltering complete:")
    print(f"  Records removed: {removed_count:,}")
    print(f"  Records kept: {kept_count:,}")
    print(f"  Output written to: {args.output_sam}")
    
    # Show what would remain
    remaining_tc_counts = {}
    for tc_value, count in tc_counts.items():
        if tc_value == 'none' or tc_value not in args.remove_tc:
            remaining_tc_counts[tc_value] = count
    
    print(f"\nRemaining TC:i: distribution:")
    print("-" * 30)
    for tc_value in sorted(remaining_tc_counts.keys(), key=lambda x: (x != 'none', x)):
        if tc_value == 'none':
            print(f"  No TC:i: (single records): {remaining_tc_counts[tc_value]:,}")
        else:
            print(f"  TC:i:{tc_value}: {remaining_tc_counts[tc_value]:,}")

if __name__ == "__main__":
    main()