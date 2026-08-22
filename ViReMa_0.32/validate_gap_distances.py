#!/usr/bin/env python3

import sys
import re
import argparse
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Validate gap distances in SAM file against reference FASTA')
    parser.add_argument('sam_file', help='Input SAM file (e.g., output_mm2.sam)')
    parser.add_argument('fasta_file', help='Reference FASTA file')
    parser.add_argument('--output', '-o', help='Output file for results (default: stdout)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    return parser.parse_args()

def load_fasta_reference(fasta_file):
    """Load reference sequences from FASTA file"""
    references = {}
    current_seq_name = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_seq_name:
                    references[current_seq_name] = ''.join(current_seq)
                # Start new sequence
                current_seq_name = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)
    
    # Save last sequence
    if current_seq_name:
        references[current_seq_name] = ''.join(current_seq)
    
    return references

def parse_cigar_for_segments(cigar):
    """
    Parse CIGAR string to identify sequence segments separated by N gaps
    Returns list of (operation_type, length) tuples
    """
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    return [(op, int(length)) for length, op in cigar_ops]

def extract_sequence_segments(sequence, cigar_ops):
    """
    Extract sequence segments from read sequence based on CIGAR operations
    Returns list of sequence segments that correspond to M/=/X operations
    """
    segments = []
    seq_pos = 0
    
    for op, length in cigar_ops:
        if op in 'MIS=X':  # Operations that consume query sequence
            if op in 'M=X':  # Operations that represent alignment to reference
                segments.append(sequence[seq_pos:seq_pos + length])
            seq_pos += length
        # N, D, H, P don't consume query sequence
    
    return segments

def calculate_reference_positions(start_pos, cigar_ops):
    """
    Calculate reference positions for each aligned segment
    Returns list of (start_pos, end_pos) for each segment
    """
    positions = []
    ref_pos = start_pos - 1  # Convert to 0-based indexing
    
    for op, length in cigar_ops:
        if op in 'M=X':  # Alignment operations
            positions.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in 'DN':  # Operations that consume reference but not query
            ref_pos += length
        # I, S, H, P don't consume reference sequence
    
    return positions

def validate_read_against_reference(read_fields, references, verbose=False):
    """
    Validate a single read against the reference
    Returns validation results as a dictionary
    """
    read_name = read_fields[0]
    flag = int(read_fields[1])
    ref_name = read_fields[2]
    start_pos = int(read_fields[3])
    cigar = read_fields[5]
    sequence = read_fields[9]
    
    # Skip unmapped reads
    if flag & 4 or cigar == '*':
        return None
    
    # Only process reads with N gaps
    if 'N' not in cigar:
        return None
    
    # Check if reference exists
    if ref_name not in references:
        return {
            'read_name': read_name,
            'error': f'Reference {ref_name} not found in FASTA file'
        }
    
    reference_seq = references[ref_name]
    cigar_ops = parse_cigar_for_segments(cigar)
    
    # Extract sequence segments from read
    read_segments = extract_sequence_segments(sequence, cigar_ops)
    
    # Calculate reference positions for each segment
    ref_positions = calculate_reference_positions(start_pos, cigar_ops)
    
    # Validate each segment
    results = {
        'read_name': read_name,
        'ref_name': ref_name,
        'start_pos': start_pos,
        'cigar': cigar,
        'segments': [],
        'all_match': True
    }
    
    segment_idx = 0
    for i, (op, length) in enumerate(cigar_ops):
        if op in 'M=X':  # Aligned segments
            if segment_idx < len(read_segments) and segment_idx < len(ref_positions):
                read_seg = read_segments[segment_idx]
                ref_start, ref_end = ref_positions[segment_idx]
                
                # Extract reference sequence for this segment
                if ref_end <= len(reference_seq):
                    ref_seg = reference_seq[ref_start:ref_end]
                    match = (read_seg.upper() == ref_seg.upper())
                    
                    segment_result = {
                        'segment_num': segment_idx + 1,
                        'read_sequence': read_seg,
                        'ref_sequence': ref_seg,
                        'ref_start': ref_start + 1,  # Convert back to 1-based
                        'ref_end': ref_end,
                        'match': match
                    }
                    
                    if not match:
                        results['all_match'] = False
                        if verbose:
                            segment_result['mismatch_details'] = compare_sequences(read_seg, ref_seg)
                    
                    results['segments'].append(segment_result)
                else:
                    results['segments'].append({
                        'segment_num': segment_idx + 1,
                        'error': f'Reference position {ref_end} exceeds reference length {len(reference_seq)}'
                    })
                    results['all_match'] = False
                
                segment_idx += 1
    
    return results

def compare_sequences(seq1, seq2):
    """Compare two sequences and return mismatch details"""
    mismatches = []
    for i, (c1, c2) in enumerate(zip(seq1.upper(), seq2.upper())):
        if c1 != c2:
            mismatches.append({
                'position': i + 1,
                'read_base': c1,
                'ref_base': c2
            })
    return mismatches

def main():
    args = parse_arguments()
    
    # Load reference sequences
    print(f"Loading reference sequences from {args.fasta_file}...")
    references = load_fasta_reference(args.fasta_file)
    print(f"Loaded {len(references)} reference sequences")
    
    # Process SAM file
    print(f"Processing SAM file {args.sam_file}...")
    
    # Open output file if specified
    output_file = open(args.output, 'w') if args.output else sys.stdout
    
    try:
        total_reads = 0
        reads_with_gaps = 0
        matching_reads = 0
        mismatching_reads = 0
        
        with open(args.sam_file, 'r') as sam_file:
            for line in sam_file:
                if line.startswith('@'):
                    continue  # Skip header lines
                
                total_reads += 1
                fields = line.strip().split('\t')
                
                if len(fields) < 10:
                    continue
                
                result = validate_read_against_reference(fields, references, args.verbose)
                
                if result is None:
                    continue  # Skip reads without N gaps or unmapped reads
                
                reads_with_gaps += 1
                
                if 'error' in result:
                    output_file.write(f"ERROR in {result['read_name']}: {result['error']}\n")
                    continue
                
                # Output results
                if result['all_match']:
                    matching_reads += 1
                    if args.verbose:
                        output_file.write(f"MATCH: {result['read_name']} - All segments match reference\n")
                        for seg in result['segments']:
                            output_file.write(f"  Segment {seg['segment_num']}: {seg['ref_start']}-{seg['ref_end']} = {seg['read_sequence']}\n")
                else:
                    mismatching_reads += 1
                    output_file.write(f"MISMATCH: {result['read_name']} ({result['ref_name']}:{result['start_pos']}) CIGAR: {result['cigar']}\n")
                    for seg in result['segments']:
                        if 'error' in seg:
                            output_file.write(f"  Segment {seg['segment_num']}: ERROR - {seg['error']}\n")
                        else:
                            status = "MATCH" if seg['match'] else "MISMATCH"
                            output_file.write(f"  Segment {seg['segment_num']} ({status}): Ref {seg['ref_start']}-{seg['ref_end']}\n")
                            output_file.write(f"    Read: {seg['read_sequence']}\n")
                            output_file.write(f"    Ref:  {seg['ref_sequence']}\n")
                            
                            if not seg['match'] and args.verbose and 'mismatch_details' in seg:
                                output_file.write(f"    Mismatches: ")
                                for mm in seg['mismatch_details']:
                                    output_file.write(f"pos{mm['position']}({mm['read_base']}->{mm['ref_base']}) ")
                                output_file.write("\n")
                    output_file.write("\n")
        
        # Summary
        summary = f"""
VALIDATION SUMMARY:
Total reads processed: {total_reads:,}
Reads with N gaps: {reads_with_gaps:,}
Matching reads: {matching_reads:,}
Mismatching reads: {mismatching_reads:,}
Match rate: {(matching_reads / reads_with_gaps * 100) if reads_with_gaps > 0 else 0:.2f}%
"""
        output_file.write(summary)
        print(summary)
        
    finally:
        if args.output:
            output_file.close()

if __name__ == "__main__":
    main()