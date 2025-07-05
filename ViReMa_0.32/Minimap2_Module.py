#!/usr/bin/env python3

import subprocess
import re
import sys
from collections import defaultdict

# Global variables for file paths and parameters
VIRUS_INDEX = None
INPUT_DATA = None
OUTPUT_SAM = None
SEED_THRESHOLD = 25  # Default threshold

def run_minimap2_workflow(config):
    """
    Main entry point for minimap2 workflow when called from ViReMa
    Uses ViReMa configuration parameters
    """
    global VIRUS_INDEX, INPUT_DATA, OUTPUT_SAM, SEED_THRESHOLD
    
    # Set parameters from ViReMa config
    VIRUS_INDEX = config.Lib1
    INPUT_DATA = config.File1
    OUTPUT_SAM = config.Output_Dir + config.File3
    SEED_THRESHOLD = int(config.Seed) if config.Seed else 25
    
    # Determine long read technology (empty string means short read mode)
    long_read_tech = config.LongReadTech if config.LongReadTech else None
    
    print(f"Minimap2 parameters:")
    print(f"  Virus Index: {VIRUS_INDEX}")
    print(f"  Input Data: {INPUT_DATA}")
    print(f"  Output SAM: {OUTPUT_SAM}")
    print(f"  Seed Threshold: {SEED_THRESHOLD}")
    print(f"  Long Read Tech: {long_read_tech if long_read_tech else 'short read'}")
    
    # Run the minimap2 workflow
    try:
        # Step 1: Run initial minimap2 alignment
        cmd = build_minimap2_command(INPUT_DATA, long_read_tech, is_initial=True)
        with open('output.sam', 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        print("Initial minimap2 alignment completed")
        
        # Step 2: Parse SAM for softclipped reads with only primary alignment
        softclipped_reads = parse_sam_for_softclipped(long_read_tech)
        
        if softclipped_reads:
            # Step 3: Save to multiRound file (debugging)
            with open('multiRound', 'w') as f:
                for read_name, data in softclipped_reads.items():
                    f.write(data['line'] + '\n')
            print(f"Saved {len(softclipped_reads)} softclipped reads to multiRound")
            
            # Step 4: Extract softclipped sequences
            temp_sequences = []
            for read_name, data in softclipped_reads.items():
                # Use sequence directly from SAM file (stored in data['sequence'])
                # This ensures we use the exact sequence that corresponds to the CIGAR
                original_seq = data['sequence']
                cigar = data['cigar']

                # Parse CIGAR to find softclipped regions with their positions relative to main alignment
                softclipped_seqs = extract_softclipped_from_cigar_with_positions(original_seq, cigar)
                # Filter sequences >= SEED_THRESHOLD bp and store with position-based naming
                for position_type, seq in softclipped_seqs:
                    if len(seq) >= SEED_THRESHOLD:
                        temp_sequences.append((f"{read_name}_softclip_{position_type}", seq))

            # Save to TEMP_READS.txt
            with open('./Test_Data/TEMP_READS.txt', 'w') as f:
                for read_name, seq in temp_sequences:
                    f.write(f">{read_name}\n{seq}\n")
            
            # Step 5: Run second minimap2 alignment (iterative for ONT)
            if long_read_tech == 'ont':
                # For ONT, run iterative alignment
                run_iterative_ont_alignment()
            else:
                # For other technologies, run single second round
                cmd = build_minimap2_command('./Test_Data/TEMP_READS.txt', long_read_tech, is_initial=False)
                with open('./TEMP_SAM', 'w') as f:
                    subprocess.run(cmd, stdout=f, check=True)
                print("Round 2 minimap2 alignment completed")
            
            # Step 6: Merge results (only for non-ONT modes)
            if long_read_tech != 'ont':
                merge_temp_sam_to_output('./TEMP_SAM')
                print("Results merged back into output.sam with grouped reads")
        
        # Step 7: Convert to ViReMa format
        convert_to_virema_format()
        print(f"Minimap2 workflow completed. Output saved to: {OUTPUT_SAM}")
        
    except Exception as e:
        print(f"Error in minimap2 workflow: {e}")
        raise

def build_minimap2_command(input_file, long_read_tech=None, is_initial=True):
    """Build minimap2 command based on technology and round"""
    if long_read_tech == 'ont':
        if is_initial:
            return ['minimap2', '-ax', 'map-ont', '-Y', VIRUS_INDEX, input_file]
        else:
            return ['minimap2', '-ax', 'sr', '-k', '10', '-w', '5', '-m', '10',
                    '-n', '2', '-A', '2', '-B', '2', '-O', '2,4', '-E', '2,1',
                    '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
                    VIRUS_INDEX, input_file]
    elif long_read_tech == 'pb':
        return ['minimap2', '-ax', 'map-pb', VIRUS_INDEX, input_file]
    elif long_read_tech == 'hifi':
        return ['minimap2', '-ax', 'map-hifi', VIRUS_INDEX, input_file]
    else:
        if is_initial:
            return ['minimap2', '-ax', 'sr', '-k', '20', '-A', '1', '-B', '2',
                    '-O', '1,1', '-E', '1,1', '-r', '100', '-g', '2000',
                    '-z', '2000,1000', '-f', '0.0001', '-n', '1', '-p', '0.05',
                    '-N', '5', '-s', '5', '-t', '8', '--end-bonus', '0',
                    VIRUS_INDEX, input_file]
        else:
            return ['minimap2', '-ax', 'sr', '-k', '10', '-w', '5', '-m', '10',
                    '-n', '2', '-A', '2', '-B', '2', '-O', '2,4', '-E', '2,1',
                    '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
                    VIRUS_INDEX, input_file]


def parse_sam_file(sam_file, filter_func=None):
    """Generic SAM file parser with optional filtering"""
    results = []
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if filter_func is None or filter_func(fields):
                results.append(fields)
    return results

def parse_sam_for_softclipped(long_read_tech=None):
    """Find reads with primary alignment and softclipping (ignore supplemental alignments)"""
    # Parse SAM file for primary alignments with softclips
    def is_primary_with_softclip(fields):
        flag = int(fields[1])
        cigar = fields[5]
        return not (flag & 2048) and 'S' in cigar  # Not supplemental and has softclip
    
    alignments = parse_sam_file('output.sam', is_primary_with_softclip)
    
    # Group alignments by read name
    read_alignments = defaultdict(list)
    for fields in alignments:
        read_name = fields[0]
        read_alignments[read_name].append(fields)
    
    # Process alignments based on technology
    softclipped_reads = {}
    for read_name, alignment_list in read_alignments.items():
        if long_read_tech == 'ont':
            # For ONT, keep alignment with highest AS score
            def get_as_score(fields):
                for field in fields[11:]:
                    if field.startswith('AS:i:'):
                        return int(field.split(':')[2])
                return 0
            best_alignment = max(alignment_list, key=get_as_score)
            softclipped_reads[read_name] = {
                'line': '\t'.join(best_alignment),
                'cigar': best_alignment[5],
                'sequence': best_alignment[9]
            }
        else:
            # For other technologies, use the first alignment
            fields = alignment_list[0]
            softclipped_reads[read_name] = {
                'line': '\t'.join(fields),
                'cigar': fields[5],
                'sequence': fields[9]
            }
    
    return softclipped_reads



def extract_softclipped_from_cigar_with_positions(sequence, cigar):
    """Extract softclipped portions with position-based indices (_softclip_#, 0=before main, 1=after main)"""
    softclipped_seqs = []

    # Parse CIGAR operations
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)

    # Find the main alignment position (first M/I/=X operation)
    main_alignment_idx = None
    for i, (length, op) in enumerate(cigar_ops):
        if op in 'MI=X':
            main_alignment_idx = i
            break

    if main_alignment_idx is None:
        # No main alignment found, shouldn't happen but handle gracefully
        return []

    pos = 0
    for i, (length, op) in enumerate(cigar_ops):
        length = int(length)

        if op == 'S':  # Softclip
            seq = sequence[pos:pos + length]
            # Assign position type based on relation to main alignment
            if i < main_alignment_idx:
                position_type = 0  # Before main alignment
            else:
                position_type = 1  # After main alignment
            softclipped_seqs.append((position_type, seq))
            pos += length
        elif op in 'MI=X':  # Operations that consume query sequence
            pos += length
        # Skip operations that don't consume query sequence (D, N, H, P)

    return softclipped_seqs



def run_iterative_ont_alignment():
    """Run iterative alignment rounds for ONT reads until convergence"""
    round_num = 2
    max_rounds = 10  # Safety limit to prevent infinite loops
    current_temp_file = './Test_Data/TEMP_READS.txt'
    temp_sam_file = './TEMP_SAM'

    print("Starting iterative ONT alignment process...")

    # Run initial second round alignment
    cmd = build_minimap2_command(current_temp_file, 'ont', is_initial=False)
    with open(temp_sam_file, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    print("Round 2 minimap2 alignment completed")

    # Process alignment results and merge into output.sam
    merge_temp_sam_to_output(temp_sam_file)

    while round_num < max_rounds:
        # Extract softclips from current results
        def is_mapped_with_softclip(fields):
            flag = int(fields[1])
            cigar = fields[5]
            return not (flag & 4) and 'S' in cigar  # Mapped and has softclip
        
        alignments = parse_sam_file(temp_sam_file, is_mapped_with_softclip)
        
        new_softclips = []
        for fields in alignments:
            read_name = fields[0]
            sequence = fields[9]
            cigar = fields[5]
            
            # Parse CIGAR to find softclipped regions with positions
            softclipped_seqs = extract_softclipped_from_cigar_with_positions(sequence, cigar)
            
            # Filter sequences >= SEED_THRESHOLD bp and create nested naming
            for position_type, seq in softclipped_seqs:
                if len(seq) >= SEED_THRESHOLD:
                    nested_name = f"{read_name}_softclip_{position_type}"
                    new_softclips.append((nested_name, seq))

        if len(new_softclips) == 0:
            print(f"No new softclips found in round {round_num}, alignment complete")
            break

        # Prepare next round's input file
        round_num += 1
        next_temp_file = f'./Test_Data/TEMP_READS_R{round_num}.txt'

        # Write new softclips to next round's input file
        with open(next_temp_file, 'w') as f:
            for read_name, seq in new_softclips:
                f.write(f">{read_name}\n{seq}\n")

        print(f"Round {round_num-1} generated {len(new_softclips)} new softclips for round {round_num}")

        # Clean up previous temp files before next round
        import os
        for file_path in [temp_sam_file, current_temp_file]:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"Cleaned up {file_path}")
            except Exception as e:
                print(f"Warning: Could not remove {file_path}: {e}")
        current_temp_file = next_temp_file

        # Run alignment for this round
        print(f"Starting alignment round {round_num}")
        cmd = build_minimap2_command(current_temp_file, 'ont', is_initial=False)
        with open(temp_sam_file, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        print(f"Round {round_num} minimap2 alignment completed")

        # Process and merge results
        merge_temp_sam_to_output(temp_sam_file)

    # Final cleanup
    import os
    for file_path in [temp_sam_file, current_temp_file]:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Cleaned up {file_path}")
        except Exception as e:
            print(f"Warning: Could not remove {file_path}: {e}")
    print(f"Iterative alignment completed after {round_num-1} rounds")


def parse_softclip_path_to_position(read_name, base_read_name):
    """
    Parse softclip path into a position tuple for sequence ordering.
    
    Examples:
    - read_softclip_0 -> (0,)
    - read_softclip_1 -> (1,) 
    - read_softclip_0_softclip_1 -> (0, 1)
    - read_softclip_0_softclip_1_softclip_0 -> (0, 1, 0)
    
    Returns: tuple representing the path through the sequence
    """
    if '_softclip_' not in read_name:
        return None  # Not a softclip read
    
    # Remove base read name to get just the softclip path
    softclip_path = read_name.replace(base_read_name + '_', '')
    
    # Split by 'softclip_' and extract the indices
    parts = softclip_path.split('softclip_')
    if len(parts) < 2:
        return None
    
    # Extract position indices from the path
    position_tuple = []
    for i in range(1, len(parts)):  # Skip the first empty part
        try:
            if parts[i]:  # Make sure the part is not empty
                # Handle cases where there might be additional content after the index
                index_part = parts[i].split('_')[0] if '_' in parts[i] else parts[i]
                position_tuple.append(int(index_part))
        except ValueError:
            # If we can't parse an index, skip it
            continue
    
    return tuple(position_tuple) if position_tuple else None

def get_sequence_position_for_read(read_name, base_read_name, primary_cigar):
    """
    Determine the sequence position for a read based on softclip path or main alignment.
    
    Returns: tuple that can be used for sorting reads in sequence order
    """
    if '_softclip_' in read_name:
        return parse_softclip_path_to_position(read_name, base_read_name)
    else:
        # This is the main alignment - assign it position (0.5,) to place it between
        # softclip_0 (0,) and softclip_1 (1,)
        return (0.5,)

def select_best_alignment_by_as(alignments):
    """Select the best alignment from a list based on AS (alignment score)"""
    if not alignments:
        return None
    
    def get_as_score(alignment_line):
        fields = alignment_line.strip().split('\t')
        for field in fields[11:]:
            if field.startswith('AS:i:'):
                return int(field.split(':')[2])
        return 0
    
    # Select alignment with highest AS score
    best_alignment = max(alignments, key=get_as_score)
    return best_alignment

def merge_temp_sam_to_output(temp_sam_file):
    """Merge a single TEMP_SAM file into output.sam with AS-based selection for softclips"""
    # Load original reads
    original_headers = []
    original_reads = {}
    
    with open('output.sam', 'r') as f:
        for line in f:
            if line.startswith('@'):
                original_headers.append(line)
            else:
                read_name = line.split('\t')[0]
                if read_name not in original_reads:
                    original_reads[read_name] = []
                original_reads[read_name].append(line)
    
    # Load temp reads with base name extraction and AS-based selection
    temp_reads = {}
    try:
        with open(temp_sam_file, 'r') as f:
            for line in f:
                if not line.startswith('@'):
                    fields = line.strip().split('\t')
                    read_name = fields[0]
                    # Extract base read name from nested softclip names
                    base_read_name = read_name.split('_softclip_')[0]
                    
                    # If this is a softclip read, we need to consider all alignments (including supplementary)
                    if '_softclip_' in read_name:
                        # Extract the softclip identifier (e.g., "read1_softclip_0")
                        softclip_identifier = read_name.rsplit('_softclip_', 1)[0] + '_softclip_' + read_name.rsplit('_softclip_', 1)[1].split('_')[0]
                        
                        if base_read_name not in temp_reads:
                            temp_reads[base_read_name] = {}
                        
                        if softclip_identifier not in temp_reads[base_read_name]:
                            temp_reads[base_read_name][softclip_identifier] = []
                        
                        temp_reads[base_read_name][softclip_identifier].append(line)
                    else:
                        # Regular read (not softclip)
                        if base_read_name not in temp_reads:
                            temp_reads[base_read_name] = {}
                        if 'regular' not in temp_reads[base_read_name]:
                            temp_reads[base_read_name]['regular'] = []
                        temp_reads[base_read_name]['regular'].append(line)
    except FileNotFoundError:
        print(f"Warning: {temp_sam_file} not found")
        return

    # Select best alignments based on AS score for each softclip
    selected_temp_reads = {}
    for base_read_name, read_groups in temp_reads.items():
        selected_temp_reads[base_read_name] = []
        
        for identifier, alignments in read_groups.items():
            if identifier == 'regular':
                # For regular reads, just add all alignments
                selected_temp_reads[base_read_name].extend(alignments)
            else:
                # For softclip reads, select the best alignment based on AS score
                best_alignment = select_best_alignment_by_as(alignments)
                if best_alignment:
                    selected_temp_reads[base_read_name].append(best_alignment)

    # Write merged output
    with open('output.sam', 'w') as f:
        # Write headers first
        for line in original_headers:
            f.write(line)

        # Write reads grouped with their softclipped sequences
        for read_name in original_reads:
            # Write original alignments for this read
            for line in original_reads[read_name]:
                f.write(line)

            # Write selected softclipped alignments for this read (if any)
            if read_name in selected_temp_reads:
                for line in selected_temp_reads[read_name]:
                    f.write(line)

        # Write any remaining temp reads that don't have original alignments
        for read_name in selected_temp_reads:
            if read_name not in original_reads:
                for line in selected_temp_reads[read_name]:
                    f.write(line)






def convert_to_virema_format():
    """Convert minimap2 output to ViReMa-like format"""
    import re

    # First pass: collect all reads and identify which have softclip alignments
    all_reads = {}
    headers = []
    reads_with_softclips = set()

    with open('output.sam', 'r') as f:
        for line in f:
            if line.startswith('@'):
                headers.append(line)
                continue

            fields = line.strip().split('\t')
            read_name = fields[0]
            flag = int(fields[1])

            # Extract base read name (remove _softclip_N suffix)
            base_read_name = read_name.split('_softclip_')[0]

            if base_read_name not in all_reads:
                all_reads[base_read_name] = []
            all_reads[base_read_name].append(fields)

            # Track reads that have softclip alignments
            if '_softclip_' in read_name:
                reads_with_softclips.add(base_read_name)

    # Keep all reads for proper softclip merging
    processed_reads = all_reads

    # Second pass: filter reads based on softclip presence
    read_groups = {}
    for base_read_name, reads in processed_reads.items():
        read_groups[base_read_name] = []

        for fields in reads:
            read_name = fields[0]
            flag = int(fields[1])

            # If read has softclip alignments, keep all records (including supplemental)
            if base_read_name in reads_with_softclips:
                read_groups[base_read_name].append(fields)
            # If read has no softclip alignments, drop supplemental alignments
            elif not (flag & 2048):  # Keep primary alignments only
                read_groups[base_read_name].append(fields)

    converted_reads = []

    for base_read_name, reads in read_groups.items():
        if len(reads) == 1:
            # Single alignment - preserve softclips since they didn't map anywhere else or fell below threshold
            read = reads[0]
            converted_read = convert_single_read(read, preserve_softclips=True)
            if converted_read:
                converted_reads.append(converted_read)
        else:
            # Multiple alignments - merge into single record with gaps or paired records
            merged_result = merge_split_alignments(base_read_name, reads)
            if merged_result:
                # Check if merged_result is a list of records (paired records) or a single record
                if isinstance(merged_result, list) and len(merged_result) > 0 and isinstance(merged_result[0], list):
                    # Paired records case - list of lists
                    converted_reads.extend(merged_result)
                else:
                    # Single merged record case - list of fields
                    converted_reads.append(merged_result)

    # Final post-processing: convert internal softclips to insertions
    final_converted_reads = []
    for read in converted_reads:
        # Apply internal softclip to insertion conversion
        read_copy = read.copy()
        if read_copy[5] != '*':  # Only process reads with valid CIGAR
            cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', read_copy[5])
            if len(cigar_ops) > 2:  # Only if internal operations possible
                new_ops = []
                for i, (length, op) in enumerate(cigar_ops):
                    if op == 'S' and 0 < i < len(cigar_ops) - 1:
                        # Convert internal softclip to insertion
                        new_ops.append(f"{length}I")
                    else:
                        new_ops.append(f"{length}{op}")
                read_copy[5] = ''.join(new_ops)


        final_converted_reads.append(read_copy)

    # Write converted SAM file
    with open(OUTPUT_SAM, 'w') as f:
        # Write headers
        for header in headers:
            f.write(header)

        # Write converted reads
        for read in final_converted_reads:
            f.write('\t'.join(read) + '\n')

    print(f"Converted {len(converted_reads)} reads to ViReMa format in {OUTPUT_SAM}")

def convert_single_read(read_fields, preserve_softclips=False):
    """Convert a single minimap2 read to ViReMa format"""
    cigar = read_fields[5]
    flag = int(read_fields[1])

    # Handle unmapped reads (flag 4) - pass them through unchanged
    if flag & 4:
        return read_fields.copy()

    # Skip reads with no CIGAR information (but not unmapped)
    if cigar == '*':
        return None

    # Convert CIGAR: preserve softclips for single-location reads, otherwise convert normally
    new_cigar = convert_cigar_to_virema(cigar, preserve_softclips)

    # Create new read with simplified tags
    new_read = read_fields.copy()
    new_read[1] = '0'    # Set flag to 0 for mapped reads
    new_read[4] = '255'  # Set MAPQ to 255 for mapped reads
    new_read[5] = new_cigar

    # Keep only essential tags
    new_tags = []
    for field in read_fields[11:]:
        if field.startswith(('NM:', 'FI:', 'TC:')):
            new_tags.append(field)

    # Reconstruct read with essential fields + essential tags
    result = new_read[:11] + new_tags
    return result

def convert_cigar_to_virema(cigar, preserve_softclips=False, is_softclip_mapping=False):
    """Convert minimap2 CIGAR to ViReMa-like CIGAR"""
    # Parse CIGAR operations
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)

    new_ops = []
    for i, (length, op) in enumerate(cigar_ops):
        length = int(length)

        # Convert operations
        if op == 'S':
            if preserve_softclips or is_softclip_mapping:
                # Keep softclips as-is for single-location reads or softclip mappings
                new_ops.append(f"{length}S")
            else:
                # Skip softclips at ends, convert internal ones to skips
                if new_ops and len([o for o in cigar_ops if o[1] in 'MI=X']) > 0:
                    new_ops.append(f"{length}N")
                # Otherwise skip (leading/trailing softclips)
        elif op == 'I':
            # Preserve insertions (these are important for accurate representation)
            new_ops.append(f"{length}I")
        elif op in 'M=X':
            new_ops.append(f"{length}M")
        elif op in 'DN':
            new_ops.append(f"{length}{op}")

    return ''.join(new_ops) if new_ops else '*'


def get_as_score_from_read_fields(read_fields):
    """Extract AS score from read fields"""
    for field in read_fields[11:]:
        if field.startswith('AS:i:'):
            return int(field.split(':')[2])
    return 0

def calculate_genomic_end_position(ref_pos, cigar):
    """Calculate genomic end position from reference position and CIGAR"""
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    genomic_span = sum(int(l) for l, o in cigar_ops if o in 'MDN=X')
    return ref_pos + genomic_span - 1

def group_contiguous_segments(segments):
    """Group segments that can be merged with N gaps based on genomic positions"""
    if not segments:
        return []
    
    contiguous_groups = []
    current_group = [segments[0]]
    
    for i in range(1, len(segments)):
        curr_seg = segments[i-1]
        next_seg = segments[i]
        
        # Calculate genomic end position of current segment
        curr_end_pos = calculate_genomic_end_position(curr_seg['ref_pos'], curr_seg['cigar'])
        
        # Check merging criteria: same reference + no overlap
        can_merge = (
            curr_seg['ref_name'] == next_seg['ref_name'] and  # Same reference sequence
            curr_end_pos < next_seg['ref_pos']  # No overlap (genomically contiguous)
        )
        
        if can_merge:
            current_group.append(next_seg)
        else:
            contiguous_groups.append(current_group)
            current_group = [next_seg]
    
    contiguous_groups.append(current_group)
    return contiguous_groups

def create_single_merged_record_new(base_read_name, primary_read, segments):
    """Create single merged record with N gaps for contiguous segments"""
    import re
    
    # Build merged CIGAR by processing segments in sequence order
    merged_cigar_parts = []
    merged_ref_pos = segments[0]['ref_pos']
    merged_ref_name = segments[0]['ref_name']
    
    for i, segment in enumerate(segments):
        if i == 0:
            # First segment - add its aligned portion
            aligned_ops = extract_aligned_operations(segment['cigar'])
            merged_cigar_parts.extend(aligned_ops)
            last_end_pos = calculate_genomic_end_position(segment['ref_pos'], segment['cigar'])
        else:
            # Add N gap to next segment
            gap_size = segment['ref_pos'] - last_end_pos - 1
            if gap_size > 0:
                merged_cigar_parts.append(f"{gap_size}N")
            
            # Add segment's aligned portion
            aligned_ops = extract_aligned_operations(segment['cigar'])
            merged_cigar_parts.extend(aligned_ops)
            last_end_pos = calculate_genomic_end_position(segment['ref_pos'], segment['cigar'])
    
    # Create merged record
    merged_read = primary_read.copy()
    merged_read[0] = base_read_name
    merged_read[1] = '0'  # Primary alignment
    merged_read[2] = merged_ref_name
    merged_read[3] = str(merged_ref_pos)
    merged_read[4] = '255'  # High MAPQ
    merged_read[5] = ''.join(merged_cigar_parts)
    
    # Keep essential tags
    essential_tags = []
    for field in primary_read[11:]:
        if field.startswith(('NM:', 'FI:', 'TC:')) and not field.startswith('SA:'):
            essential_tags.append(field)
    
    return merged_read[:11] + essential_tags

def extract_aligned_operations(cigar):
    """Extract only the reference-aligned operations from CIGAR (M, D, N, =, X)"""
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    aligned_ops = []
    
    for length, op in cigar_ops:
        if op in 'M=X':
            aligned_ops.append(f"{length}M")  # Normalize to M
        elif op in 'DN':
            aligned_ops.append(f"{length}{op}")
    
    return aligned_ops

def create_paired_records_new(base_read_name, primary_read, contiguous_groups):
    """Create paired records from multiple contiguous groups"""
    records = []
    
    for group_idx, group in enumerate(contiguous_groups):
        # Build CIGAR for this group
        group_cigar_parts = []
        group_ref_pos = group[0]['ref_pos']
        group_ref_name = group[0]['ref_name']
        
        for seg_idx, segment in enumerate(group):
            if seg_idx == 0:
                # First segment in group
                aligned_ops = extract_aligned_operations(segment['cigar'])
                group_cigar_parts.extend(aligned_ops)
                last_end_pos = calculate_genomic_end_position(segment['ref_pos'], segment['cigar'])
            else:
                # Add N gap to next segment
                gap_size = segment['ref_pos'] - last_end_pos - 1
                if gap_size > 0:
                    group_cigar_parts.append(f"{gap_size}N")
                
                aligned_ops = extract_aligned_operations(segment['cigar'])
                group_cigar_parts.extend(aligned_ops)
                last_end_pos = calculate_genomic_end_position(segment['ref_pos'], segment['cigar'])
        
        # Calculate hard clip size (other groups' sequence lengths)
        other_groups_length = 0
        for other_idx, other_group in enumerate(contiguous_groups):
            if other_idx != group_idx:
                for other_seg in other_group:
                    # Estimate sequence length from aligned operations
                    other_groups_length += sum(int(l) for l, o in re.findall(r'(\d+)([MI=X])', other_seg['cigar']))
        
        # Create record for this group
        record = primary_read.copy()
        record[0] = base_read_name
        record[1] = '0' if group_idx == 0 else '2048'  # Primary vs supplementary
        record[2] = group_ref_name
        record[3] = str(group_ref_pos)
        record[4] = '255'
        record[5] = ''.join(group_cigar_parts)
        
        # Add hard clips if there are other groups
        if other_groups_length > 0:
            if group_idx == 0:
                # First group: hard clip at end
                record[5] += f"{other_groups_length}H"
            else:
                # Later groups: hard clip at beginning
                record[5] = f"{other_groups_length}H" + record[5]
        
        # Set mate information for paired records
        if len(contiguous_groups) > 1:
            if group_idx == 0:
                # Primary points to first supplementary
                next_group = contiguous_groups[1]
                record[6] = next_group[0]['ref_name']
                record[7] = str(next_group[0]['ref_pos'])
                record[8] = '0'
            else:
                record[6] = '*'
                record[7] = '0'
                record[8] = '0'
        
        # Add tags
        fi_value = group_idx + 1
        tc_value = len(contiguous_groups)
        tags = [f'FI:i:{fi_value}', 'NM:i:0', f'TC:i:{tc_value}']
        
        records.append(record[:11] + tags)
    
    return records

def merge_split_alignments(base_read_name, reads):
    """Merge multiple alignment records using softclip merging logic"""
    import re

    # Check if all reads are unmapped
    all_unmapped = all(int(read[1]) & 4 for read in reads)
    if all_unmapped:
        representative_read = reads[0].copy()
        representative_read[0] = base_read_name
        return representative_read

    # Step 1: Parse and Position - classify reads
    primary_read = None
    primary_candidates = []
    softclip_reads = {}
    
    for read in reads:
        read_name = read[0]
        flag = int(read[1])

        if '_softclip_' in read_name:
            # Parse softclip path to get sequence position
            sequence_position = get_sequence_position_for_read(read_name, base_read_name, None)
            if sequence_position is not None and flag & 4 == 0:  # Only if aligned
                # Use AS score selection for multiple alignments at same position
                if sequence_position not in softclip_reads:
                    softclip_reads[sequence_position] = read
                else:
                    current_as = get_as_score_from_read_fields(softclip_reads[sequence_position])
                    new_as = get_as_score_from_read_fields(read)
                    if new_as > current_as:
                        softclip_reads[sequence_position] = read
        elif not (flag & 4) and not (flag & 2048):  # Primary alignment candidate
            primary_candidates.append(read)
    
    # Select primary read with highest AS score
    if primary_candidates:
        primary_read = max(primary_candidates, key=get_as_score_from_read_fields)

    if not primary_read:
        return None

    # Step 2: Create segments for sequence ordering
    segments = []
    
    # Add main alignment segment
    main_ref_pos = int(primary_read[3])
    main_ref_name = primary_read[2]
    main_seq_pos = (0.5,)  # Between softclip_0 and softclip_1
    segments.append({
        'ref_pos': main_ref_pos,
        'ref_name': main_ref_name,
        'seq_pos': main_seq_pos,
        'type': 'main',
        'read': primary_read,
        'cigar': primary_read[5]
    })
    
    # Add softclip segments
    for seq_pos, softclip_read in softclip_reads.items():
        softclip_ref_pos = int(softclip_read[3])
        softclip_ref_name = softclip_read[2]
        segments.append({
            'ref_pos': softclip_ref_pos,
            'ref_name': softclip_ref_name,
            'seq_pos': seq_pos,
            'type': 'softclip',
            'read': softclip_read,
            'cigar': softclip_read[5]
        })
    
    # Step 3: Sequence ordering - sort by sequence position tuples
    segments.sort(key=lambda x: x['seq_pos'])
    
    # Step 4: Genomic merging - group contiguous segments
    contiguous_groups = group_contiguous_segments(segments)
    
    # Step 5: Output generation
    if len(contiguous_groups) == 1:
        # Single merged record
        return create_single_merged_record_new(base_read_name, primary_read, contiguous_groups[0])
    else:
        # Paired records
        return create_paired_records_new(base_read_name, primary_read, contiguous_groups)



