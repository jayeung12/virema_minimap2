#!/usr/bin/env python3

import subprocess
import re
import sys
from collections import defaultdict

# Global variables for file paths and parameters
VIRUS_INDEX = None
INPUT_DATA = None
OUTPUT_SAM = None
SEED_THRESHOLD = 15  # Default threshold

def run_minimap2_workflow(config):
    """
    Main entry point for minimap2 workflow when called from ViReMa
    Uses ViReMa configuration parameters instead of command line arguments
    """
    global VIRUS_INDEX, INPUT_DATA, OUTPUT_SAM, SEED_THRESHOLD
    
    # Set parameters from ViReMa config
    VIRUS_INDEX = config.Lib1
    INPUT_DATA = config.File1
    OUTPUT_SAM = config.Output_Dir + config.File3
    SEED_THRESHOLD = int(config.Seed) if config.Seed else 15
    
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
        run_minimap2_initial(long_read_tech)
        
        # Step 2: Parse SAM for softclipped reads with only primary alignment
        softclipped_reads = parse_sam_for_softclipped(long_read_tech)
        
        if softclipped_reads:
            # Step 3: Save to multiRound file
            save_multiround_file(softclipped_reads)
            
            # Step 4: Extract softclipped sequences
            extract_softclipped_sequences(softclipped_reads)
            
            # Step 5: Run second minimap2 alignment (iterative for ONT)
            run_minimap2_second_round(long_read_tech)
            
            # Step 6: Merge results (only for non-ONT modes)
            if long_read_tech != 'ont':
                merge_results()
        
        # Step 7: Convert to ViReMa format
        convert_to_virema_format()
        print(f"Minimap2 workflow completed. Output saved to: {OUTPUT_SAM}")
        
    except Exception as e:
        print(f"Error in minimap2 workflow: {e}")
        raise

def run_minimap2_initial(long_read_tech=None):
    """Run initial minimap2 alignment"""
    if long_read_tech == 'ont':
        cmd = [
            'minimap2', '-ax', 'map-ont', '-Y',
            VIRUS_INDEX, INPUT_DATA
        ]
    elif long_read_tech == 'pb':
        cmd = [
            'minimap2', '-ax', 'map-pb',
            VIRUS_INDEX, INPUT_DATA
        ]
    elif long_read_tech == 'hifi':
        cmd = [
            'minimap2', '-ax', 'map-hifi',
            VIRUS_INDEX, INPUT_DATA
        ]
    else:
        cmd = [
            'minimap2', '-ax', 'sr', '-k', '20', '-A', '1', '-B', '2',
            '-O', '1,1', '-E', '1,1', '-r', '100', '-g', '2000',
            '-z', '2000,1000', '-f', '0.0001', '-n', '1', '-p', '0.05',
            '-N', '5', '-s', '5', '-t', '8', '--end-bonus', '0',
            VIRUS_INDEX, INPUT_DATA
        ]

    with open('output.sam', 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    print("Initial minimap2 alignment completed")

def parse_sam_for_softclipped(long_read_tech=None):
    """Find reads with primary alignment and softclipping (ignore supplemental alignments)"""
    sam_lines = []
    primary_read_counts = defaultdict(int)
    primary_read_lines = defaultdict(str)
    all_primary_alignments = defaultdict(list)  # Store all alignments per read for ONT filtering

    with open('output.sam', 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            read_name = fields[0]
            flag = int(fields[1])
            cigar = fields[5]

            # Only count and store primary alignments (ignore supplemental alignments from initial mapping)
            if not (flag & 2048):  # Not supplemental
                primary_read_counts[read_name] += 1
                primary_read_lines[read_name] = line.strip()
                # For ONT mode, store all alignments to filter by AS score later
                if long_read_tech == 'ont':
                    all_primary_alignments[read_name].append(line.strip())

    # Find reads with primary alignment and softclipping
    softclipped_reads = {}
    for read_name, count in primary_read_counts.items():
        # For long reads, allow multiple primary alignments
        if count >= 1:  # At least one primary alignment
            line = primary_read_lines[read_name]
            fields = line.split('\t')
            cigar = fields[5]

            # Check for softclipping (S in CIGAR)
            if 'S' in cigar:
                softclipped_reads[read_name] = {
                    'line': line,
                    'cigar': cigar,
                    'sequence': fields[9]
                }

    # For ONT mode: filter to keep only the alignment with highest AS score per read
    if long_read_tech == 'ont':
        filtered_softclipped_reads = {}
        for read_name in primary_read_counts.keys():
            if read_name in all_primary_alignments:
                alignments = all_primary_alignments[read_name]
                best_alignment = None
                best_as_score = float('-inf')

                # Find alignment with highest AS score among those with softclips
                for line in alignments:
                    fields = line.split('\t')
                    cigar = fields[5]

                    # Only consider alignments with softclips
                    if 'S' in cigar:
                        # Extract AS score from optional fields
                        as_score = 0
                        for field in fields[11:]:
                            if field.startswith('AS:i:'):
                                as_score = int(field.split(':')[2])
                                break

                        # Keep alignment with highest AS score
                        if as_score > best_as_score:
                            best_as_score = as_score
                            best_alignment = line

                # Store the best alignment if found
                if best_alignment:
                    fields = best_alignment.split('\t')
                    filtered_softclipped_reads[read_name] = {
                        'line': best_alignment,
                        'cigar': fields[5],
                        'sequence': fields[9]
                    }

        return filtered_softclipped_reads

    return softclipped_reads

def save_multiround_file(softclipped_reads):
    """Save filtered SAM output to multiRound file"""
    with open('multiRound', 'w') as f:
        for read_name, data in softclipped_reads.items():
            f.write(data['line'] + '\n')
    print(f"Saved {len(softclipped_reads)} softclipped reads to multiRound")

def extract_softclipped_sequences(softclipped_reads):
    """Extract softclipped portions from SAM sequences and save to TEMP_READS.txt"""
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



def extract_softclipped_from_cigar_with_positions(sequence, cigar):
    """Extract softclipped portions with position-based indices (0=before main, 1=after main)"""
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

def run_minimap2_alignment_round(input_file, output_file, long_read_tech=None, round_num=2):
    """Run minimap2 alignment on softclipped sequences for a specific round"""
    if long_read_tech == 'ont':
        cmd = [
            'minimap2', '-ax', 'sr', '-k', '10', '-w', '5', '-m', '10',
            '-n', '2', '-A', '2', '-B', '2', '-O', '2,4', '-E', '2,1',
            '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
            VIRUS_INDEX, input_file
        ]
    elif long_read_tech == 'pb':
        cmd = [
            'minimap2', '-ax', 'map-pb',
            VIRUS_INDEX, input_file
        ]
    elif long_read_tech == 'hifi':
        cmd = [
            'minimap2', '-ax', 'map-hifi',
            VIRUS_INDEX, input_file
        ]
    else:
        cmd = [
            'minimap2', '-ax', 'sr', '-k', '10', '-w', '5', '-m', '10',
            '-n', '2', '-A', '2', '-B', '2', '-O', '2,4', '-E', '2,1',
            '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
            VIRUS_INDEX, input_file
        ]

    with open(output_file, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    print(f"Round {round_num} minimap2 alignment completed")

def extract_softclips_from_sam(sam_file, round_num):
    """Extract softclipped sequences from SAM alignment results for next round"""
    temp_sequences = []

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            read_name = fields[0]
            flag = int(fields[1])
            cigar = fields[5]
            sequence = fields[9]

            # Skip unmapped reads
            if flag & 4:
                continue

            # Check for softclipping (S in CIGAR)
            if 'S' in cigar:
                # Parse CIGAR to find softclipped regions with positions
                softclipped_seqs = extract_softclipped_from_cigar_with_positions(sequence, cigar)

                # Filter sequences >= SEED_THRESHOLD bp and create nested naming
                for position_type, seq in softclipped_seqs:
                    if len(seq) >= SEED_THRESHOLD:
                        # Create nested softclip name: original_softclip_X becomes original_softclip_X_softclip_Y
                        nested_name = f"{read_name}_softclip_{position_type}"
                        temp_sequences.append((nested_name, seq))

    return temp_sequences

def cleanup_temp_files(*files):
    """Clean up temporary files"""
    import os
    for file_path in files:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Cleaned up {file_path}")
        except Exception as e:
            print(f"Warning: Could not remove {file_path}: {e}")

def run_iterative_ont_alignment():
    """Run iterative alignment rounds for ONT reads until convergence"""
    round_num = 2
    max_rounds = 10  # Safety limit to prevent infinite loops
    current_temp_file = './Test_Data/TEMP_READS.txt'
    temp_sam_file = './TEMP_SAM'

    print("Starting iterative ONT alignment process...")

    # Run initial second round alignment
    run_minimap2_alignment_round(current_temp_file, temp_sam_file, 'ont', 2)

    # Process alignment results and merge into output.sam
    merge_temp_sam_to_output(temp_sam_file)

    while round_num < max_rounds:
        # Extract softclips from current results
        new_softclips = extract_softclips_from_sam(temp_sam_file, round_num)

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
        cleanup_temp_files(temp_sam_file, current_temp_file)
        current_temp_file = next_temp_file

        # Run alignment for this round
        print(f"Starting alignment round {round_num}")
        run_minimap2_alignment_round(current_temp_file, temp_sam_file, 'ont', round_num)

        # Process and merge results
        merge_temp_sam_to_output(temp_sam_file)

    # Final cleanup
    cleanup_temp_files(temp_sam_file, current_temp_file)
    print(f"Iterative alignment completed after {round_num-1} rounds")

def merge_temp_sam_to_output(temp_sam_file):
    """Merge a single TEMP_SAM file into output.sam"""
    # Read current output.sam content
    original_lines = []
    original_reads = {}

    with open('output.sam', 'r') as f:
        for line in f:
            if line.startswith('@'):
                original_lines.append(line)
            else:
                read_name = line.split('\t')[0]
                if read_name not in original_reads:
                    original_reads[read_name] = []
                original_reads[read_name].append(line)

    # Read new alignments from temp SAM file
    temp_reads = {}
    try:
        with open(temp_sam_file, 'r') as f:
            for line in f:
                if not line.startswith('@'):
                    read_name = line.split('\t')[0]
                    # Extract base read name (handle nested softclip names)
                    base_read_name = extract_base_read_name(read_name)
                    if base_read_name not in temp_reads:
                        temp_reads[base_read_name] = []
                    temp_reads[base_read_name].append(line)
    except FileNotFoundError:
        print(f"Warning: {temp_sam_file} not found")
        return

    # Write merged output
    with open('output.sam', 'w') as f:
        # Write headers first
        for line in original_lines:
            f.write(line)

        # Write reads grouped with their softclipped sequences
        for read_name in original_reads:
            # Write original alignments for this read
            for line in original_reads[read_name]:
                f.write(line)

            # Write softclipped alignments for this read (if any)
            if read_name in temp_reads:
                for line in temp_reads[read_name]:
                    f.write(line)

        # Write any remaining temp reads that don't have original alignments
        for read_name in temp_reads:
            if read_name not in original_reads:
                for line in temp_reads[read_name]:
                    f.write(line)

def extract_base_read_name(read_name):
    """Extract the original base read name from nested softclip names"""
    # Handle nested softclip names like: SRR5085928.1_softclip_0_softclip_1
    # Should return: SRR5085928.1
    parts = read_name.split('_softclip_')
    return parts[0]

def run_minimap2_second_round(long_read_tech=None):
    """Run second minimap2 alignment on softclipped sequences"""
    if long_read_tech == 'ont':
        # For ONT, run iterative alignment
        run_iterative_ont_alignment()
    else:
        # For other technologies, run single second round
        run_minimap2_alignment_round('./Test_Data/TEMP_READS.txt', './TEMP_SAM', long_read_tech, 2)

def merge_results():
    """Add TEMP_SAM results back to output.sam with proper grouping"""
    # Read original alignments
    original_lines = []
    original_reads = {}
    with open('output.sam', 'r') as f:
        for line in f:
            if line.startswith('@'):
                original_lines.append(line)
            else:
                read_name = line.split('\t')[0]
                if read_name not in original_reads:
                    original_reads[read_name] = []
                original_reads[read_name].append(line)

    # Read new alignments
    temp_reads = {}
    with open('./TEMP_SAM', 'r') as f:
        for line in f:
            if not line.startswith('@'):
                read_name = line.split('\t')[0]
                # Extract base read name (remove _softclip_N suffix)
                base_read_name = '_'.join(read_name.split('_')[:-2]) if '_softclip_' in read_name else read_name
                if base_read_name not in temp_reads:
                    temp_reads[base_read_name] = []
                temp_reads[base_read_name].append(line)

    # Write merged output with grouped reads
    with open('output.sam', 'w') as f:
        # Write headers first
        for line in original_lines:
            f.write(line)

        # Write reads grouped with their softclipped sequences
        for read_name in original_reads:
            # Write original alignments for this read
            for line in original_reads[read_name]:
                f.write(line)

            # Write softclipped alignments for this read (if any)
            if read_name in temp_reads:
                for line in temp_reads[read_name]:
                    f.write(line)

        # Write any remaining temp reads that don't have original alignments
        for read_name in temp_reads:
            if read_name not in original_reads:
                for line in temp_reads[read_name]:
                    f.write(line)

    print("Results merged back into output.sam with grouped reads")

def integrate_nested_softclips_properly(all_reads):
    """Properly integrate nested softclips into their direct parents only"""
    import re

    processed_reads = {}

    for base_read_name, reads in all_reads.items():
        # Group reads by hierarchy level
        original_reads = []
        softclip_alignments = {}  # full_name -> alignment

        for fields in reads:
            read_name = fields[0]
            flag = int(fields[1])

            if '_softclip_' not in read_name:
                # Original read (including supplemental alignments)
                original_reads.append(fields)
            else:
                # Store all softclip alignments by their full name
                softclip_alignments[read_name] = fields

        # Process nested softclips: only integrate children into their direct parents
        final_softclips = {}

        # First, identify which softclips have children
        parent_child_map = {}  # parent_name -> [child_names]

        for softclip_name in softclip_alignments.keys():
            # Count softclip levels
            softclip_count = softclip_name.count('_softclip_')

            if softclip_count > 1:
                # This is a nested softclip - find its direct parent
                parts = softclip_name.split('_softclip_')
                parent_name = '_softclip_'.join(parts[:-1])  # Remove last level
                child_idx = int(parts[-1])

                if parent_name not in parent_child_map:
                    parent_child_map[parent_name] = {}
                parent_child_map[parent_name][child_idx] = softclip_alignments[softclip_name]

        # Now process each first-level softclip
        for softclip_name, alignment in softclip_alignments.items():
            softclip_count = softclip_name.count('_softclip_')

            if softclip_count == 1:
                # This is a first-level softclip
                if softclip_name in parent_child_map:
                    # This softclip has children - integrate them
                    integrated_alignment = integrate_children_into_softclip(alignment, parent_child_map[softclip_name])
                    final_softclips[softclip_name] = integrated_alignment
                else:
                    # No children - keep as is
                    final_softclips[softclip_name] = alignment

        # Build final read list: original reads + processed first-level softclips only
        final_reads = original_reads.copy()
        for softclip_name, alignment in final_softclips.items():
            final_reads.append(alignment)

        processed_reads[base_read_name] = final_reads

    return processed_reads

def integrate_children_into_softclip(parent_softclip, child_alignments):
    """Integrate child alignments into a parent softclip alignment"""
    import re

    if not child_alignments:
        return parent_softclip

    parent_fields = parent_softclip.copy()
    parent_cigar = parent_fields[5]
    parent_pos = int(parent_fields[3])

    # Parse parent CIGAR to find where softclips would be
    parent_ops = re.findall(r'(\d+)([MIDNSHPX=])', parent_cigar)

    # Find the main alignment in parent CIGAR
    main_alignment_idx = None
    for i, (length, op) in enumerate(parent_ops):
        if op in 'MI=X':
            main_alignment_idx = i
            break

    if main_alignment_idx is None:
        return parent_softclip  # No main alignment found

    # Build new CIGAR by potentially replacing end softclips with child alignments
    new_cigar_parts = []
    current_ref_pos = parent_pos

    for i, (length, op) in enumerate(parent_ops):
        length = int(length)

        if op == 'S':
            # Determine if this is a leading (0) or trailing (1) softclip
            if i < main_alignment_idx:
                softclip_position = 0  # Leading
            else:
                softclip_position = 1  # Trailing

            # Check if we have a child alignment for this position
            if softclip_position in child_alignments:
                child_fields = child_alignments[softclip_position]
                child_cigar = child_fields[5]
                child_pos = int(child_fields[3])

                # Skip unmapped children
                if int(child_fields[1]) & 4:
                    new_cigar_parts.append(f"{length}S")
                    continue

                # Calculate gap to child alignment
                if softclip_position == 0:
                    # Leading softclip: gap from child end to parent start
                    child_aligned_length = sum(int(l) for l, o in re.findall(r'(\d+)([MI=X])', child_cigar))
                    child_end = child_pos + child_aligned_length
                    gap_size = current_ref_pos - child_end
                else:
                    # Trailing softclip: gap from parent end to child start
                    gap_size = child_pos - current_ref_pos

                # Add child alignment
                child_virema_cigar = convert_softclip_cigar_to_virema(child_cigar)
                if child_virema_cigar and child_virema_cigar != '*':
                    if softclip_position == 0:
                        # For leading: child + gap + (parent continues)
                        new_cigar_parts.insert(0, child_virema_cigar)
                        if gap_size > 0:
                            new_cigar_parts.insert(1, f"{gap_size}N")
                    else:
                        # For trailing: (parent continues) + gap + child
                        if gap_size > 0:
                            new_cigar_parts.append(f"{gap_size}N")
                        new_cigar_parts.append(child_virema_cigar)
                else:
                    # Child didn't map properly, keep original softclip
                    new_cigar_parts.append(f"{length}S")
            else:
                # No child alignment for this softclip position
                new_cigar_parts.append(f"{length}S")
        else:
            # Non-softclip operation
            if op in 'M=X':
                new_cigar_parts.append(f"{length}M")
                current_ref_pos += length
            elif op in 'DN':
                new_cigar_parts.append(f"{length}{op}")
                if op == 'D':
                    current_ref_pos += length
            else:
                new_cigar_parts.append(f"{length}{op}")

    # Update the parent alignment
    parent_fields[5] = ''.join(new_cigar_parts)

    return parent_fields

def convert_internal_softclips_to_insertions(cigar):
    """Convert internal softclips (not at ends) to insertions"""
    import re

    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    if len(cigar_ops) <= 2:
        return cigar  # No internal operations possible

    new_ops = []
    for i, (length, op) in enumerate(cigar_ops):
        if op == 'S':
            # Check if this is an internal softclip (not first or last)
            if 0 < i < len(cigar_ops) - 1:
                # Convert internal softclip to insertion
                new_ops.append(f"{length}I")
            else:
                # Keep end softclips as-is
                new_ops.append(f"{length}S")
        else:
            new_ops.append(f"{length}{op}")

    return ''.join(new_ops)


def convert_to_virema_format():
    """Convert minimap2 output.sam to ViReMa-like format in output_conv.sam"""
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
            base_read_name = extract_base_read_name(read_name)

            if base_read_name not in all_reads:
                all_reads[base_read_name] = []
            all_reads[base_read_name].append(fields)

            # Track reads that have softclip alignments
            if '_softclip_' in read_name:
                reads_with_softclips.add(base_read_name)

    # Process nested softclip integration BEFORE the main conversion
    processed_reads = integrate_nested_softclips_properly(all_reads)

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
            read_copy[5] = convert_internal_softclips_to_insertions(read_copy[5])


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

def convert_cigar_to_virema(cigar, preserve_softclips=False):
    """Convert minimap2 CIGAR to ViReMa-like CIGAR"""
    # Parse CIGAR operations
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)

    new_ops = []
    for length, op in cigar_ops:
        length = int(length)

        # Convert operations
        if op == 'S':
            if preserve_softclips:
                # Keep softclips as-is for single-location reads
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

def convert_softclip_cigar_to_virema(softclip_cigar):
    """Convert softclip mapping CIGAR to ViReMa format, preserving trailing softclips and insertions"""
    import re

    # Parse CIGAR operations
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', softclip_cigar)

    new_ops = []
    for length, op in cigar_ops:
        length = int(length)

        # Convert operations
        if op == 'I':
            # Preserve insertions (these are important for accurate representation)
            new_ops.append(f"{length}I")
        elif op in 'M=X':
            new_ops.append(f"{length}M")
        elif op == 'S':
            # Preserve softclips (these are remaining unmapped portions)
            new_ops.append(f"{length}S")
        elif op in 'DN':
            new_ops.append(f"{length}{op}")

    return ''.join(new_ops) if new_ops else '*'

def merge_split_alignments(base_read_name, reads):
    """Merge multiple alignment records into ViReMa-style record(s)"""
    import re

    # Check if all reads are unmapped - if so, return the first one as representative
    all_unmapped = all(int(read[1]) & 4 for read in reads)
    if all_unmapped:
        representative_read = reads[0].copy()
        representative_read[0] = base_read_name  # Use base read name
        return representative_read

    # Classify alignment records
    primary_read = None
    supplemental_reads = []
    softclip_reads = {}

    for read in reads:
        read_name = read[0]
        flag = int(read[1])

        if '_softclip_' in read_name:
            # Extract softclip index (0 for beginning, 1 for end, etc.)
            softclip_idx = int(read_name.split('_softclip_')[1])
            if flag & 4 == 0:  # Only if aligned
                softclip_reads[softclip_idx] = read
        elif flag & 2048:  # Supplemental alignment
            supplemental_reads.append(read)
        elif not (flag & 4):  # Primary alignment
            primary_read = read

    if not primary_read:
        return None



    # Parse primary alignment CIGAR to understand softclip structure
    primary_cigar = primary_read[5]
    primary_pos = int(primary_read[3])

    # Parse CIGAR operations to find softclips
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', primary_cigar)

    # Track which softclips have alignments and which don't
    unmapped_softclips = {}

    # Check if we can merge into a single record or need paired records
    segments = []
    current_ref_pos = primary_pos

    for i, (length, op) in enumerate(cigar_ops):
        length = int(length)

        if op == 'S':
            # Find the main alignment position (first M/I/=X operation)
            main_alignment_idx = None
            for j, (_, op_check) in enumerate(cigar_ops):
                if op_check in 'MI=X':
                    main_alignment_idx = j
                    break

            # Determine softclip index based on position relative to main alignment
            if main_alignment_idx is not None:
                if i < main_alignment_idx:
                    softclip_idx = 0  # Before main alignment
                else:
                    softclip_idx = 1  # After main alignment
            else:
                # No main alignment found - this shouldn't happen for reads that passed filtering
                # Skip this softclip since we can't determine its relative position
                continue

            if softclip_idx in softclip_reads:
                softclip_read = softclip_reads[softclip_idx]
                softclip_pos = int(softclip_read[3])
                softclip_cigar = softclip_read[5]

                # Calculate aligned length from softclip CIGAR (for sorting/gap calculation)
                softclip_aligned = sum(int(l) for l, o in re.findall(r'(\d+)([MI=X])', softclip_cigar))

                # Store both aligned length and full CIGAR for proper reconstruction
                segments.append((softclip_pos, softclip_aligned, 'softclip', i, softclip_idx, softclip_cigar))
            else:
                # Track unmapped softclips to preserve them
                unmapped_softclips[i] = (length, 'S')

        elif op in 'MI=X':
            # Main alignment segment
            segments.append((current_ref_pos, length, 'main', i, -1, f"{length}M"))
            current_ref_pos += length

        elif op == 'D':
            current_ref_pos += length
        elif op == 'N':
            current_ref_pos += length

    # Check if segments map to different reference sequences
    ref_sequences = set()
    ref_sequences.add(primary_read[2])  # Add primary read's reference

    for read in reads:
        if not read[0].startswith(base_read_name + '_softclip'):  # Skip softclip reads for now
            ref_sequences.add(read[2])

    # Also check softclip reads for their reference sequences
    for softclip_idx, softclip_read in softclip_reads.items():
        ref_sequences.add(softclip_read[2])

    different_references = len(ref_sequences) > 1

    # Check if segments are in genomic order, accounting for overlaps (only relevant for same reference)
    in_genomic_order = True
    if len(segments) >= 2:
        # Sort segments by CIGAR position to identify primary vs softclip order
        segments_with_cigar_pos = []
        for seg in segments:
            seg_pos, seg_length, seg_type, cigar_pos, softclip_idx, seg_cigar = seg
            segments_with_cigar_pos.append((seg_pos, seg_length, seg_type, cigar_pos, softclip_idx, seg_cigar))

        # Sort by CIGAR position to get read order
        segments_with_cigar_pos.sort(key=lambda x: x[3])  # Sort by cigar_pos

        # Check each adjacent pair for duplication/overlap
        for i in range(len(segments_with_cigar_pos) - 1):
            curr_pos, curr_len, curr_type, curr_cigar_pos, curr_softclip_idx, curr_cigar = segments_with_cigar_pos[i]
            next_pos, next_len, next_type, next_cigar_pos, next_softclip_idx, next_cigar = segments_with_cigar_pos[i + 1]

            # Calculate end positions using MDN=X operations (same as duplication detection)
            curr_end = curr_pos + sum(int(l) for l, o in re.findall(r'(\d+)([MDN=X])', curr_cigar)) - 1
            next_end = next_pos + sum(int(l) for l, o in re.findall(r'(\d+)([MDN=X])', next_cigar)) - 1

            # Check for overlap/duplication patterns
            if curr_type == 'softclip' and curr_softclip_idx == 0 and next_type == 'main':
                # Leading softclip before main - duplication if softclip end >= main start
                if curr_end >= next_pos:
                    in_genomic_order = False
            elif curr_type == 'main' and next_type == 'softclip' and next_softclip_idx == 1:
                # Main before trailing softclip - duplication if softclip start <= main end
                if next_pos <= curr_end:
                    in_genomic_order = False

    # Prioritize different references - always create paired records when softclips map to different references
    if different_references or (not in_genomic_order and len(segments) >= 2):
        # Create paired records for different references or out-of-order segments
        return create_paired_records(base_read_name, primary_read, segments, softclip_reads, unmapped_softclips)
    elif unmapped_softclips or any(seg[2] == 'softclip' for seg in segments):
        # If there are unmapped softclips OR any mapped softclips, create a single merged record that handles them
        return create_single_merged_record_with_softclips(base_read_name, primary_read, segments, unmapped_softclips)
    else:
        # Create single merged record for same reference and in-order segments
        return create_single_merged_record(base_read_name, primary_read, segments)

def create_paired_records(base_read_name, primary_read, segments, softclip_reads, unmapped_softclips=None):
    """Create paired records like ViReMa format for out-of-order segments or different references"""
    import re

    # Identify which segment corresponds to main alignment vs softclip
    main_segment = None
    softclip_segment = None

    for seg in segments:
        pos, aligned_length, seg_type, cigar_idx, softclip_idx = seg[:5]  # Take only first 5 elements
        if seg_type == 'main':
            main_segment = seg
        elif seg_type == 'softclip':
            softclip_segment = seg

    if not main_segment or not softclip_segment:
        # Fallback to original logic if we can't identify both segments
        segments_by_pos = sorted(segments, key=lambda x: x[0])
        return create_fallback_paired_records(base_read_name, primary_read, segments_by_pos, softclip_reads, unmapped_softclips)

    # Determine order based on CIGAR index (position in read sequence)
    main_cigar_idx = main_segment[3]
    softclip_cigar_idx = softclip_segment[3]

    # If softclip comes before main alignment in read sequence, swap the roles
    if softclip_cigar_idx < main_cigar_idx:
        # Softclip segment should be the primary record, main segment should be supplementary
        first_segment = softclip_segment
        second_segment = main_segment
        first_is_softclip = True
    else:
        # Main segment should be the primary record, softclip segment should be supplementary
        first_segment = main_segment
        second_segment = softclip_segment
        first_is_softclip = False

    records = []

    # Get softclip read info
    softclip_idx = softclip_segment[4]
    softclip_read = softclip_reads.get(softclip_idx)

    if not softclip_read:
        return create_fallback_paired_records(base_read_name, primary_read, segments, softclip_reads, unmapped_softclips)

    # Create records based on read sequence order (first_segment gets primary flag, second gets supplementary)
    if first_is_softclip:
        # Softclip comes first in read sequence - it becomes the primary record
        first_read_data = softclip_read
        first_ref = softclip_read[2]
        first_pos = softclip_read[3]
        second_read_data = primary_read
        second_ref = primary_read[2]
        second_pos = primary_read[3]
    else:
        # Main alignment comes first in read sequence - it becomes the primary record
        first_read_data = primary_read
        first_ref = primary_read[2]
        first_pos = primary_read[3]
        second_read_data = softclip_read
        second_ref = softclip_read[2]
        second_pos = softclip_read[3]

    # Build CIGAR and extract sequences for both records
    primary_cigar = primary_read[5]
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', primary_cigar)
    if unmapped_softclips is None:
        unmapped_softclips = {}

    # Calculate sequence segments from original read
    seq_segments = []
    seq_pos = 0

    for i, (length, op) in enumerate(cigar_ops):
        length = int(length)
        if op in 'SMI=X':
            segment_seq = primary_read[9][seq_pos:seq_pos + length]
            segment_qual = primary_read[10][seq_pos:seq_pos + length] if primary_read[10] != '*' else ''

            if op == 'S':
                if i == softclip_cigar_idx:
                    # This is the mapped softclip segment
                    seg_type = 'mapped_softclip'
                elif i in unmapped_softclips:
                    # This is an unmapped softclip segment
                    seg_type = 'unmapped_softclip'
                else:
                    seg_type = 'unmapped_softclip'  # Default for softclips
            else:
                # This is main alignment
                seg_type = 'main'

            seq_segments.append((seg_type, length, segment_seq, segment_qual, i, op))
            seq_pos += length
        elif op in 'DN':
            # These operations don't consume query sequence but are part of main alignment
            seq_segments.append(('main', length, '', '', i, op))

    # Build first record (primary)
    record1 = primary_read.copy()
    record1[0] = base_read_name
    record1[1] = '0'  # Primary alignment flag
    record1[2] = first_ref
    record1[3] = first_pos
    record1[4] = '255'  # MAPQ
    record1[6] = second_ref  # RNEXT
    record1[7] = second_pos  # PNEXT
    record1[8] = '0'  # TLEN

    # Build second record (supplementary)
    record2 = primary_read.copy()
    record2[0] = base_read_name
    record2[1] = '2048'  # Supplementary alignment flag
    record2[2] = second_ref
    record2[3] = second_pos
    record2[4] = '255'  # MAPQ
    record2[6] = '*'
    record2[7] = '0'
    record2[8] = '0'

    # Build CIGAR and sequences based on segment order
    if first_is_softclip:
        # First record: softclip alignment + hard clips for other parts
        first_cigar_parts = []
        first_seq_parts = []
        first_qual_parts = []
        second_cigar_parts = []
        second_seq_parts = []
        second_qual_parts = []

        # Get the actual CIGAR from the softclip alignment record
        softclip_cigar = softclip_read[5]
        softclip_seq = softclip_read[9]
        softclip_qual = softclip_read[10] if softclip_read[10] != '*' else ''
        
        # Use the softclip alignment CIGAR directly for the first record
        first_cigar_parts.append(softclip_cigar)
        first_seq_parts.append(softclip_seq)
        if softclip_qual:
            first_qual_parts.append(softclip_qual)
        
        # Calculate lengths for hard clipping
        softclip_seq_length = len(softclip_seq)
        main_seq_length = len(primary_read[9]) - softclip_seq_length
        
        # Process other segments for the second record
        for seg_type, length, seq, qual, idx, op in seq_segments:
            if seg_type == 'mapped_softclip':
                # Skip - already handled above with actual softclip CIGAR
                continue
                        
            elif seg_type == 'main':
                # This goes to second record - preserve all original CIGAR operations
                second_cigar_parts.append(f"{length}{op}")
                if op in 'SMI=X':
                    second_seq_parts.append(seq)
                    if qual:
                        second_qual_parts.append(qual)
                        
            elif seg_type == 'unmapped_softclip':
                # Preserve in first record, but determine which record based on position
                if idx < softclip_cigar_idx:
                    # This unmapped softclip belongs to the main alignment portion
                    second_cigar_parts.append(f"{length}S")
                    second_seq_parts.append(seq)
                    if qual:
                        second_qual_parts.append(qual)
                else:
                    # This unmapped softclip belongs to the softclip alignment portion  
                    first_cigar_parts.append(f"{length}S")
                    first_seq_parts.append(seq)
                    if qual:
                        first_qual_parts.append(qual)

        # Add hard clips representing the other portion of the read
        if main_seq_length > 0:
            first_cigar_parts.append(f"{main_seq_length}H")
        if softclip_seq_length > 0:
            second_cigar_parts.insert(0, f"{softclip_seq_length}H")

        record1[5] = ''.join(first_cigar_parts)
        record1[9] = ''.join(first_seq_parts)
        record1[10] = ''.join(first_qual_parts) if first_qual_parts else '*'

        record2[5] = ''.join(second_cigar_parts)
        record2[9] = ''.join(second_seq_parts)
        record2[10] = ''.join(second_qual_parts) if second_qual_parts else '*'
    else:
        # Original logic - main alignment first
        first_cigar_parts = []
        first_seq_parts = []
        first_qual_parts = []
        second_cigar_parts = []
        second_seq_parts = []
        second_qual_parts = []

        # Get softclip data from actual alignment record
        softclip_cigar = softclip_read[5]
        softclip_seq = softclip_read[9]
        softclip_qual = softclip_read[10] if softclip_read[10] != '*' else ''
        
        # Calculate lengths for hard clipping
        softclip_seq_length = len(softclip_seq)
        main_seq_length = len(primary_read[9]) - softclip_seq_length

        for seg_type, length, seq, qual, idx, op in seq_segments:
            if seg_type == 'main':
                # This goes to first record - preserve all original CIGAR operations
                first_cigar_parts.append(f"{length}{op}")
                if op in 'SMI=X':
                    first_seq_parts.append(seq)
                    if qual:
                        first_qual_parts.append(qual)
                        
            elif seg_type == 'mapped_softclip':
                # This goes to second record - use actual softclip CIGAR
                if not second_cigar_parts:  # Only add once
                    second_cigar_parts.append(softclip_cigar)
                    second_seq_parts.append(softclip_seq)
                    if softclip_qual:
                        second_qual_parts.append(softclip_qual)
                        
            elif seg_type == 'unmapped_softclip':
                # Preserve in appropriate record based on position
                if idx < softclip_cigar_idx:
                    # This unmapped softclip belongs to the main alignment portion
                    first_cigar_parts.append(f"{length}S")
                    first_seq_parts.append(seq)
                    if qual:
                        first_qual_parts.append(qual)
                else:
                    # This unmapped softclip belongs to the softclip alignment portion  
                    second_cigar_parts.append(f"{length}S")
                    second_seq_parts.append(seq)
                    if qual:
                        second_qual_parts.append(qual)

        # Add hard clips representing the other portion of the read
        if softclip_seq_length > 0:
            first_cigar_parts.append(f"{softclip_seq_length}H")
        if main_seq_length > 0:
            second_cigar_parts.insert(0, f"{main_seq_length}H")

        record1[5] = ''.join(first_cigar_parts)
        record1[9] = ''.join(first_seq_parts)
        record1[10] = ''.join(first_qual_parts) if first_qual_parts else '*'

        record2[5] = ''.join(second_cigar_parts)
        record2[9] = ''.join(second_seq_parts)
        record2[10] = ''.join(second_qual_parts) if second_qual_parts else '*'

    records.append(record1[:11] + ['FI:i:1', 'NM:i:0', 'TC:i:2'])
    records.append(record2[:11] + ['FI:i:2', 'NM:i:0', 'TC:i:2'])

    return records

def create_fallback_paired_records(base_read_name, primary_read, segments_by_pos, softclip_reads, unmapped_softclips=None):
    """Fallback function for creating paired records when segments can't be clearly identified"""
    records = []

    for i, seg in enumerate(segments_by_pos):
        pos, aligned_length, seg_type, cigar_idx, softclip_idx = seg[:5]
        record = primary_read.copy()
        record[0] = base_read_name
        record[3] = str(pos)
        record[4] = '255'  # MAPQ like ViReMa

        if i == 0:
            # First record
            record[1] = '0'  # Primary alignment flag
            record[5] = f"{aligned_length}M{sum(s[1] for s in segments_by_pos[1:])}H"  # Main + hard clip
            if len(segments_by_pos) > 1:
                second_seg = segments_by_pos[1]
                if second_seg[2] == 'softclip' and second_seg[4] in softclip_reads:
                    record[6] = softclip_reads[second_seg[4]][2]  # Reference from softclip read
                else:
                    record[6] = record[2]  # Same reference
                record[7] = str(second_seg[0])  # PNEXT
            else:
                record[6] = '*'
                record[7] = '0'
            record[8] = '0'  # TLEN
            record[10] = record[10][:aligned_length] if len(record[10]) >= aligned_length else record[10]
            new_tags = ['FI:i:1', 'NM:i:0', 'TC:i:2']
        else:
            # Second record
            record[1] = '2048'  # Supplementary alignment flag
            record[5] = f"{sum(s[1] for s in segments_by_pos[:-1])}H{aligned_length}M"  # Hard clip + main
            record[6] = '*'
            record[7] = '0'
            record[8] = '0'
            if seg_type == 'softclip' and softclip_idx in softclip_reads:
                softclip_read = softclip_reads[softclip_idx]
                record[2] = softclip_read[2]  # Update reference sequence
                record[9] = softclip_read[9]
                record[10] = softclip_read[10] if softclip_read[10] != '*' else record[10][-aligned_length:]
            new_tags = ['FI:i:2', 'NM:i:0', 'TC:i:2']

        result = record[:11] + new_tags
        records.append(result)

    return records

def create_single_merged_record(base_read_name, primary_read, segments):
    """Create single merged record with gaps for in-order segments"""
    # Sort segments by read sequence order
    segments.sort(key=lambda x: x[3])

    # Build merged CIGAR
    merged_cigar = ""
    last_end = 0
    start_pos = segments[0][0] if segments else int(primary_read[3])

    for i, seg in enumerate(segments):
        pos, aligned_length, seg_type, cigar_idx, softclip_idx = seg[:5]
        if i == 0:
            merged_cigar += f"{aligned_length}M"
            last_end = pos + aligned_length
        else:
            gap_size = pos - last_end
            if gap_size > 0:
                merged_cigar += f"{gap_size}N"
            merged_cigar += f"{aligned_length}M"
            last_end = pos + aligned_length

    if not merged_cigar:
        return None

    # Create merged read
    merged_read = primary_read.copy()
    merged_read[0] = base_read_name
    merged_read[1] = '0'    # Set flag to 0 for mapped reads
    merged_read[3] = str(start_pos)
    merged_read[4] = '255'  # Set MAPQ to 255 for mapped reads
    merged_read[5] = merged_cigar

    # Keep only essential tags
    new_tags = []
    for field in primary_read[11:]:
        if field.startswith(('NM:', 'FI:', 'TC:')) and not field.startswith('SA:'):
            new_tags.append(field)

    result = merged_read[:11] + new_tags
    return result

def calculate_genomic_span(cigar):
    """Calculate the genomic span of a CIGAR string (reference bases consumed)"""
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    genomic_span = 0
    for length, op in cigar_ops:
        length = int(length)
        if op in 'MD=X':  # Operations that consume reference bases
            genomic_span += length
    return genomic_span

def create_single_merged_record_with_softclips(base_read_name, primary_read, segments, unmapped_softclips):
    """Create single merged record preserving unmapped softclips"""
    import re


    # Parse original CIGAR to understand structure
    primary_cigar = primary_read[5]
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', primary_cigar)
    primary_pos = int(primary_read[3])

    # Build mapping of CIGAR index to segment info
    mapped_segments = {}
    for seg in segments:
        if seg[2] == 'softclip':
            cigar_idx = seg[3]  # Use the CIGAR index directly from the segment
            mapped_segments[cigar_idx] = seg

    # Build new CIGAR by processing operations in order
    new_cigar_ops = []
    start_pos = None
    main_alignment_length = 0

    # First pass: calculate full genomic span of main alignment
    for length, op in cigar_ops:
        length = int(length)
        if op in 'MD=X':  # Operations that consume reference bases
            main_alignment_length += length

    for i, (length, op) in enumerate(cigar_ops):
        length = int(length)

        if op == 'S':
            if i in unmapped_softclips:
                # Preserve unmapped softclip as-is
                new_cigar_ops.append(f"{length}S")
            elif i in mapped_segments:
                # This softclip has an alignment - use its full CIGAR with gap
                seg = mapped_segments[i]
                softclip_pos = seg[0]
                softclip_length = seg[1]
                softclip_cigar = seg[5]  # Full CIGAR from softclip mapping

                # Set overall start position
                if start_pos is None:
                    if i == 0:  # Leading softclip
                        start_pos = softclip_pos
                    else:  # Trailing softclip
                        start_pos = primary_pos

                if i == 0:  # Leading softclip
                    # Leading softclip maps first - convert its CIGAR
                    converted_softclip_cigar = convert_softclip_cigar_to_virema(softclip_cigar)
                    new_cigar_ops.append(converted_softclip_cigar)
                else:  # Trailing softclip
                    # Calculate gap between end of primary and start of softclip
                    gap_size = softclip_pos - (primary_pos + main_alignment_length)
                    if gap_size > 0:
                        new_cigar_ops.append(f"{gap_size}N")
                    # Convert the softclip CIGAR to preserve internal softclips
                    converted_softclip_cigar = convert_softclip_cigar_to_virema(softclip_cigar)
                    new_cigar_ops.append(converted_softclip_cigar)
            else:
                # Unmapped softclip - preserve as S
                new_cigar_ops.append(f"{length}S")

        elif op in 'MI=XDN':
            # Main alignment operations - preserve complex CIGAR structure
            if start_pos is None:
                start_pos = primary_pos

            # Check if we need a gap before the main alignment (leading softclip case)
            if new_cigar_ops and 0 in mapped_segments and op in 'MI=X':
                # Only add gap before first match/insert operation, not for every operation
                first_main_op = True
                for prev_i in range(i):
                    if cigar_ops[prev_i][1] in 'MI=X':
                        first_main_op = False
                        break

                if first_main_op:
                    seg = mapped_segments[0]
                    softclip_pos = seg[0]
                    softclip_cigar = seg[5]
                    # Calculate actual genomic span of softclip
                    softclip_genomic_span = calculate_genomic_span(softclip_cigar)
                    gap_size = primary_pos - (softclip_pos + softclip_genomic_span)
                    if gap_size > 0:
                        new_cigar_ops.append(f"{gap_size}N")

            # Preserve the operation as-is for main alignment
            if op in 'M=X':
                new_cigar_ops.append(f"{length}M")
            else:
                new_cigar_ops.append(f"{length}{op}")

    # Use calculated start position or fall back to original
    if start_pos is None:
        start_pos = primary_pos


    # Create merged read
    merged_read = primary_read.copy()
    merged_read[0] = base_read_name
    merged_read[1] = '0'    # Set flag to 0 for mapped reads
    merged_read[3] = str(start_pos)
    merged_read[4] = '255'  # Set MAPQ to 255 for mapped reads
    merged_read[5] = ''.join(new_cigar_ops) if new_cigar_ops else '*'

    # Keep only essential tags
    new_tags = []
    for field in primary_read[11:]:
        if field.startswith(('NM:', 'FI:', 'TC:')) and not field.startswith('SA:'):
            new_tags.append(field)

    result = merged_read[:11] + new_tags
    return result

def create_single_record_with_gaps(base_read_name, primary_read, supplemental_read, softclip_reads):
    """Create single record with N gaps for primary + supplemental + softclip alignments"""
    import re


    # Parse primary CIGAR to understand softclip structure
    primary_cigar = primary_read[5]
    primary_pos = int(primary_read[3])
    supplemental_pos = int(supplemental_read[3])

    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', primary_cigar)

    # Find softclip positions and extract main alignment block
    leading_softclip = None
    trailing_softclip = None
    main_alignment_cigar = []
    main_alignment_ref_length = 0

    for i, (length, op) in enumerate(cigar_ops):
        length = int(length)
        if op == 'S':
            if i == 0:  # Leading softclip
                leading_softclip = length
            else:  # Trailing softclip (assume it's at the end)
                trailing_softclip = length
        elif op in 'MI=XDN':
            # This is part of the main alignment - preserve the full structure
            if op in 'MI=X':
                main_alignment_cigar.append(f"{length}M")
            else:
                main_alignment_cigar.append(f"{length}{op}")

            # Track reference length for gap calculation (only M, D, N consume reference)
            if op in 'MDN=X':
                main_alignment_ref_length += length

    main_alignment_cigar_str = ''.join(main_alignment_cigar)

    # Determine which softclip mapped and calculate gap
    new_cigar_parts = []
    start_pos = primary_pos

    # Check if we have a leading softclip that mapped
    if leading_softclip and 0 in softclip_reads:
        softclip_read = softclip_reads[0]
        softclip_pos = int(softclip_read[3])
        softclip_cigar = softclip_read[5]
        softclip_genomic_span = calculate_genomic_span(softclip_cigar)
        softclip_converted_cigar = convert_softclip_cigar_to_virema(softclip_cigar)

        # Calculate gap between softclip end and primary start
        gap_size = primary_pos - (softclip_pos + softclip_genomic_span)
        if gap_size > 0:
            new_cigar_parts = [softclip_converted_cigar, f"{gap_size}N", main_alignment_cigar_str]
            start_pos = softclip_pos
        else:
            # No gap needed, sequences are adjacent
            new_cigar_parts = [softclip_converted_cigar, main_alignment_cigar_str]
            start_pos = softclip_pos

        # Handle trailing softclip if present but unmapped
        if trailing_softclip:
            new_cigar_parts.append(f"{trailing_softclip}S")

    # Check if we have a trailing softclip that mapped
    elif trailing_softclip and 1 in softclip_reads:
        softclip_read = softclip_reads[1]
        softclip_pos = int(softclip_read[3])
        softclip_cigar = softclip_read[5]
        softclip_genomic_span = calculate_genomic_span(softclip_cigar)
        softclip_converted_cigar = convert_softclip_cigar_to_virema(softclip_cigar)

        # Calculate gap between primary end and softclip start
        gap_size = softclip_pos - (primary_pos + main_alignment_ref_length)
        if gap_size > 0:
            new_cigar_parts = [main_alignment_cigar_str, f"{gap_size}N", softclip_converted_cigar]
        else:
            # No gap needed, sequences are adjacent
            new_cigar_parts = [main_alignment_cigar_str, softclip_converted_cigar]

        # Handle leading softclip if present but unmapped
        if leading_softclip:
            new_cigar_parts.insert(0, f"{leading_softclip}S")

    else:
        # Fallback: just use primary alignment with preserved softclips
        return convert_single_read(primary_read, preserve_softclips=True)

    # Create merged record
    merged_read = primary_read.copy()
    merged_read[0] = base_read_name
    merged_read[1] = '0'    # Set flag to 0 for mapped reads
    merged_read[3] = str(start_pos)
    merged_read[4] = '255'  # Set MAPQ to 255 for mapped reads
    merged_read[5] = ''.join(new_cigar_parts)

    # Keep only essential tags
    new_tags = []
    for field in primary_read[11:]:
        if field.startswith(('NM:', 'FI:', 'TC:')) and not field.startswith('SA:'):
            new_tags.append(field)

    result = merged_read[:11] + new_tags
    return result

