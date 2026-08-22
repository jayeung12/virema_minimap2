#!/usr/bin/env python3

import subprocess
import re
import sys
from collections import defaultdict

# Using Python implementation with cached optimizations
print("Using Python implementation with cached optimizations")

# Global variables for file paths and parameters
VIRUS_INDEX = None
INPUT_DATA = None
OUTPUT_SAM = None
SEED_THRESHOLD = 25  # Default threshold
MICROINDEL_THRESHOLD = 0  # Default threshold (below this, an indel is ordinary noise)
CHUNK_SIZE = 1000000  # Default chunk size for memory-efficient processing (1M reads, matching ViReMa)
THREADS = '1'  # Default thread count

# Long-read technologies that get the full multi-round rescue treatment (embedded-insertion
# rewrite every round, iterative softclip realignment up to max_rounds, segment-pair
# preservation, redundant-segment dedup). All of this logic is itself tech-agnostic -- it
# operates on CIGAR/SAM content, not on which aligner preset produced it -- so extending it
# beyond 'ont' is just a matter of routing 'pb'/'hifi' through the same path instead of the
# single-round fallback used for plain short-read mode.
ITERATIVE_RESCUE_TECHS = ('ont', 'pb', 'hifi')

def run_minimap2_workflow(config):
    """
    Main entry point for minimap2 workflow when called from ViReMa
    Uses ViReMa configuration parameters
    """
    global VIRUS_INDEX, INPUT_DATA, OUTPUT_SAM, SEED_THRESHOLD, MICROINDEL_THRESHOLD, CHUNK_SIZE, THREADS

    # Set parameters from ViReMa config
    VIRUS_INDEX = config.Lib1
    INPUT_DATA = config.File1
    OUTPUT_SAM = config.Output_Dir + config.File3
    SEED_THRESHOLD = int(config.Seed) if config.Seed else 25
    MICROINDEL_THRESHOLD = int(config.MicroInDel_Length) if config.MicroInDel_Length else 0
    THREADS = config.Threads if config.Threads else '1'
    
    # Use ViReMa's chunk size if available, otherwise use default
    if hasattr(config, 'Chunk') and config.Chunk:
        CHUNK_SIZE = int(config.Chunk)
    else:
        CHUNK_SIZE = 1000000  # Default 1M reads
    
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

        # Step 1b: Convert any large embedded insertions into explicit soft-clips so the
        # existing rescue pipeline below finds them like any other soft-clipped read (see
        # rewrite_embedded_insertions_as_softclips docstring for why this is needed).
        if long_read_tech in ITERATIVE_RESCUE_TECHS:
            rewrite_embedded_insertions_as_softclips('output.sam', MICROINDEL_THRESHOLD, SEED_THRESHOLD)

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
                ins_len = data.get('ins_len')

                # Parse CIGAR to find softclipped regions with their positions relative to main alignment
                softclipped_seqs = extract_softclipped_from_cigar_with_positions(original_seq, cigar)
                # Filter sequences >= SEED_THRESHOLD bp and store with position-based naming
                for position_type, seq in softclipped_seqs:
                    if len(seq) >= SEED_THRESHOLD:
                        temp_sequences.extend(
                            _softclip_candidates_for_record(read_name, position_type, seq, ins_len)
                        )

            # Save to TEMP_READS.txt
            with open('./Test_Data/TEMP_READS.txt', 'w') as f:
                for read_name, seq in temp_sequences:
                    f.write(f">{read_name}\n{seq}\n")
            
            # Step 5: Run second minimap2 alignment (iterative for ont/pb/hifi)
            if long_read_tech in ITERATIVE_RESCUE_TECHS:
                run_iterative_ont_alignment(long_read_tech)
            else:
                # Plain short-read mode: run single second round
                cmd = build_minimap2_command('./Test_Data/TEMP_READS.txt', long_read_tech, is_initial=False)
                with open('./TEMP_SAM', 'w') as f:
                    subprocess.run(cmd, stdout=f, check=True)
                print("Round 2 minimap2 alignment completed")

            # Step 6: Merge results (only for the single-round fallback -- the iterative path
            # merges internally after every round)
            if long_read_tech not in ITERATIVE_RESCUE_TECHS:
                merge_temp_sam_to_output('./TEMP_SAM')
                print("Results merged back into output.sam with grouped reads")
        
        # Step 7: Convert to ViReMa format
        convert_to_virema_format()
        print(f"Minimap2 workflow completed. Output saved to: {OUTPUT_SAM}")
        
        # Clear cache to free memory after processing
        clear_softclip_cache()
        
    except Exception as e:
        print(f"Error in minimap2 workflow: {e}")
        # Clear cache even on error to free memory
        clear_softclip_cache()
        raise

def build_minimap2_command(input_file, long_read_tech=None, is_initial=True):
    """Build minimap2 command based on technology and round"""
    if long_read_tech == 'ont':
        if is_initial:
            return ['minimap2', '-a', '-k', '15', '-w', '5',
                    '-A', '1', '-B', '2', '-O', '2,32', '-E', '1,0',
                    '-z', '200', '-g', '2000', '-Y',
                    '-t', THREADS, VIRUS_INDEX, input_file]
        else:
            return ['minimap2', '-ax', 'sr', '-k', '10', '-w', '5', '-m', '10',
                    '-n', '2', '-A', '2', '-B', '2', '-O', '2,4', '-E', '2,1',
                    '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
                    '-t', THREADS, VIRUS_INDEX, input_file]
    elif long_read_tech == 'pb':
        if is_initial:
            return ['minimap2', '-ax', 'map-pb', '-t', THREADS, VIRUS_INDEX, input_file]
        else:
            # Round 2+ rescue: realigning a short extracted fragment, not a whole read --
            # map-pb (tuned for whole long reads) fails outright on short queries (confirmed
            # empirically: a 167bp fragment came back completely unmapped under map-pb but
            # aligned cleanly, MAPQ 60, under this preset). Same base as ONT's round-2+ preset,
            # since CLR and ONT are both ~85%-accuracy long-read technologies facing the same
            # short-fragment-sensitivity problem -- reusing already-validated parameters here
            # rather than inventing new, untested ones for a comparable error regime.
            return ['minimap2', '-ax', 'sr', '-k', '10', '-w', '5', '-m', '10',
                    '-n', '2', '-A', '2', '-B', '2', '-O', '2,4', '-E', '2,1',
                    '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
                    '-t', THREADS, VIRUS_INDEX, input_file]
    elif long_read_tech == 'hifi':
        if is_initial:
            return ['minimap2', '-ax', 'map-hifi', '-t', THREADS, VIRUS_INDEX, input_file]
        else:
            # Round 2+ rescue, same rationale as the 'pb' branch above -- map-hifi fails
            # outright on short fragments for the same reason map-pb does. Same short-read
            # base preset as ONT/PB's round-2+, but with a larger k (15 vs 10): HiFi's much
            # lower error rate (~96.5%+ vs ~85%) doesn't need ONT/CLR's extra seed sensitivity,
            # and a larger, more specific minimizer reduces spurious short-repeat matches.
            return ['minimap2', '-ax', 'sr', '-k', '15', '-w', '5', '-m', '10',
                    '-n', '2', '-A', '2', '-B', '2', '-O', '2,4', '-E', '2,1',
                    '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
                    '-t', THREADS, VIRUS_INDEX, input_file]
    else:
        if is_initial:
            return ['minimap2', '-ax', 'sr', '-k', '20', '-A', '1', '-B', '2',
                    '-O', '2,8', '-g', '2000',
                    '-z', '800,400', '-n', '1', '-p', '0.3',
                    '-N', '3', '-s', '20', '-t', THREADS, '--end-bonus', '0',
                    VIRUS_INDEX, input_file]
        else:
            return ['minimap2', '-ax', 'sr', '-k', '10', '-w', '5', '-m', '10',
                    '-n', '2', '-A', '1', '-B', '2', '-O', '12,32', '-E', '2,1',
                    '--end-bonus', '5', '-s', '20', '-z', '200', '-r', '50',
                    '-t', THREADS, VIRUS_INDEX, input_file]


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
        if long_read_tech in ITERATIVE_RESCUE_TECHS:
            # For ont/pb/hifi, keep alignment with highest AS score
            def get_as_score(fields):
                for field in fields[11:]:
                    if field.startswith('AS:i:'):
                        return int(field.split(':')[2])
                return 0
            best_alignment = max(alignment_list, key=get_as_score)
            softclipped_reads[read_name] = {
                'line': '\t'.join(best_alignment),
                'cigar': best_alignment[5],
                'sequence': best_alignment[9],
                'ins_len': _get_ins_len_tag(best_alignment)
            }
        else:
            # For other technologies, use the first alignment
            fields = alignment_list[0]
            softclipped_reads[read_name] = {
                'line': '\t'.join(fields),
                'cigar': fields[5],
                'sequence': fields[9],
                'ins_len': _get_ins_len_tag(fields)
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


def _get_ins_len_tag(fields):
    """Read back the ZI:i: tag rewrite_embedded_insertions_as_softclips stamps on a record it
    rewrote -- the length of the originally-embedded insertion, before it got folded into a
    trailing soft-clip together with everything after it. None if absent (record wasn't
    rewritten, or predates this tag)."""
    for field in fields[11:]:
        if field.startswith('ZI:i:'):
            return int(field.split(':')[2])
    return None


def _softclip_candidates_for_record(read_name, position_type, seq, ins_len):
    """Round-2+ candidates to emit for one soft-clipped piece of a record.

    When `ins_len` is set (this soft-clip came from rewrite_embedded_insertions_as_softclips,
    which always clips everything from the insertion point through the end of the read), split
    it into two INDEPENDENT candidates -- the insertion content alone, and everything after it
    alone -- instead of realigning them glued together as one fragment.

    Why this matters: the glued fragment is only unambiguous when the insertion's true origin
    and the flank that follows it fall in increasing reference order (confirmed empirically for
    cross-locus insertions where the donor sits upstream of the acceptor -- "earlier" case).
    When the donor sits downstream of the acceptor ("later" case), the glued fragment contains a
    real backward jump in the middle: minimap2 aligns it as one confident, single placement
    anchored by whichever piece is longer/more informative (usually the flank), and the shorter
    piece gets dragged along and mis-anchored -- often smeared into ordinary-looking small
    indels rather than left as a soft-clip, so the existing recursive rescue (which watches for
    exactly that signal) never notices anything is wrong. Aligning the insertion by itself
    removes the contamination in both directions regardless of which way it points.

    Falls back to the original single-glued-candidate behavior when there's no tagged insertion
    length, the clip isn't the trailing (position_type == 1) kind rewrite always produces, or
    either resulting piece would be shorter than SEED_THRESHOLD (nothing to gain by splitting a
    piece too short to independently seed anyway).
    """
    if ins_len is None or position_type != 1:
        return [(f"{read_name}_softclip_{position_type}", seq)]
    ins_seq = seq[:ins_len]
    rest_seq = seq[ins_len:]
    if len(ins_seq) < SEED_THRESHOLD or len(rest_seq) < SEED_THRESHOLD:
        return [(f"{read_name}_softclip_{position_type}", seq)]
    return [
        (f"{read_name}_softclip_{position_type}_ins", ins_seq),
        (f"{read_name}_softclip_{position_type}", rest_seq),
    ]


def rewrite_embedded_insertions_as_softclips(sam_file, microindel_threshold, seed_threshold):
    """Convert one large embedded CIGAR 'I' op per record into an explicit trailing soft-clip,
    so the existing soft-clip rescue pipeline (parse_sam_for_softclipped /
    run_iterative_ont_alignment) discovers it on its own without any other changes.

    Why: minimap2 often represents a genuine two-locus chimeric read (e.g. a fragment spliced
    in from a different reference segment) as ONE alignment with the whole fragment folded
    into a single embedded 'I' op, rather than splitting into primary+supplementary records --
    there is then no soft-clip left for the existing rescue mechanism to ever notice. Turning
    "...M<bigI>M..." into "...M<S covering bigI + everything after it>" makes it look exactly
    like an ordinary soft-clipped read: the clipped piece starts right at the insertion
    boundary, and when minimap2 realigns just that piece fresh, empirically (see loc2_S1_16 in
    the FHV benchmark this was built for) it reliably produces a correct primary+supplementary
    split, because there's now ~0 flanking context on one side to tempt it into re-embedding.

    Deletions ('D'/'N' ops) are never touched -- only 'I' -- so this cannot affect deletion
    detection.

    Threshold: an insertion must be strictly longer than `microindel_threshold` (so it's not
    ordinary small-scale ONT noise -- MicroInDel_Length is ViReMa's existing "how small an
    indel counts as noise" parameter, reused here for the same purpose) AND at least
    `seed_threshold` long (so the clipped-off piece is actually long enough to be worth
    re-aligning at all -- otherwise this would truncate a perfectly good alignment for a
    fragment that the existing `len(seq) >= SEED_THRESHOLD` gate would just discard anyway).

    Only the first qualifying insertion in a record is rewritten; if more remain in whatever
    gets clipped off, they get caught the same way when that piece's own alignment is scanned
    in a later round.
    """
    min_len = max(microindel_threshold + 1, seed_threshold)

    try:
        with open(sam_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return

    out_lines = []
    for line in lines:
        if line.startswith('@'):
            out_lines.append(line)
            continue

        fields = line.rstrip('\n').split('\t')
        flag = int(fields[1])
        cigar = fields[5]
        if (flag & 4) or cigar == '*':
            out_lines.append(line)
            continue

        ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
        split_idx = None
        for i, (length, op) in enumerate(ops):
            if op == 'I' and int(length) >= min_len:
                split_idx = i
                break

        if split_idx is None:
            out_lines.append(line)
            continue

        kept_ops = ops[:split_idx]
        # Need real aligned content before the split to keep as a meaningful primary alignment
        if not any(op in 'M=X' for _, op in kept_ops):
            out_lines.append(line)
            continue

        ins_len = int(ops[split_idx][0])
        clip_len = sum(int(length) for length, op in ops[split_idx:] if op in 'MIS=X')
        new_cigar = ''.join(f'{length}{op}' for length, op in kept_ops) + f'{clip_len}S'

        fields[5] = new_cigar
        # Tag with the qualifying insertion's own length so the round-2+ extraction step can
        # test it in isolation instead of only ever gluing it to what follows it (see
        # _softclip_candidates_for_record).
        fields.append(f'ZI:i:{ins_len}')
        out_lines.append('\t'.join(fields) + '\n')

    with open(sam_file, 'w') as f:
        f.writelines(out_lines)


def run_iterative_ont_alignment(long_read_tech='ont'):
    """Run iterative alignment rounds (ont/pb/hifi) until convergence"""
    round_num = 2
    max_rounds = 6  # Safety limit to prevent infinite loops -- successful rescues observed to
                     # converge by round 3-4, so this bounds worst-case cost for reads that
                     # never resolve without affecting reads that do
    current_temp_file = './Test_Data/TEMP_READS.txt'
    temp_sam_file = './TEMP_SAM'

    print(f"Starting iterative {long_read_tech} alignment process...")

    # Run initial second round alignment
    cmd = build_minimap2_command(current_temp_file, long_read_tech, is_initial=False)
    with open(temp_sam_file, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    print("Round 2 minimap2 alignment completed")
    rewrite_embedded_insertions_as_softclips(temp_sam_file, MICROINDEL_THRESHOLD, SEED_THRESHOLD)

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
            ins_len = _get_ins_len_tag(fields)

            # Parse CIGAR to find softclipped regions with positions
            softclipped_seqs = extract_softclipped_from_cigar_with_positions(sequence, cigar)

            # Filter sequences >= SEED_THRESHOLD bp and create nested naming
            for position_type, seq in softclipped_seqs:
                if len(seq) >= SEED_THRESHOLD:
                    new_softclips.extend(
                        _softclip_candidates_for_record(read_name, position_type, seq, ins_len)
                    )

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
        cmd = build_minimap2_command(current_temp_file, long_read_tech, is_initial=False)
        with open(temp_sam_file, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        print(f"Round {round_num} minimap2 alignment completed")
        rewrite_embedded_insertions_as_softclips(temp_sam_file, MICROINDEL_THRESHOLD, SEED_THRESHOLD)

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


def calculate_sequence_position(read_name, base_read_name):
    """
    Calculate sequence position for a softclip read based on its path.
    
    The position represents where in the sequence this segment should appear:
    - Main alignment is at position 0.5
    - softclip_0 variants appear before their parent (subtract offset)
    - softclip_1 variants appear after their parent (add offset)
    - Each nesting level uses smaller offsets to maintain proper ordering
    
    Examples:
    - main: 0.5
    - main_softclip_0: 0.4 (before main)
    - main_softclip_0_softclip_1: 0.45 (between softclip_0 and main)
    - main_softclip_1: 0.6 (after main)
    - main_softclip_1_softclip_0: 0.55 (between main and softclip_1)
    
    Returns: float representing sequence position, or None if not a softclip read
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
    path_indices = []
    for i in range(1, len(parts)):  # Skip the first empty part
        try:
            if parts[i]:  # Make sure the part is not empty
                # Handle cases where there might be additional content after the index
                index_part = parts[i].split('_')[0] if '_' in parts[i] else parts[i]
                path_indices.append(int(index_part))
        except ValueError:
            # If we can't parse an index, skip it
            continue
    
    if not path_indices:
        return None
    
    # Calculate position by walking the path
    # Start from main alignment position
    position = 0.5
    
    # Each step in the path adjusts position
    # Use decreasing offsets for each level to maintain proper ordering
    base_offset = 0.1
    
    for level, direction in enumerate(path_indices):
        # Calculate offset for this level (smaller for deeper nesting)
        level_offset = base_offset / (2 ** level)
        
        if direction == 0:
            # softclip_0: move to earlier position (before parent)
            position -= level_offset
        elif direction == 1:
            # softclip_1: move to later position (after parent)
            position += level_offset
    
    return position

def get_sequence_position_for_read(read_name, base_read_name, primary_cigar):
    """
    Determine the sequence position for a read based on softclip path or main alignment.
    
    Returns: float that can be used for sorting reads in sequence order
    """
    if '_softclip_' in read_name:
        return calculate_sequence_position(read_name, base_read_name)
    else:
        # This is the main alignment - assign it position 0.5 to place it between
        # softclip_0 variants (< 0.5) and softclip_1 variants (> 0.5)
        return 0.5

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

def select_alignments_for_identifier(alignments):
    """Select which SAM lines to keep for one softclip identifier.

    Keeps the best-scoring PRIMARY record plus any SUPPLEMENTARY records (flag 2048) for the
    same identifier -- those are complementary parts of the SAME chimeric alignment (e.g. when
    rewrite_embedded_insertions_as_softclips exposes a genuine two-locus junction and this
    identifier's own realignment splits into primary+supplementary), not competing
    alternatives, so collapsing them down to a single "best" record would silently discard a
    real segment. Genuinely competing SECONDARY records (flag 256) are dropped, same as before.
    """
    def flag_of(line):
        return int(line.strip().split('\t')[1])

    primaries = [a for a in alignments if not (flag_of(a) & (256 | 2048))]
    supplementaries = [a for a in alignments if flag_of(a) & 2048]

    if not primaries:
        # No primary among the candidates (shouldn't normally happen) -- fall back to the
        # original single-best-AS behavior rather than guess.
        best = select_best_alignment_by_as(alignments)
        return [best] if best else []

    best_primary = select_best_alignment_by_as(primaries)
    return [best_primary] + supplementaries

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
                        # Extract the softclip identifier (e.g., "read1_softclip_0"). An "_ins"
                        # suffix (see _softclip_candidates_for_record) marks the isolated-insertion
                        # half of a split candidate pair -- keep it as its own identifier, distinct
                        # from its "rest" sibling that shares the same digit, so AS-based selection
                        # doesn't collapse the pair down to a single winner (both must survive to
                        # merge_split_alignments, which pairs them explicitly).
                        _suffix = read_name.rsplit('_softclip_', 1)[1]
                        _digit = _suffix.split('_')[0]
                        _ins_marker = '_ins' if _suffix.endswith('_ins') else ''
                        softclip_identifier = read_name.rsplit('_softclip_', 1)[0] + '_softclip_' + _digit + _ins_marker
                        
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
                # For softclip reads, keep the best primary plus any supplementary companions
                # (see select_alignments_for_identifier docstring for why)
                selected_temp_reads[base_read_name].extend(select_alignments_for_identifier(alignments))

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
    """Convert minimap2 output to ViReMa-like format with chunked processing"""
    import re
    import os
    
    # Check file size for memory optimization guidance
    try:
        file_size = os.path.getsize('output.sam')
        if file_size > 500 * 1024 * 1024:  # 500MB
            print(f"Processing large SAM file ({file_size // (1024*1024):,} MB) using chunked processing...")
    except OSError:
        pass

    # First pass: collect all reads and identify which have softclip alignments
    all_reads = {}
    headers = []
    reads_with_softclips = set()

    print("Loading SAM file and identifying reads...")
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

    print(f"Loaded {len(all_reads):,} base reads.")

    # Keep all reads for proper softclip merging (but can clear the original reference)
    processed_reads = all_reads
    all_reads = None  # Clear reference to help with memory

    # Second pass: filter reads based on softclip presence
    print("Filtering reads and preparing for conversion...")
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

    # Clear processed_reads after filtering to save memory
    processed_reads = None
    
    # Memory-efficient chunked processing and writing
    print("Converting reads to ViReMa format...")
    final_converted_reads = []
    total_converted = 0
    chunk_count = 0
    
    # Open output file for writing
    with open(OUTPUT_SAM, 'w') as output_file:
        # Write headers first
        for header in headers:
            output_file.write(header)
        
        # Process reads in chunks
        current_chunk = []
        base_read_names = list(read_groups.keys())
        
        for i, base_read_name in enumerate(base_read_names):
            reads = read_groups[base_read_name]
            
            if len(reads) == 1:
                # Single alignment - preserve softclips since they didn't map anywhere else or fell below threshold
                read = reads[0]
                converted_read = convert_single_read(read, preserve_softclips=True)
                if converted_read:
                    current_chunk.append(converted_read)
            else:
                # Multiple alignments - merge into single record with gaps or paired records
                merged_result = merge_split_alignments(base_read_name, reads)
                if merged_result:
                    # Check if merged_result is a list of records (paired records) or a single record
                    if isinstance(merged_result, list) and len(merged_result) > 0 and isinstance(merged_result[0], list):
                        # Paired records case - list of lists
                        current_chunk.extend(merged_result)
                    else:
                        # Single merged record case - list of fields
                        current_chunk.append(merged_result)
            
            # Process chunk when it reaches the specified size or at the end
            if len(current_chunk) >= CHUNK_SIZE or i == len(base_read_names) - 1:
                if current_chunk:
                    chunk_count += 1
                    print(f"Processing chunk {chunk_count} with {len(current_chunk):,} reads...")
                    
                    # Apply final post-processing: convert internal softclips to insertions
                    for read in current_chunk:
                        # Update hardclipped sequences and quality scores first
                        update_hardclipped_values(read)
                        
                        # Apply internal softclip to insertion conversion
                        if read[5] != '*':  # Only process reads with valid CIGAR
                            cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', read[5])
                            if len(cigar_ops) > 2:  # Only if internal operations possible
                                new_ops = []
                                for j, (length, op) in enumerate(cigar_ops):
                                    if op == 'S' and 0 < j < len(cigar_ops) - 1:
                                        # Convert internal softclip to insertion
                                        new_ops.append(f"{length}I")
                                    else:
                                        new_ops.append(f"{length}{op}")
                                read[5] = ''.join(new_ops)
                        
                        # Write read immediately to avoid memory accumulation
                        output_file.write('\t'.join(read) + '\n')
                        total_converted += 1
                    
                    # Clear chunk to free memory
                    current_chunk = []

    print(f"Converted {total_converted:,} reads to ViReMa format in {OUTPUT_SAM}")

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
    """Create single merged record with N gaps for contiguous segments, using SAM-based accounting"""
    import re
    
    # Build hierarchical accounting map based on actual SAM file data
    segment_accounting = build_hierarchical_accounting_map(base_read_name, segments)
    
    # Build merged CIGAR by processing segments in sequence order
    merged_cigar_parts = []
    merged_ref_pos = segments[0]['ref_pos']
    merged_ref_name = segments[0]['ref_name']
    
    for i, segment in enumerate(segments):
        is_first_segment = (i == 0)
        is_last_segment = (i == len(segments) - 1)
        
        # Get accounting info for this specific segment
        segment_read_name = segment['read'][0]
        accounted_positions = segment_accounting.get(segment_read_name, set())
        
        if i == 0:
            # First segment
            aligned_ops = extract_aligned_operations(segment['cigar'], accounted_positions, is_first_segment, is_last_segment)
            merged_cigar_parts.extend(aligned_ops)
            last_end_pos = calculate_genomic_end_position(segment['ref_pos'], segment['cigar'])
        else:
            # Add N gap to next segment
            gap_size = segment['ref_pos'] - last_end_pos - 1
            if gap_size > 0:
                merged_cigar_parts.append(f"{gap_size}N")
            
            # Add segment's aligned portion
            aligned_ops = extract_aligned_operations(segment['cigar'], accounted_positions, is_first_segment, is_last_segment)
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

def extract_aligned_operations(cigar, accounted_softclip_positions=None, is_first_segment=False, is_last_segment=False):
    """Extract reference-aligned operations from CIGAR, handling softclips based on global sequence position"""
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    aligned_ops = []
    
    # Find main alignment position to determine softclip positions
    main_alignment_idx = None
    for i, (length, op) in enumerate(cigar_ops):
        if op in 'MI=X':
            main_alignment_idx = i
            break
    
    if accounted_softclip_positions is None:
        accounted_softclip_positions = set()
    
    for i, (length, op) in enumerate(cigar_ops):
        if op in 'M=X':
            aligned_ops.append(f"{length}M")  # Normalize to M
        elif op in 'DN':
            aligned_ops.append(f"{length}{op}")
        elif op == 'I':
            # Always preserve insertion operations
            aligned_ops.append(f"{length}I")
        elif op == 'S':
            # Determine softclip position type using scalar positions
            if main_alignment_idx is not None:
                if i < main_alignment_idx:
                    position_scalar = 0.4  # Before main alignment (softclip_0 position)
                else:
                    position_scalar = 0.6  # After main alignment (softclip_1 position)
                
                # Only process softclips that are NOT accounted for by softclip alignments
                if position_scalar not in accounted_softclip_positions:
                    # Determine if this is a global edge or internal softclip
                    is_global_edge = False
                    
                    if i == 0 and is_first_segment:
                        # Beginning of first segment in global sequence
                        is_global_edge = True
                    elif i == len(cigar_ops) - 1 and is_last_segment:
                        # End of last segment in global sequence
                        is_global_edge = True
                    
                    if is_global_edge:
                        # Global edge softclip - preserve as softclip
                        aligned_ops.append(f"{length}S")
                    else:
                        # Internal softclip - convert to insertion
                        aligned_ops.append(f"{length}I")
                # If softclip is accounted for by a softclip alignment, skip it completely
    
    return aligned_ops

def create_paired_records_new(base_read_name, primary_read, contiguous_groups):
    """Create paired records from multiple contiguous groups using SAM-based accounting"""
    records = []
    
    # Build hierarchical accounting map based on actual SAM file data
    all_segments = []
    for group in contiguous_groups:
        all_segments.extend(group)
    segment_accounting = build_hierarchical_accounting_map(base_read_name, all_segments)
    
    # Pre-process all groups to get their final CIGAR operations
    processed_groups = []
    
    for group_idx, group in enumerate(contiguous_groups):
        # Build CIGAR for this group
        group_cigar_parts = []
        group_ref_pos = group[0]['ref_pos']
        group_ref_name = group[0]['ref_name']
        
        for seg_idx, segment in enumerate(group):
            # Determine global position: first segment overall and last segment overall
            global_first = (group_idx == 0 and seg_idx == 0)
            global_last = (group_idx == len(contiguous_groups) - 1 and seg_idx == len(group) - 1)
            
            # Get accounting info for this specific segment
            segment_read_name = segment['read'][0]
            accounted_positions = segment_accounting.get(segment_read_name, set())
            
            if seg_idx == 0:
                # First segment in group
                aligned_ops = extract_aligned_operations(segment['cigar'], accounted_positions, global_first, global_last)
                group_cigar_parts.extend(aligned_ops)
                last_end_pos = calculate_genomic_end_position(segment['ref_pos'], segment['cigar'])
            else:
                # Add N gap to next segment
                gap_size = segment['ref_pos'] - last_end_pos - 1
                if gap_size > 0:
                    group_cigar_parts.append(f"{gap_size}N")
                
                # Add segment's aligned portion
                aligned_ops = extract_aligned_operations(segment['cigar'], accounted_positions, global_first, global_last)
                group_cigar_parts.extend(aligned_ops)
                last_end_pos = calculate_genomic_end_position(segment['ref_pos'], segment['cigar'])
        
        # Store processed group information
        processed_groups.append({
            'cigar_parts': group_cigar_parts,
            'ref_pos': group_ref_pos,
            'ref_name': group_ref_name
        })
    
    # Now build records with correct hard clips calculated from processed operations
    for group_idx, processed_group in enumerate(processed_groups):
        # Calculate directional hard clip sizes from processed CIGAR operations
        preceding_groups_length = 0
        following_groups_length = 0
        
        # Calculate length of groups that come before this one (for beginning hard clip)
        for idx in range(group_idx):
            processed_cigar = ''.join(processed_groups[idx]['cigar_parts'])
            # Calculate query-consuming operations from processed CIGAR
            preceding_groups_length += sum(int(l) for l, o in re.findall(r'(\d+)([MIS=X])', processed_cigar))
        
        # Calculate length of groups that come after this one (for end hard clip)
        for idx in range(group_idx + 1, len(processed_groups)):
            processed_cigar = ''.join(processed_groups[idx]['cigar_parts'])
            # Calculate query-consuming operations from processed CIGAR
            following_groups_length += sum(int(l) for l, o in re.findall(r'(\d+)([MIS=X])', processed_cigar))
        
        # Create record for this group
        record = primary_read.copy()
        record[0] = base_read_name
        record[1] = '0' if group_idx == 0 else '2048'  # Primary vs supplementary
        record[2] = processed_group['ref_name']
        record[3] = str(processed_group['ref_pos'])
        record[4] = '255'
        
        # Build CIGAR with directional hard clips
        cigar_parts = []
        
        # Add beginning hard clip if there are preceding groups
        if preceding_groups_length > 0:
            cigar_parts.append(f"{preceding_groups_length}H")
        
        # Add the main CIGAR for this group
        cigar_parts.extend(processed_group['cigar_parts'])
        
        # Add end hard clip if there are following groups
        if following_groups_length > 0:
            cigar_parts.append(f"{following_groups_length}H")
        
        record[5] = ''.join(cigar_parts)
        
        # Set mate information for paired records
        if len(processed_groups) > 1:
            if group_idx == 0:
                # Primary points to first supplementary
                next_group = processed_groups[1]
                record[6] = next_group['ref_name']
                record[7] = str(next_group['ref_pos'])
                record[8] = '0'
            else:
                record[6] = '*'
                record[7] = '0'
                record[8] = '0'
        
        # Add tags
        fi_value = group_idx + 1
        tc_value = len(processed_groups)
        tags = [f'FI:i:{fi_value}', 'NM:i:0', f'TC:i:{tc_value}']
        
        records.append(record[:11] + tags)
    
    return records

def get_mapped_softclip_positions(base_read_name, sam_file='output.sam'):
    """Parse intermediate SAM file to determine which nested softclip read names are actually mapped"""
    mapped_read_names = set()
    
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            read_name = fields[0]
            flag = int(fields[1])
            
            # Check if this is a nested softclip read for our base read
            if read_name.startswith(base_read_name + '_softclip_'):
                # Check if it's mapped (flag bit 4 not set)
                if not (flag & 4):
                    # Store the full mapped read name
                    mapped_read_names.add(read_name)
    
    return mapped_read_names


# Global cache for Python implementation
_softclip_cache = None

def get_mapped_softclip_positions_cached(base_read_name, sam_file='output.sam'):
    """Python cached version of get_mapped_softclip_positions - reads file once and caches all results"""
    global _softclip_cache
    
    # Check if cache exists
    if _softclip_cache is not None:
        return _softclip_cache.get(base_read_name, set())
    
    # Build cache - read file once
    print("Building softclip cache from SAM file (Python)...")
    _softclip_cache = {}
    
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
                
            read_name = fields[0]
            flag = int(fields[1])
            
            # Check if this is a softclip read
            if '_softclip_' in read_name:
                # Extract base read name
                base_name = read_name.split('_softclip_')[0]
                # Check if it's mapped (flag bit 4 not set)
                if not (flag & 4):
                    if base_name not in _softclip_cache:
                        _softclip_cache[base_name] = set()
                    _softclip_cache[base_name].add(read_name)
    
    print(f"Softclip cache built with {len(_softclip_cache)} base reads (Python)")
    return _softclip_cache.get(base_read_name, set())

def clear_softclip_cache():
    """Clear the Python softclip cache"""
    global _softclip_cache
    _softclip_cache = None
    print("Python softclip cache cleared")

def update_hardclipped_values(read):
    """
    Remove hardclipped nucleotides from sequence and corresponding base quality scores for reads with hardclips.
    Modifies read[9] (sequence) and read[10] (quality scores) to contain only the portions corresponding to the non-hardclipped regions.
    """
    import re
    
    sequence = read[9]
    quality_scores = read[10]
    cigar = read[5]
    
    # Check if there are hardclips in the CIGAR string
    if 'H' not in cigar:
        return  # No hardclips, nothing to do
    
    # Parse CIGAR operations
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    
    # Calculate how many nucleotides to remove from the start and end
    start_hardclip = 0
    end_hardclip = 0
    
    # Check for hardclip at the beginning
    if cigar_ops and cigar_ops[0][1] == 'H':
        start_hardclip = int(cigar_ops[0][0])
    
    # Check for hardclip at the end
    if cigar_ops and cigar_ops[-1][1] == 'H':
        end_hardclip = int(cigar_ops[-1][0])
    
    # Remove hardclipped regions from both sequence and quality scores
    if start_hardclip > 0 or end_hardclip > 0:
        if end_hardclip > 0:
            # Remove from both ends
            updated_sequence = sequence[start_hardclip:-end_hardclip]
            updated_quality = quality_scores[start_hardclip:-end_hardclip]
        else:
            # Remove only from the start
            updated_sequence = sequence[start_hardclip:]
            updated_quality = quality_scores[start_hardclip:]
        
        # Update both sequence and quality scores in the read
        read[9] = updated_sequence
        read[10] = updated_quality

def build_hierarchical_accounting_map(base_read_name, segments, sam_file='output.sam'):
    """Build accounting map that considers the hierarchical softclip structure"""
    # Get all mapped softclip read names from SAM file - use cached version
    mapped_read_names = get_mapped_softclip_positions_cached(base_read_name, sam_file)
    
    # Build accounting map for each segment based on which nested softclips are mapped
    segment_accounting = {}
    
    for segment in segments:
        segment_read_name = segment['read'][0]  # Get the read name from the segment
        segment_accounting[segment_read_name] = set()
        
        if segment['type'] == 'main':
            # For main alignment, check if direct softclip reads are mapped
            # If base_softclip_0 or base_softclip_1 are mapped, account for those positions
            direct_softclip_0 = f"{base_read_name}_softclip_0"
            direct_softclip_1 = f"{base_read_name}_softclip_1"
            
            if direct_softclip_0 in mapped_read_names:
                segment_accounting[segment_read_name].add(0.4)  # softclip_0 position
            if direct_softclip_1 in mapped_read_names:
                segment_accounting[segment_read_name].add(0.6)  # softclip_1 position
        
        elif segment['type'] == 'softclip':
            # For softclip segments, check if nested softclips from this segment are mapped
            # If they are, then we need to account for the corresponding positions
            for mapped_read_name in mapped_read_names:
                # Check if this mapped read is a nested softclip from the current segment
                if mapped_read_name.startswith(segment_read_name + '_softclip_'):
                    # Extract what position this nested softclip represents in the current segment
                    nested_suffix = mapped_read_name[len(segment_read_name):]
                    if nested_suffix.startswith('_softclip_0'):
                        # This is a softclip_0 from the current segment, so position 0 is accounted for
                        segment_accounting[segment_read_name].add(0.4)
                    elif nested_suffix.startswith('_softclip_1'):
                        # This is a softclip_1 from the current segment, so position 1 is accounted for  
                        segment_accounting[segment_read_name].add(0.6)
    
    return segment_accounting


def _true_query_start(read_fields):
    """Start offset of this record's aligned portion within the ORIGINAL query orientation
    (i.e. as it appeared before any strand flip). For a reverse-strand record (flag & 16), its
    own CIGAR/SEQ are expressed relative to the revcomp of the original query, so its own
    leading clip corresponds to the original query's TRAILING side, and vice versa -- this
    flips the two before comparing, so two segments from the same identifier (one forward, one
    reverse) sort into their true left-to-right order along the original read.
    """
    flag = int(read_fields[1])
    ops = re.findall(r'(\d+)([MIDNSHPX=])', read_fields[5])
    lead = 0
    for length, op in ops:
        if op in 'HS':
            lead += int(length)
        else:
            break
    trail = 0
    for length, op in reversed(ops):
        if op in 'HS':
            trail += int(length)
        else:
            break
    return trail if (flag & 16) else lead


def _own_frame_side_for_true_side(flag, true_side):
    """Map a side ('leading' or 'trailing') in the ORIGINAL query orientation to which side
    that is in this record's OWN CIGAR/SEQ frame, accounting for strand (see
    _true_query_start)."""
    if not (flag & 16):
        return true_side
    return 'trailing' if true_side == 'leading' else 'leading'


def _clip_side_to_hardclip(read_fields, own_frame_side):
    """Convert this record's own leading or trailing S (if present) to H. Used when two
    segments from the same softclip identifier split into primary+supplementary: the side of
    each that faces its sibling is already explained by that sibling, so it must become a hard
    clip -- otherwise extract_aligned_operations (keyed on read name, which both segments
    share) can't tell it's already accounted for, and would re-add it as a duplicate
    operation."""
    ops = re.findall(r'(\d+)([MIDNSHPX=])', read_fields[5])
    if not ops:
        return read_fields
    idx = 0 if own_frame_side == 'leading' else -1
    if ops[idx][1] != 'S':
        return read_fields
    new_ops = list(ops)
    new_ops[idx] = (new_ops[idx][0], 'H')
    new_read = list(read_fields)
    new_read[5] = ''.join(f'{length}{op}' for length, op in new_ops)
    return new_read


def _ref_span(cigar):
    """Reference bases consumed by a CIGAR (M/D/N/=/X)."""
    ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    return sum(int(length) for length, op in ops if op in 'MDN=X')


def _aligned_query_len(cigar):
    """Query bases actually aligned by a CIGAR (M/I/=/X) -- used as an alignment-quality proxy
    when choosing which of two overlapping segments to keep."""
    ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    return sum(int(length) for length, op in ops if op in 'MI=X')


def _is_ins_split_sibling(name_a, name_b):
    """True if name_a/name_b are the two deliberately-complementary halves of one
    _softclip_candidates_for_record split (isolated insertion + the rest that follows it) --
    they're EXPECTED to overlap heavily in reference space (the insertion's true origin can sit
    anywhere the flank also naturally covers) and must never be treated as duplicates of each
    other, unlike a genuine coincidental overlap between two unrelated segments."""
    return name_a == name_b + '_ins' or name_b == name_a + '_ins'


def _dedupe_redundant_segments(segments):
    """Drop softclip segments that substantially overlap another segment already kept, on the
    same reference -- keeping whichever has more aligned query bases. Only fires on genuinely
    high reciprocal overlap (>50% of the smaller segment's reference span); segments that only
    partially/coincidentally overlap are left untouched, since in practice ViReMa's own
    independently-rescued nested softclips sometimes carry real, different information from a
    same-locus minimap2 supplementary rather than being a pure duplicate of it. Never fires
    between the two halves of an insertion/rest split pair (see _is_ins_split_sibling) -- those
    are deliberately complementary, not competing alternatives."""
    kept = []
    for seg in segments:
        if seg['type'] != 'softclip':
            kept.append(seg)
            continue
        seg_span = _ref_span(seg['cigar'])
        seg_start, seg_end = seg['ref_pos'], seg['ref_pos'] + seg_span
        duplicate_of = None
        for other in kept:
            if other['type'] != 'softclip' or other['ref_name'] != seg['ref_name']:
                continue
            if _is_ins_split_sibling(seg['read'][0], other['read'][0]):
                continue
            other_span = _ref_span(other['cigar'])
            other_start, other_end = other['ref_pos'], other['ref_pos'] + other_span
            overlap = min(seg_end, other_end) - max(seg_start, other_start)
            if seg_span and other_span and overlap > 0 and overlap / min(seg_span, other_span) > 0.5:
                duplicate_of = other
                break
        if duplicate_of is None:
            kept.append(seg)
        elif _aligned_query_len(seg['cigar']) > _aligned_query_len(duplicate_of['cigar']):
            kept.remove(duplicate_of)
            kept.append(seg)
        # else: seg is the worse duplicate of an already-kept segment -- drop it
    return kept


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
    softclip_reads = {}  # sequence_position -> list of reads at that position (usually 1)

    for read in reads:
        read_name = read[0]
        flag = int(read[1])

        if '_softclip_' in read_name:
            # Parse softclip path to get sequence position
            sequence_position = get_sequence_position_for_read(read_name, base_read_name, None)
            if sequence_position is not None and flag & 4 == 0:  # Only if aligned
                softclip_reads.setdefault(sequence_position, []).append(read)
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
    main_seq_pos = 0.5  # Between softclip_0 and softclip_1
    segments.append({
        'ref_pos': main_ref_pos,
        'ref_name': main_ref_name,
        'seq_pos': main_seq_pos,
        'type': 'main',
        'read': primary_read,
        'cigar': primary_read[5]
    })
    
    # Add softclip segments. Usually exactly one read per position, but when
    # rewrite_embedded_insertions_as_softclips exposes a genuine two-locus junction, this
    # identifier's own realignment can itself split into primary+supplementary -- both must be
    # kept as their own segments (not collapsed to one), correctly ordered and with the query
    # territory each of them claims made non-overlapping.
    for seq_pos, reads_at_pos in softclip_reads.items():
        if len(reads_at_pos) == 1:
            r = reads_at_pos[0]
            segments.append({
                'ref_pos': int(r[3]), 'ref_name': r[2], 'seq_pos': seq_pos,
                'type': 'softclip', 'read': r, 'cigar': r[5]
            })
        elif len(reads_at_pos) == 2:
            ins_reads = [r for r in reads_at_pos if r[0].endswith('_ins')]
            if len(ins_reads) == 1:
                # One is the isolated-insertion half of a split candidate pair (see
                # _softclip_candidates_for_record) -- it always precedes its "rest" sibling in
                # the original fragment BY CONSTRUCTION, not by inference. _true_query_start
                # infers order from each record's own soft-clip amount, which only reflects
                # position in a SHARED parent frame when both records came from splitting one
                # minimap2 alignment (the natural primary+supplementary case); here the two were
                # aligned independently as disjoint substrings, so their own clip amounts don't
                # mean that and can't be trusted to order them.
                r0 = ins_reads[0]
                r1 = [r for r in reads_at_pos if r is not r0][0]
            else:
                r0, r1 = sorted(reads_at_pos, key=_true_query_start)
            r0 = _clip_side_to_hardclip(r0, _own_frame_side_for_true_side(int(r0[1]), 'trailing'))
            r1 = _clip_side_to_hardclip(r1, _own_frame_side_for_true_side(int(r1[1]), 'leading'))
            for sub_idx, r in enumerate((r0, r1)):
                segments.append({
                    'ref_pos': int(r[3]), 'ref_name': r[2],
                    'seq_pos': seq_pos + (sub_idx - 0.5) * 1e-6,
                    'type': 'softclip', 'read': r, 'cigar': r[5]
                })
        else:
            # 3+ candidates at one identifier -- fall back to the original single-best-AS
            # behavior rather than guess at an N-way ordering we haven't observed/verified.
            best = max(reads_at_pos, key=get_as_score_from_read_fields)
            segments.append({
                'ref_pos': int(best[3]), 'ref_name': best[2], 'seq_pos': seq_pos,
                'type': 'softclip', 'read': best, 'cigar': best[5]
            })
    
    # Step 2b: Drop softclip segments that are near-duplicates of another segment already
    # present -- can happen when a nested softclip-of-a-softclip independently rediscovers
    # essentially the same locus that its own parent's supplementary alignment already found
    # (they come from different softclip identifiers/positions, so nothing earlier catches
    # this). Only applied to genuinely high-overlap pairs; low/partial overlap is left alone.
    segments = _dedupe_redundant_segments(segments)

    # Step 3: Sequence ordering - sort by sequence position
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



