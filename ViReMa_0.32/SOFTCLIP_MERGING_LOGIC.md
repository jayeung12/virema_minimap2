# Softclip Merging Logic for ViReMa Minimap2 Module

## Overview

This document describes the detailed logic for merging softclip alignments with their corresponding main alignments in the ViReMa Minimap2 module. The merging process uses a numerical position-based system for sequence ordering and genomic positioning rules to determine whether segments can be merged with N gaps or must be split into paired records.

## Softclip Naming Convention

### Basic Convention
- `read_softclip_0`: Softclip derived from the **beginning** of an alignment (before main alignment)
- `read_softclip_1`: Softclip derived from the **end** of an alignment (after main alignment)

### Nested Convention
- `read_softclip_0_softclip_1`: Softclip derived from the **end** of the `read_softclip_0` alignment
- `read_softclip_1_softclip_0`: Softclip derived from the **beginning** of the `read_softclip_1` alignment
- `read_softclip_0_softclip_1_softclip_0`: Further nesting (beginning of the `read_softclip_0_softclip_1` alignment)

### Position Interpretation
- All `softclip_0` paths (regardless of nesting) are **before** the main alignment in sequence order
- All `softclip_1` paths (regardless of nesting) are **after** the main alignment in sequence order
- Nested paths indicate the actual contiguity relationships

## Sequence Ordering Rules

### 1. Numerical Position-Based System
Each softclip read gets a calculated numerical position based on its path:
- Main alignment → `0.5` (reference point)
- `read_softclip_0` → `0.4` (0.5 - 0.1)
- `read_softclip_1` → `0.6` (0.5 + 0.1)
- `read_softclip_0_softclip_1` → `0.45` (0.4 + 0.05)
- `read_softclip_1_softclip_0` → `0.55` (0.6 - 0.05)
- `read_softclip_0_softclip_1_softclip_0` → `0.425` (0.45 - 0.025)

### 2. Position Calculation Algorithm
```python
def calculate_sequence_position(read_name, base_read_name):
    # Start from main alignment position
    position = 0.5
    base_offset = 0.1
    
    for level, direction in enumerate(path_indices):
        # Calculate offset for this level (smaller for deeper nesting)
        level_offset = base_offset / (2 ** level)
        
        if direction == 0:  # softclip_0: move earlier
            position -= level_offset
        elif direction == 1:  # softclip_1: move later
            position += level_offset
    
    return position
```

### 3. Sequence Order Determination
Sort all segments by their numerical positions to get the correct sequence order:
1. All segments with position < 0.5 (softclip_0 variants)
2. Main alignment at position 0.5
3. All segments with position > 0.5 (softclip_1 variants)

### Example Sequence Orders

**Case 1: Simple nesting**
- `softclip_0` (0.4) → `softclip_0_softclip_1` (0.45) → main (0.5) → `softclip_1` (0.6)

**Case 2: Complex nesting**
- `softclip_0` (0.4) → `softclip_0_softclip_1` (0.45) → `softclip_0_softclip_1_softclip_0` (0.425) → main (0.5) → `softclip_1_softclip_0` (0.55) → `softclip_1` (0.6)

Note: The numerical system automatically handles proper nesting order through the decreasing offset calculation.

## Genomic Merging Rules

### 1. Overlap Detection
For any two adjacent segments in sequence order, calculate:
- **Start position**: From SAM field 4 (POS)
- **End position**: Start + sum of reference-consuming CIGAR operations (M/D/N/=X) - 1
- **Overlap check**: If `end_pos_segment1 >= start_pos_segment2`, segments overlap

### 2. Merging Criteria
Two adjacent segments can be merged with an N gap if:
1. They map to the **same reference sequence** (RNAME)
2. They are **genomically contiguous** (no overlap): `end_pos_segment1 < start_pos_segment2`

### 3. Paired Records Criteria
Segments must be split into paired records if:
1. They map to **different reference sequences**, OR
2. They **overlap genomically**, OR
3. They are **out of genomic order**

## Detailed Examples

### Example 1: Simple Contiguous Case

**Input segments (after sequence ordering):**
- `SRR123_softclip_0`: position 100, CIGAR `50M`, end = 149
- Main alignment: position 200, CIGAR `50S100M30S`, end = 299
- `SRR123_softclip_1`: position 350, CIGAR `30M`, end = 379

**Genomic analysis:**
- softclip_0 (149) < main (200): ✓ Can merge with N gap
- main (299) < softclip_1 (350): ✓ Can merge with N gap

**Output:** Single merged record
- Position: 100
- CIGAR: `50M50N100M50N30M` (with appropriate N gaps)

### Example 2: Overlap Case (Paired Records)

**Input segments (after sequence ordering):**
- `SRR123_softclip_0`: position 200, CIGAR `200M`, end = 399
- Main alignment: position 300, CIGAR `200S100M`, end = 399

**Genomic analysis:**
- softclip_0 ends at 399, main starts at 300
- 399 >= 300: ❌ Overlap detected

**Output:** Paired records
- **Record 1** (Primary): softclip_0 only + hard clip for main
- **Record 2** (Supplementary): hard clip for softclip_0 + main alignment

### Example 3: Complex Nested Case (From Original Issue)

**Input segments:**
- `SRR5085928.4191_softclip_0`: genomic position 1000, CIGAR `200M100S`, end = 1299, sequence position = 0.4
- `SRR5085928.4191_softclip_0_softclip_1`: genomic position 500, CIGAR `100M`, end = 599, sequence position = 0.45
- Main alignment: genomic position 1200, CIGAR `300S200M`, end = 1399, sequence position = 0.5

**Sequence order (by numerical position):** softclip_0 (0.4) → softclip_0_softclip_1 (0.45) → main (0.5)

**Genomic analysis:**
- softclip_0 (genomic end 1299) vs softclip_0_softclip_1 (genomic start 500): Out of order (1299 >= 500) ❌
- softclip_0_softclip_1 (genomic end 599) vs main (genomic start 1200): 599 < 1200 ✓ Can merge

**First merge contiguous segments:**
- Group 1: `softclip_0` (alone)
- Group 2: `softclip_0_softclip_1` + main (merged with N gap)

**Output:** Paired records
- **Record 1** (Primary): 
  - Position: 1000
  - CIGAR: `200M[hardclip for group 2]H`
- **Record 2** (Supplementary):
  - Position: 500  
  - CIGAR: `[hardclip for softclip_0]H100M601N200M`

### Example 4: Different Reference Sequences

**Input segments:**
- `SRR123_softclip_0`: chr1:1000, CIGAR `100M`
- Main alignment: chr1:1200, CIGAR `100S150M75S`  
- `SRR123_softclip_1`: chr2:500, CIGAR `75M`

**Analysis:**
- softclip_0 and main: same reference, contiguous ✓
- main and softclip_1: different references ❌

**First merge contiguous segments:**
- Group 1: softclip_0 + main (merged)
- Group 2: softclip_1 (alone)

**Output:** Paired records
- **Record 1**: chr1:1000 with merged softclip_0 + main
- **Record 2**: chr2:500 with softclip_1

## Implementation Algorithm

### Step 1: Parse and Position
1. For each read in the read group:
   - If it's a softclip read (`_softclip_` in name): calculate numerical sequence position
   - If it's the main alignment: assign position 0.5
   - Select best alignment by AS score if multiple alignments exist for same position

### Step 2: Sequence Ordering
1. Sort all segments by their numerical positions (float values)
2. This gives the correct sequence order regardless of genomic positions

### Step 3: Genomic Merging
1. For each adjacent pair in sequence order:
   - Calculate genomic spans (start to end positions)
   - Check for same reference and no overlap: `end_pos_segment1 < start_pos_segment2`
   - Group contiguous segments for N-gap merging

### Step 4: Output Generation
1. If all segments can be merged: create single merged record with N gaps
2. If segments must be split: create paired records with appropriate hard clips

### Step 5: AS Score Selection
- When multiple alignments map to the same sequence position, select the one with the highest AS (alignment score)
- This handles supplementary alignments for the same softclip automatically during read collection

## Key Considerations

1. **Numerical sequence positions** determine the logical flow of the read through calculated float values
2. **Genomic merging rules determine** the technical feasibility of merging (same reference + no overlap)
3. **AS scores resolve conflicts** when multiple alignments exist for the same sequence position
4. **Hierarchical accounting** ensures proper softclip handling in merged records through SAM-based analysis
5. **N-gap merging happens first**, then paired record creation if needed
6. **Hard clips in paired records** represent the portions of the sequence present in other records

This logic ensures that both the biological sequence relationships and the technical genomic constraints are properly respected in the final output, using a robust numerical positioning system that handles arbitrary nesting complexity.
