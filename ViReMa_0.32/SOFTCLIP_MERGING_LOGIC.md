# Softclip Merging Logic for ViReMa Minimap2 Module

## Overview

This document describes the detailed logic for merging softclip alignments with their corresponding main alignments in the ViReMa Minimap2 module. The merging process must respect both sequence ordering (based on softclip paths) and genomic positioning rules to determine whether segments can be merged with N gaps or must be split into paired records.

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

### 1. Path-Based Positioning
Each softclip read gets a position tuple based on its path:
- `read_softclip_0` → `(0,)`
- `read_softclip_1` → `(1,)`
- `read_softclip_0_softclip_1` → `(0, 1)`
- `read_softclip_1_softclip_0` → `(1, 0)`
- `read_softclip_0_softclip_1_softclip_0` → `(0, 1, 0)`
- Main alignment → `(0.5,)` (conceptually between 0 and 1)

### 2. Sequence Order Determination
Sort all segments by their path tuples to get the correct sequence order:
1. All `(0,...)` paths (most nested first)
2. Main alignment `(0.5,)`
3. All `(1,...)` paths (most nested first)

### Example Sequence Orders

**Case 1: Simple nesting**
- Segments: `softclip_0`, `softclip_0_softclip_1`, main, `softclip_1`
- Sequence order: `softclip_0` → `softclip_0_softclip_1` → main → `softclip_1`

**Case 2: Complex nesting**
- Segments: `softclip_0`, `softclip_0_softclip_1`, `softclip_0_softclip_1_softclip_0`, main, `softclip_1_softclip_0`, `softclip_1`
- Sequence order: `softclip_0` → `softclip_0_softclip_1` → `softclip_0_softclip_1_softclip_0` → main → `softclip_1_softclip_0` → `softclip_1`

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
- `SRR5085928.4191_softclip_0`: position 1000, CIGAR `200M100S`, end = 1299
- `SRR5085928.4191_softclip_0_softclip_1`: position 500, CIGAR `100M`, end = 599
- Main alignment: position 1200, CIGAR `300S200M`, end = 1399

**Sequence order:** softclip_0 → softclip_0_softclip_1 → main

**Genomic analysis:**
- softclip_0 (1299) vs softclip_0_softclip_1 (500): Out of order (1299 >= 500) ❌
- softclip_0_softclip_1 (599) vs main (1200): 599 < 1200 ✓ Can merge

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
1. Parse each read name to extract softclip path
2. Convert path to position tuple
3. Assign sequence positions to all segments

### Step 2: Sequence Ordering
1. Sort all segments by position tuples
2. This gives the correct sequence order regardless of genomic positions

### Step 3: Genomic Merging
1. For each adjacent pair in sequence order:
   - Calculate genomic spans (start to end positions)
   - Check for same reference and no overlap
   - Group contiguous segments for N-gap merging

### Step 4: Output Generation
1. If all segments can be merged: create single merged record
2. If segments must be split: create paired records with appropriate hard clips

### Step 5: AS Score Selection
- When multiple alignments map to the same sequence position, select the one with the highest AS (alignment score)
- This handles supplementary alignments for the same softclip

## Key Considerations

1. **Sequence order takes precedence** for determining the logical flow of the read
2. **Genomic merging rules determine** the technical feasibility of merging
3. **N-gap merging happens first**, then paired record creation
4. **Hard clips in paired records** represent the portions of the sequence present in other records
5. **AS scores resolve conflicts** when multiple alignments exist for the same sequence position

This logic ensures that both the biological sequence relationships and the technical genomic constraints are properly respected in the final output.
