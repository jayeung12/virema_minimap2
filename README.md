# Minimap2 Multi-Round Alignment with ViReMa Format Conversion

ViReMa but with minimap2. Multi-round alignment approach with Minimap2, followed by conversion to ViReMa-compatible format. Everything downstream of SAM compilation is mostly the same.

## Overview

1. Performs initial alignment with Minimap2
2. Identifies and re-aligns softclipped regions to detect split alignments
3. Converts alignment patterns to ViReMa-compatible format
4. Preserves full CIGAR complexity (I and Ds) while representing recombination events

## Minimap2_Module

**Usage:**
```bash
python ViReMa.py --Aligner minimap2  --Seed 25 ./Test_Data/FHV_Genome.txt ./Test_Data/FHV_small.txt  output_mm2.sam --MicroInDel_Length 20

python ViReMa.py --Aligner minimap2 -lr ont  --Seed 25 Test_Data/SARS2_Genome.fasta Test_Data/combined_1000bp_duplications.fastq output_mm2_ont_dups.sam --MicroInDel_Length 20
```

**New parameter:**
- `-lr`: Long read technology (ont for Oxford Nanopore, pb for PacBio CLR, hifi for PacBio HiFi)

## Logic

### 1. Initial Alignment
- **Tool**: Minimap2 with technology-specific parameters:
  - **Short reads (default)**: `-ax sr -k 20 -A 1 -B 2`
  - **Oxford Nanopore**: `-ax map-ont`
  - **PacBio CLR**: `-ax map-pb` (not implemented yet)
  - **PacBio HiFi**: `-ax map-hifi` (not implemented yet)
- **Purpose**: Broad alignment to identify primary mappings and softclipped regions
- **Output**: Primary and supplemental alignments with softclipped portions

Short reads that map often have supplemental alignments. ONT reads occassionally have supplemental alignments or multiple primary alignments in different orientations. Only the alignment with the highest alignment score has its softclips sent for additional rounds of alignment.

### 2. Softclip Extraction and Re-alignment
- **Extraction**: Identifies softclipped sequences ≥ threshold length (Seed parameter) from primary alignments
- **Naming Convention**: 
  - `softclip_0`: Softclips occurring **before** the main alignment
  - `softclip_1`: Softclips occurring **after** the main alignment
- **Re-alignment**: Uses technology-specific parameters:
  - **Short reads**: More sensitive parameters (`-k 10 -w 5 -m 10`) to map softclipped regions
  - **Long reads**: Same technology-specific presets as initial alignment
- **Purpose**: Detect split alignments indicating potential recombination events

Sometimes, mapping softclips will have further softclips that exceed threshold length. These are sent for another round of mapping where the read name has an additional softclip_0 or softclip_1 appended to the read name in the intermediate file (output.sam). If these softclips map, they are stitched back to the primary alignment in the correct order. Softclips that remain unmapped but are internal (not at the ends of the final, stitched together mapping) are turned into I events in the CIGAR string.

### 3. Result Merging and Classification

#### Single Alignments
- Reads with only primary alignment
- **Action**: Preserve original softclips as unmapped regions
- **Output**: Standard SAM record with softclips intact

#### Multiple Alignments
Reads with primary + supplemental + softclip alignments are processed based on genomic positioning:

##### In-Order Segments (Recombination Events)
- Softclip mappings occur in expected genomic order relative to primary alignment
- **Action**: Create single merged record with N gaps representing genomic distances
- **Example**: 
  ```
  Original: 50M40S (primary) + 40M (softclip at distant location)
  Converted: 50M150N40M (gap of 150bp between segments)
  ```

##### Out-of-Order Segments Or Inter-Segment Events
- Softclip mappings occur before the primary alignment genomically
- **Action**: Create paired records with hard clips (ViReMa-style)
- **Example**:
  ```
  Record 1: 59M31H (primary + hard clip for out-of-order softclip)
  Record 2: 59H31M (hard clip for primary + mapped softclip)
  Tags: FI:i:1, FI:i:2, TC:i:2
  ```

## Output Format Specifications

### Single Merged Records
```
QNAME  FLAG  RNAME  POS  MAPQ  CIGAR           RNEXT  PNEXT  TLEN  SEQ  QUAL  TAGS
read1  0     ref    100  255   50M150N40M7S    *      0      0     ...  ...   NM:i:0
```

### Paired Records (Out-of-order)
```
read1  0     ref    100  255   59M31H          ref    50     0     ...  ...   FI:i:1 NM:i:0 TC:i:2
read1  2048  ref    50   255   59H31M          *      0      0     ...  ...   FI:i:2 NM:i:0 TC:i:2
```

## Intermediate Files

- **`output.sam`**: Complete alignment results (initial + merged softclip alignments)
- **`multiRound`**: Filtered reads with qualifying softclipped regions 
- **`TEMP_SAM`**: Secondary alignment results for extracted softclips
- **`Test_Data/TEMP_READS.txt`**: FASTA format extracted softclipped sequences

## System Requirements

- **Python**: 3.x with standard libraries (`subprocess`, `re`, `argparse`, `collections`)
- **Minimap2**: Must be available in system PATH
- **Memory**: Proportional to input size (typically modest requirements)
- **Storage**: ~3x input file size for intermediate files
