# Minimap2 Multi-Round Alignment with ViReMa Format Conversion

This repository contains a comprehensive pipeline for detecting viral recombination events using a multi-round alignment approach with Minimap2, followed by conversion to ViReMa-compatible format.

## Overview

The project implements a sophisticated alignment strategy that:
1. Performs initial sensitive alignment with Minimap2
2. Identifies and re-aligns softclipped regions to detect split alignments
3. Converts complex alignment patterns to ViReMa-compatible format
4. Preserves full CIGAR complexity while representing recombination events

## Key Features

- **Multi-round alignment**: Initial broad alignment followed by targeted softclip re-alignment
- **Intelligent format conversion**: Handles complex CIGAR operations, insertions, deletions, and gaps
- **Recombination detection**: Identifies in-order and out-of-order alignment patterns
- **ViReMa compatibility**: Produces output compatible with viral recombination analysis tools

## Main Script: `multiround_alignment.py`

### Core Functionality
A complete pipeline that processes reads through multiple alignment rounds and converts the results to ViReMa format.

**Usage:**
```bash
python multiround_alignment.py <virus_index> <input_fastq> <output_sam> [--Seed <threshold>] [-lr <technology>]
```

**Parameters:**
- `virus_index`: Path to virus genome reference file
- `input_fastq`: Input FASTQ file with reads
- `output_sam`: Output SAM file in ViReMa format
- `--Seed`: Minimum softclip length threshold (default: 15bp)
- `-lr`: Long read technology (ont for Oxford Nanopore, pb for PacBio CLR, hifi for PacBio HiFi)

**Examples:**
```bash
# Short read analysis (default)
python multiround_alignment.py ./Test_Data/FHV_Genome.txt ./Test_Data/FHV_small.txt ./output_small.sam --Seed 20

# Oxford Nanopore data
python multiround_alignment.py ./Test_Data/FHV_Genome.txt ./Test_Data/FHV_nanopore.fastq ./output_ont.sam -lr ont

# PacBio CLR data
python multiround_alignment.py ./Test_Data/FHV_Genome.txt ./Test_Data/FHV_pacbio.fastq ./output_pb.sam -lr pb

# PacBio HiFi data
python multiround_alignment.py ./Test_Data/FHV_Genome.txt ./Test_Data/FHV_hifi.fastq ./output_hifi.sam -lr hifi
```

## Pipeline Workflow

### 1. Initial Alignment
- **Tool**: Minimap2 with technology-specific parameters:
  - **Short reads (default)**: `-ax sr -k 20 -A 1 -B 2` (viral-optimized)
  - **Oxford Nanopore**: `-ax map-ont`
  - **PacBio CLR**: `-ax map-pb`
  - **PacBio HiFi**: `-ax map-hifi`
- **Purpose**: Broad alignment to identify primary mappings and softclipped regions
- **Output**: Primary and supplemental alignments with softclipped portions

### 2. Softclip Extraction and Re-alignment
- **Extraction**: Identifies softclipped sequences ≥ threshold length from primary alignments
- **Naming Convention**: 
  - `softclip_0`: Softclips occurring **before** the main alignment
  - `softclip_1`: Softclips occurring **after** the main alignment
- **Re-alignment**: Uses technology-specific parameters:
  - **Short reads**: More sensitive parameters (`-k 10 -w 5 -m 10`) to map softclipped regions
  - **Long reads**: Same technology-specific presets as initial alignment
- **Purpose**: Detect split alignments indicating potential recombination events

### 3. Result Merging and Classification

The pipeline classifies reads into several categories:

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

##### Out-of-Order Segments (Complex Events)
- Softclip mappings occur before the primary alignment genomically
- **Action**: Create paired records with hard clips (ViReMa-style)
- **Example**:
  ```
  Record 1: 59M31H (primary + hard clip for out-of-order softclip)
  Record 2: 59H31M (hard clip for primary + mapped softclip)
  Tags: FI:i:1, FI:i:2, TC:i:2
  ```

### 4. Advanced CIGAR Preservation

#### Complex Primary Alignments
- **Preserves**: Full CIGAR complexity including deletions (D), insertions (I), matches (M)
- **Example**: `40S18M9D32M` → `40M1017N18M9D32M`
- **Maintains**: Complete alignment information while adding recombination gaps

#### Softclip Mapping Complexity
- **Preserves**: Insertions and trailing softclips from softclip mappings
- **Example**: Softclip CIGAR `35M7S` → Final: `48M798N35M7S`
- **Purpose**: Retains unmapped portions that couldn't be aligned

#### CIGAR Operation Handling
- **M, =, X**: Converted to M (matches)
- **I**: Preserved as insertions
- **D**: Preserved as deletions  
- **N**: Added for recombination gaps
- **S**: Preserved for unmapped regions
- **H**: Used in paired records for hard clips

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

### Key Fields
- **FLAG**: 0 (primary), 2048 (supplementary)
- **MAPQ**: 255 for all mapped reads
- **CIGAR**: Complex preservation with N gaps for recombination
- **Tags**: 
  - `NM:i:X`: Edit distance
  - `FI:i:X`: Fragment index (1/2 for paired records)
  - `TC:i:X`: Total count of records for this read

## Sample Data and Testing

### Test Dataset
- **Reference**: `Test_Data/FHV_Genome.txt` (Flock House Virus)
- **Reads**: `Test_Data/FHV_small.txt` (Sample FASTQ)
- **Expected Output**: Various recombination patterns and complex alignments

### Validation Examples
- **Standard recombination**: `50M150N40M` (in-order segments)
- **Complex recombination**: `48M798N35M7S` (with trailing softclip)
- **Out-of-order events**: Paired records with hard clips
- **Complex CIGARs**: `48M1D7M569N35M` (deletion preservation)

## Technical Implementation

### Decision Algorithm
```
For reads with primary + supplemental + softclip:
├── Check genomic order of softclip mappings
│   ├── softclip_1 maps before primary? → Out-of-order
│   └── softclip_0 maps after primary? → Out-of-order
├── Out-of-order? → Create paired records with hard clips
└── In-order? → Create single merged record with N gaps
```

### Softclip Index Logic
- **Determination**: Based on position relative to main alignment in read sequence
- **softclip_0**: Occurs before any M/I/=/X operations in CIGAR
- **softclip_1**: Occurs after any M/I/=/X operations in CIGAR
- **Purpose**: Maintains consistent mapping of sequence portions to genomic locations

### Quality Assurance
- **Sequence validation**: Ensures CIGAR operations match sequence length
- **Reference consumption**: Validates gap calculations and positioning
- **Format compliance**: Maintains ViReMa compatibility requirements

## Intermediate Files

- **`output.sam`**: Complete alignment results (initial + merged softclip alignments)
- **`multiRound`**: Filtered reads with qualifying softclipped regions  
- **`TEMP_SAM`**: Secondary alignment results for extracted softclips
- **`Test_Data/TEMP_READS.txt`**: FASTA format extracted softclipped sequences
- **Final output**: User-specified ViReMa-compatible SAM file

## System Requirements

- **Python**: 3.x with standard libraries (`subprocess`, `re`, `argparse`, `collections`)
- **Minimap2**: Must be available in system PATH
- **Memory**: Proportional to input size (typically modest requirements)
- **Storage**: ~3x input file size for intermediate files

## Applications

This pipeline is designed for:
- **Viral recombination analysis**: Detection of recombination junctions in viral genomes
- **Complex alignment scenarios**: Handling reads spanning genomic rearrangements  
- **ViReMa workflow integration**: Preprocessing for downstream recombination analysis
- **Research applications**: Studying viral evolution and recombination patterns

## Notes and Limitations

- **Reference dependency**: Optimized for viral genomes; may need parameter adjustment for other contexts
- **Computational complexity**: Multi-round alignment increases processing time vs. single-pass approaches
- **Memory usage**: Keeps alignment data in memory; very large datasets may require modifications
- **Format assumptions**: Output specifically designed for ViReMa compatibility