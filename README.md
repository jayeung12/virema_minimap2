# ViReMa with Minimap2

ViReMa with minimap2 as an alignment option, for both short reads and long reads (ONT, PacBio CLR, PacBio HiFi). Everything downstream of SAM compilation works as in standard ViReMa.

## Requirements

- Python 3
- minimap2 on `PATH`, or in the directory given by `--Aligner_Directory`
- bowtie (`bowtie-build` and `bowtie-inspect`). The compile step reads reference data from a bowtie index even when the aligner is minimap2.
- samtools, only for `-BAM`, `-Only_BAM` and `-BED12`

## Usage

```bash
python ViReMa_0.32/ViReMa.py [options] Virus_Index Input_Data Output_SAM
```

Positional arguments:

| Argument | Description |
|---|---|
| `Virus_Index` | Virus reference genome in FASTA format. The index is built automatically. |
| `Input_Data` | Reads in FASTQ format. minimap2 also accepts FASTA directly; the `-Fasta` flag only applies to the other aligners. |
| `Output_SAM` | Name of the output SAM file. |

### Examples

Run these from `ViReMa_0.32/`. They use the files in `Test_Data/`.

Short reads (FHV, Illumina):

```bash
python ViReMa.py --Aligner minimap2 --Seed 25 --MicroInDel_Length 20 \
    Test_Data/FHV_Genome.txt Test_Data/FHV_small.txt output_mm2.sam

python ViReMa.py --Aligner minimap2 --Seed 25 --MicroInDel_Length 20 \
    Test_Data/FHV_Genome.txt Test_Data/FHV_10k.txt output_mm2_10k.sam
```

Oxford Nanopore reads (FHV):

```bash
python ViReMa.py --Aligner minimap2 -lr ont --Seed 25 --MicroInDel_Length 20 \
    Test_Data/FHV_Genome.txt Test_Data/FHV_nanopore.fastq output_mm2_ont_fhv.sam
```

Oxford Nanopore reads (SARS-CoV-2):

```bash
python ViReMa.py --Aligner minimap2 -lr ont --Seed 25 --MicroInDel_Length 20 \
    Test_Data/SARS2_Genome.fasta Test_Data/ont_reads_0001.fq output_mm2_ont.sam
```

PacBio CLR reads (SARS-CoV-2):

```bash
python ViReMa.py --Aligner minimap2 -lr pb --Seed 25 --MicroInDel_Length 20 \
    Test_Data/SARS2_Genome.fasta Test_Data/pacbio_reads_0001.fq output_mm2_pb.sam
```

## Options

### Alignment

| Flag | Description |
|---|---|
| `--Aligner` | `bowtie` (default), `bowtie2`, `bwa` or `minimap2`. |
| `-lr` | Long-read technology for minimap2: `ont`, `pb` (PacBio CLR) or `hifi` (PacBio HiFi). Omit for short-read mode. |
| `--Aligner_Directory` | Directory containing the aligner executable, if it is not on `PATH`. |
| `--Chunk` | Number of reads to process together. Default 1M. |
| `--Pad` | Number of A's to add to the 3' end of the viral genome. Not used with minimap2. |

### Error handling

| Flag | Description |
|---|---|
| `--Seed` | Default 25. With minimap2, this is the minimum softclip length that is sent for realignment. |
| `--N` | Mismatches tolerated in seeds and mapped segments. Default 1. With minimap2 it does not change the alignment, but it still affects substitution calls during compilation. |
| `--MicroInDel_Length` | Insertions and deletions up to this size are reported as MicroInDels. Default 0. |
| `--BackSplice_limit` | Minimum size of back-splice or duplication that is reported. Default 0. |
| `--Compound_Handling` | Number of nucleotides to map for compound recombination events. Must be less than Seed and greater than MicroInDel_Length. Default off. |
| `--X`, `--ThreePad`, `--FivePad`, `--ErrorDensity`, `--Internal_Pad` | Bowtie/BWA-only search parameters. Not used with minimap2. |

### Fuzz handling
May not work as intended with minimap2.

| Flag | Description |
|---|---|
| `--Defuzz` | Where to report fuzzy junctions: `5` (5' end), `3` (3' end) or `0` (centre). Default off. |
| `--MaxFuzz` | Maximum fuzzy region length. Longer events are not reported. Default Seed. |
| `-FuzzEntry` | Append the fuzz of each recombination event to the results. |

### Output

| Flag | Description |
|---|---|
| `--Output_Dir` | Directory for all output files. |
| `--Output_Tag` | Tag appended to each output file name. |
| `-Overwrite` | Allow overwriting previous output. |
| `-BAM` | Produce a sorted BAM file and BEDGraph with samtools. |

### Compilation

| Flag | Description |
|---|---|
| `-ReadNamesEntry` | Append the names of contributing reads to each result, and write `Read_Events.tsv`. |
| `-DeDup` | Remove potential PCR duplicates. |
| `--UMI` | String or delimiter in read names that marks the UMI barcode. |
| `-BED` | Write recombination data to BED files. |
| `-BED12` | Write BED12 files with coverage and nucleotides. Forces BAM output. |
| `-Stranded` | Compute coverage for + and - strands separately. Requires `-BED12`. |
| `--PseudoRef` | Length of reference segment used to build a pseudo-reference for each event in the BED file. Default off. |
| `--Coverage_Offset` | Distance from the junction at which read coverage is reported in the BED file. Default 0. |
| `-Only_Map` | Align only. Skip compilation. |
| `-Only_BAM` | Build a sorted BAM and coverage from an existing output SAM. Skips alignment and compilation. |
| `-Only_Compile` | Recompile an existing ViReMa output with new parameters. |

### Runtime

| Flag | Description |
|---|---|
| `--p` | Number of processors. Default 1. |
| `--Debug` | File to write debug output to. Also keeps minimap2's intermediate files. |
| `--Debug_Func` | Comma-separated list of functions to debug. Requires `--Debug`. |
| `--MaxIters` | Maximum number of iterations. Not used with minimap2. |
| `-Windows` | Run from a Windows/Cygwin shell. Not used with minimap2. |

### Host

| Flag | Description |
|---|---|
| `--Host_Index`, `--Host_Seed` | Host genome alignment. Not supported with minimap2: the host index is built, but no reads are aligned to it. |

For details of how minimap2 softclips are realigned and merged, see the preprint.
