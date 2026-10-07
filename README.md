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
- `-lr`: Long read technology (`ont` for Oxford Nanopore, `pb` for PacBio CLR, `hifi` for PacBio HiFi) -- all three are fully implemented and share the same iterative rescue pipeline (see below); they differ only in minimap2 preset.

Seed now acts as a threshold determining whether or not softclips are sent for alignment. Softclip length > Seed -> sent for alignment.

**Defunct flags with `--Aligner minimap2`** (short- or long-read mode): these are standard ViReMa
arguments that only apply to the bowtie/bowtie2/bwa alignment path and have no effect here --
`--X`, `--ThreePad`, `--FivePad`, `--ErrorDensity`, `--MaxIters`, `--Internal_Pad`, `-Windows`,
`-Fasta`, `--Pad`, `--Host_Index`, `--Host_Seed`. `--Host_Index` in particular is a dead end rather
than a no-op: it still builds a host bowtie index, but `Minimap2_Module.py` has no host-alignment
step at all, so no read is ever aligned against it and the Host_* output files stay empty.
`--Aligner_Directory` does work here (it redirects the `minimap2` call the same way it redirects
bowtie/bwa), and `--N` has no effect on the minimap2 alignment itself but still affects substitution
classification during compilation.

## Logic

### 1. Initial Alignment
- **Tool**: Minimap2 with technology-specific parameters:
  - **Short reads (default)**: `-ax sr -k 20 -A 1 -B 2`
  - **Oxford Nanopore**: fully hand-tuned flags (not the bare `-ax map-ont` preset) -- `-k 15 -w 5 -A 1 -B 2 -O 2,32 -E 1,0 -z 200 -g 2000 -Y`, chosen to make long, cheap gaps (deletions) preferable to short expensive ones.
  - **PacBio CLR**: `-ax map-pb`
  - **PacBio HiFi**: `-ax map-hifi`
- **Purpose**: Broad alignment to identify primary mappings and softclipped regions
- **Output**: Primary and supplemental alignments with softclipped portions

Before softclip extraction, any large embedded `I` (insertion) CIGAR operation -- minimap2's usual way of representing a chimeric/foreign-content read as one alignment rather than splitting it -- is rewritten into an explicit trailing softclip (`rewrite_embedded_insertions_as_softclips`), run on every round's output for all three long-read techs. This is what lets the softclip-rescue pipeline below discover chimeric junctions that minimap2 never actually soft-clipped in the first place (cross-locus insertions, cross-segment fusions, tandem duplications, inversions all commonly present this way).

Only the alignment with the highest alignment score has its softclips sent for additional rounds of alignment. Iterative re-alignment now runs identically for all three long-read technologies (`ITERATIVE_RESCUE_TECHS = ('ont','pb','hifi')`) -- earlier, PacBio CLR/HiFi used a single non-iterated round and reused the whole-read preset for softclip rescue, which failed outright on short extracted fragments; both are now unified with ONT's iterative, short-fragment-tuned rescue. Max rounds is 6 (successful rescues converge by round 3-4 in practice).

### 2. Softclip Extraction and Re-alignment
- **Extraction**: Identifies softclipped sequences ≥ threshold length (Seed parameter) from primary alignments
- **Naming Convention**: 
  - `softclip_0`: Softclips occurring **before** the main alignment
  - `softclip_1`: Softclips occurring **after** the main alignment
  - `softclip_1_ins` (etc.): the isolated-insertion half of a split candidate pair -- see below
- **Re-alignment**: Uses technology-specific parameters:
  - **Short reads**: More sensitive parameters (`-k 10 -w 5 -m 10`) to map softclipped regions
  - **ONT / PacBio CLR**: same short-fragment preset as short reads (`-k 10`)
  - **PacBio HiFi**: same short-fragment preset with a larger `-k 15` (lower error rate supports a more specific minimizer)
- **Purpose**: Detect split alignments indicating potential recombination events

When a softclip came from the embedded-insertion rewrite above, it is split into **two independent realignment candidates** instead of one glued fragment: the insertion content alone (`..._softclip_1_ins`) and everything after it alone (`..._softclip_1`). Gluing them together is only unambiguous when the insertion's true origin and the following flank sit in increasing reference order; when the donor sits on the other side, the glued fragment contains a real backward jump internally, and minimap2 resolves it as one confident, silently-wrong placement (anchored on whichever piece is longer) rather than leaving a softclip the rescue mechanism would notice. Realigning the two pieces independently removes that failure mode regardless of which way the jump points. A matching exemption in the redundant-segment dedup step (below) keeps the two halves from being collapsed back into one, since they're expected to overlap heavily in reference space by design.

Sometimes, mapping softclips will have further softclips that exceed threshold length. These are sent for another round of mapping where the read name has an additional softclip_0 or softclip_1 appended to the read name in the intermediate file (`<Output_SAM>_temp`). If these softclips map, they are stitched back to the primary alignment in the correct order. Softclips that remain unmapped but are internal (not at the ends of the final, stitched together mapping) are turned into I events in the CIGAR string. Segments that map to the same locus as an already-kept segment (e.g. a nested softclip rediscovering what its own parent's supplementary already found) are deduplicated, keeping whichever has more aligned query bases. See SOFTCLIP_MERGING_LOGIC.md for more details.

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

As with the `TEMPREADS`/`TEMPSAM1` files of the bowtie/bwa workflow, these are written to `--Output_Dir`
(the current directory if unset) and removed once alignment has finished. With `--Debug`, the
intermediate SAM and `multiRound` are kept.

- **`<Output_SAM>_temp`**: Complete alignment results (initial + merged softclip alignments)
- **`multiRound`**: Filtered reads with qualifying softclipped regions 
- **`TEMP_SAM`**: Secondary alignment results for extracted softclips
- **`TEMP_READS.txt`** (`TEMP_READS_R<n>.txt` in later rounds): FASTA format extracted softclipped sequences

## System Requirements

- **Python**: 3.x with standard libraries (`subprocess`, `re`, `argparse`, `collections`)
- **Minimap2**: Must be available in system PATH

## Recent Changes

Validated via a 120+ run synthetic benchmark suite (deletions, tandem duplications, FHV
cross-segment insertions, same-chromosome cross-locus insertions, and an edge-of-reference
proximity sweep, each across ONT/PacBio CLR/PacBio HiFi) -- see that suite's own README for
methodology and numbers.

- **Insertion-in-isolation rescue** (`Minimap2_Module.py`): see "Softclip Extraction and
  Re-alignment" above. Fixes confidently-wrong placements for cross-locus insertions whose donor
  sits downstream of the acceptor.
- **PacBio CLR/HiFi now share ONT's full iterative rescue pipeline** instead of a single
  non-iterated round with the whole-read preset reused for short-fragment rescue (which failed
  outright on short fragments).
- **Recombination classification fix** (`Compiler_Module.py`, long-read mode only): a Donor/
  Acceptor jump separated by a small residual (below `Seed`) is now classified as `Recombination`/
  `MicroDeletion`/`MicroInsertion` instead of being discarded as `UnknownInsertion`, which
  previously threw away the real event and reported only the residue's own length. The residual is
  also separately recorded as its own `MicroInsertion`. Includes a guard against a related
  pre-existing field-misalignment issue where Donor/Acceptor can pick up raw sequence data instead
  of a reference name for reads with 3+ segments.
- **BED output fix** (`Compiler_Module.py`): BED files were silently empty whenever
  `-ReadNamesEntry` was set (i.e. effectively always, for any run that also wants per-read
  attribution) -- the `-ReadNamesEntry` code path in `WriteFinalDict` was missing the
  `WritetoBEDFile` call the other path already had.
- **Dead `.mmi` index build removed** from `ViReMa.py` -- it was never read anywhere in
  `Minimap2_Module.py`, and reusing it would have silently overridden each round's own `-k`/`-w`
  flags (minimap2 bakes those into the index at build time).
- **`Read_Events.tsv` now written in short-read mode too** (`Compiler_Module.py`): the per-read
  table was gated on `cfg.LongReadTech`, so it was never produced for `bowtie`/`bwa`/`bowtie2`/
  short-read-`minimap2` runs at all, even with `-ReadNamesEntry` set -- discovered when building a
  short-read benchmark that needed it (nothing could score a short-read run against it before this
  fix). The classification logic that populates it (`WriteReadEvent`, called unconditionally
  whenever a read gets a real `EventType`) is generic over every aligner, so there was no reason
  for the gate; changed to `cfg.ReadNamesEntry`, matching the flag's own stated purpose.

### Fixed: strand information was discarded before classification

Traced while investigating inversion detection. `Minimap2_Module.py` correctly finds reverse-
strand alignments for inverted content once it's isolated by the softclip-rescue mechanism
(confirmed directly: an isolated inversion candidate aligned at `flag=16, MAPQ=60` at the true
boundary) -- and `Compiler_Module.py`'s classifier already has working `_RevStrand`-aware logic
for this, used correctly by the short-read pipeline. But three record-construction sites in
`Minimap2_Module.py` unconditionally hard-coded the SAM strand flag to forward before the record
ever reached classification, discarding that correctly-detected information. This is now fixed --
`convert_single_read`, `create_single_merged_record_new`, and `create_paired_records_new` all
preserve strand.

The fix required one subtlety beyond "just stop hard-coding the flag": a segment's own raw SAM
strand bit is relative to the *original as-sequenced read*, not to the shared query frame every
merged/paired record actually uses (`primary_read`'s own SEQ, via `.copy()`). Naively preserving
every segment's own absolute flag caused a real regression on ordinary (non-inversion) DUP/CISINS
events -- roughly half of all reads, since long-read sequencing is strand-agnostic, happen to have
their whole-read primary alignment reported reverse for no event-related reason, which (before the
correction) spuriously flipped the Donor/Acceptor junction side used for position calculation.
`_segment_strand` (`Minimap2_Module.py`) now treats the primary/'main' segment as always forward
in its own frame (correct by construction) and uses independently re-aligned rescue segments'
own reported flag as-is (already relative to that same frame, since their query is a literal
substring of the primary's SEQ). Verified against the DUP-100bp, CISINS-1000bp, and DEL-1000bp
benchmark cells: TP/FP now match the pre-strand-fix baseline exactly (DUP 19/64, CISINS 108/149)
or are consistent with prior validated performance (DEL 65/66).

**Separately, and not fixed by the above**: testing this against a real synthetic inversion (a
1000bp reverse-complemented block in a ~30kb reference) through the actual production entry point
(not a hand-rolled pipeline replay) shows the ONT alignment preset's very cheap gap costs (`-O
2,32 -E 1,0 -z 200 -g 2000`, tuned to prefer long deletions) make it *cheaper* for minimap2 to
plow straight through the ~1000bp reverse-complemented block as a dense run of small-indel noise
within one alignment record than to soft-clip it out. With no softclip produced, the rescue
mechanism (and therefore this strand fix) never gets a chance to act -- there's nothing to
re-align. So while strand information is no longer *discarded* once found, genuine inversion
detection for an inversion of this size is still not demonstrated to work under the current
preset; that is a separate, alignment-tuning problem.
