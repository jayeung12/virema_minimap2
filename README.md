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
  - **Short reads (default)**: `-ax sr -k 20 -A 1 -B 2 -O 2,8 -g 2000 -z 800,400 -n 1 -p 0.3 -N 3 -s 20 --end-bonus 0` (`-k` overridable via the `VIREMA_MM2_SR_K` env var, a leftover short-read tuning knob)
  - **Oxford Nanopore**: fully hand-tuned flags (not the bare `-ax map-ont` preset) -- `-k 15 -w 5 -A 1 -B 2 -O 2,32 -E 1,0 -z 200 -g 2000 -Y`, chosen to make long, cheap gaps (deletions) preferable to short expensive ones.
  - **PacBio CLR**: `-ax map-pb`
  - **PacBio HiFi**: `-ax map-hifi`
- **Purpose**: Broad alignment to identify primary mappings and softclipped regions
- **Output**: Primary and supplemental alignments with softclipped portions

Before softclip extraction, any large embedded `I` (insertion) CIGAR operation -- minimap2's usual way of representing a chimeric/foreign-content read as one alignment rather than splitting it -- is rewritten into an explicit trailing softclip (`rewrite_embedded_insertions_as_softclips`), run on every round's output for all three long-read techs. This is what lets the softclip-rescue pipeline below discover chimeric junctions that minimap2 never actually soft-clipped in the first place (cross-locus insertions, cross-segment fusions, tandem duplications, inversions all commonly present this way).

Only the alignment with the highest alignment score has its softclips sent for additional rounds of alignment (all modes, short reads included). Iterative re-alignment now runs identically for all three long-read technologies (`ITERATIVE_RESCUE_TECHS = ('ont','pb','hifi')`) -- earlier, PacBio CLR/HiFi used a single non-iterated round and reused the whole-read preset for softclip rescue, which failed outright on short extracted fragments; both are now unified with ONT's iterative, short-fragment-tuned rescue. Max rounds is 6 for long reads (`LONG_READ_MAX_ROUNDS`; successful rescues converge by round 3-4 in practice). Short reads use the same iterative loop, capped at round 4 (`SHORT_READ_MAX_ROUNDS`, i.e. up to 3 rescue rounds), with their own rescue preset and without the embedded-insertion rewrite.

### 2. Softclip Extraction and Re-alignment
- **Extraction**: Identifies softclipped sequences ≥ threshold length (Seed parameter) from primary alignments
- **Naming Convention**: 
  - `softclip_0`: Softclips occurring **before** the main alignment
  - `softclip_1`: Softclips occurring **after** the main alignment
  - `softclip_1_ins` (etc.): the isolated-insertion half of a split candidate pair -- see below
- **Re-alignment**: all four modes use an `-ax sr`-based preset for this step (even for long-read
  data -- see rationale below), but **not the same preset as each other**:
  - **Short reads**: `-k 10 -w 5 -m 10 -n 2 -A 1 -B 2 -O 12,32 -E 2,1 --end-bonus 5 -s 20 -z 200 -r 50`
  - **ONT / PacBio CLR**: `-k 10 -w 5 -m 10 -n 2 -A 2 -B 2 -O 2,8 -E 2,1 --end-bonus 5 -s 20 -z 200 -r 50`
  - **PacBio HiFi**: same as ONT/CLR with a larger `-k 15` (lower error rate supports a more specific minimizer)
  - Short reads' own gap cost (`-O 12,32`) is considerably stricter than the long-read modes'
    (`-O 2,8`) -- they share the general shape (`-k 10`, tight `-w`/`-z`/`-r`) but are not
    interchangeable; the long-read value was raised from an original `2,4` this session (see
    Recent Changes) and was never shared with the short-read preset, which was independently
    already strict.
- **Purpose**: Detect split alignments indicating potential recombination events

When a softclip came from the embedded-insertion rewrite above, it is split into **two independent realignment candidates** instead of one glued fragment: the insertion content alone (`..._softclip_1_ins`) and everything after it alone (`..._softclip_1`). Gluing them together is only unambiguous when the insertion's true origin and the following flank sit in increasing reference order; when the donor sits on the other side, the glued fragment contains a real backward jump internally, and minimap2 resolves it as one confident, silently-wrong placement (anchored on whichever piece is longer) rather than leaving a softclip the rescue mechanism would notice. Realigning the two pieces independently removes that failure mode regardless of which way the jump points. A matching exemption in the redundant-segment dedup step (below) keeps the two halves from being collapsed back into one, since they're expected to overlap heavily in reference space by design.

Sometimes, mapping softclips will have further softclips that exceed threshold length. These are sent for another round of mapping where the read name has an additional softclip_0 or softclip_1 appended to the read name in the intermediate file (`<Output_SAM>_temp`). If these softclips map, they are stitched back to the primary alignment in the correct order. Softclips that remain unmapped but are internal (not at the ends of the final, stitched together mapping) are turned into I events in the CIGAR string. Segments that map to the same locus as an already-kept segment (e.g. a nested softclip rediscovering what its own parent's supplementary already found) are deduplicated, keeping whichever has more aligned query bases. See SOFTCLIP_MERGING_LOGIC.md for more details.

### 3. Result Merging and Classification

#### Selecting Among Several Candidate Alignments
Before a read reaches the single/multiple classification below, there are three separate points where the pipeline has to pick among competing candidate alignments for the same read or softclip fragment, each with different rules:

1. **Round-1 multi-primary (`parse_sam_for_softclipped`)**: covered above (line 49) -- in every mode (short reads included), if round 1 produced more than one primary-like record for a read, only the highest-`AS:i:`-scoring one has its softclips queued for rescue.
2. **Merging a rescue round back in (`select_alignments_for_identifier`)**: when a softclip identifier comes back from a rescue round with several competing records, the best-scoring PRIMARY record is kept **together with every SUPPLEMENTARY record** (`FLAG & 2048`) -- supplementaries are never dropped in favor of a single best-AS winner, since a genuine chimeric/split alignment needs all of its pieces. Only genuinely-competing SECONDARY records (`FLAG & 256`) are discarded. If no primary exists at all among the candidates (rare), this falls back to plain single-best-AS across everything.
3. **Classifying records for the final conversion (`merge_split_alignments`)**: records are bucketed into "softclip segment" (name contains `_softclip_`) or "primary candidate" (no `_softclip_` in the name, not secondary/supplementary); the primary candidate with the highest `AS:i:` becomes the main alignment. **Known gap**: a record that IS supplementary but does NOT have `_softclip_` in its name -- i.e. a genuine round-1 primary+supplementary pair minimap2 already got right on its own -- matches neither bucket and is silently dropped. This only surfaces when the same read also needs rescue somewhere else (so round 1's output for that read is a mix of clean primary+supplementary plus still-unresolved softclips); see "Mechanism B" under Recent Changes below. Separately, if one softclip position ever collects 3+ candidate records (only possible via the insertion/rest 2-way split producing an unexpected extra fragment), the code falls back to plain single-best-AS across all of them rather than guessing at an N-way merge order.

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

- **Short-read mode now selects by best AS and iterates up to 3 rescue rounds**
  (`Minimap2_Module.py`): round-1 softclip selection used the first non-supplementary record for
  short reads and best-`AS:i:` only for long reads; it is now best-`AS:i:` everywhere (ties still
  go to minimap2's primary, which it writes first). Short reads also realigned softclips exactly
  once, so a rescued clip's own clip >= Seed was never placed; they now run through the same
  iterative loop as long reads, capped at round 4 (`SHORT_READ_MAX_ROUNDS`), without the
  long-read-only embedded-insertion rewrite. On `Test_Data/FHV_10k.txt` (`--Seed 25
  --MicroInDel_Length 20`): the selected record was unchanged for all 7,305 rescued reads, 5
  nested clips reached round 3 and the loop converged there, 4 of them were deduplicated against
  their parent's own supplementary, and 1 read gained a placed 25-nt tail
  (`40M960N25M25S` -> `40M960N25M12N25M`, one new micro-deletion). Long-read output and the
  250-read `FHV_small.txt` output are byte-identical to before.

- **Insertions sitting at a recombination junction are now reported** (`Minimap2_Module.py`,
  `Compiler_Module.py`): an insert only ever reached `Virus_Insertions.txt` and
  `Read_Events.tsv` when its flanking coordinates happened to be adjacent
  (`acceptor == donor + 1`). An insert at a junction was routed to the `UnknownRecombination`
  counter, whose output file is commented out, so the junction *and* the inserted sequence were
  both discarded. Found on a NanoLuc reporter influenza library (SRR34426338), where a 577nt
  NanoLuc+2A cassette between NS 657 and NS 1 appeared in no output file at all while an 88nt
  insert between NS 724 and 725 was reported normally. Three separate defects:
  - **The inserted bases were not in the SAM's coordinate system.** When a fragment's own
    realignment failed, `merge_split_alignments()` dropped it (`flag & 4`), and the segment that
    would have carried its bases as a soft clip could then be removed by
    `_dedupe_redundant_segments()` in favour of an equivalent record that HARD-clips them. The
    merged records described a shorter read than was sequenced, and `SEQ` stopped matching
    `CIGAR` (105 of 201 records on a 120-read fixture; samtools refused the file outright).
    Unplaced fragments are now kept and reconciled back into the merged CIGAR as an explicit soft
    clip at the position their softclip path puts them in.
  - **A trailing soft clip before a terminal hard clip was rewritten to `I`**, because the
    "internal softclip" test counted the terminal clip as a bounding operation. That produces an
    insertion hanging off the end of an alignment (invalid SAM) and splits one insert across two
    records.
  - **Two adjacent unmapped fields shifted every later field.** `ResultsSort()` reads an insert
    as `line[i+2]` and its acceptor as `line[i+3]`, i.e. it assumes at most one unmapped field
    between two mapped segments; a split insert gave it two, so raw sequence landed in the
    acceptor-reference column (607 rows on SRR34426338). Runs are now merged by
    `CoalesceUnmappedFields()`.

  A junction carrying an insert now emits **two rows** in `Read_Events.tsv` at the same
  `segment_index` -- the junction, classified by its own jump size exactly as a clean junction
  is, and the insert -- and the sequence goes to `Virus_Insertions.txt` keyed on the same
  donor/acceptor pair. No new event categories. Measured on the full SRR34426338 library
  (145,547 reads): 4,722 insert rows at the cassette junction with a **median size of 577**, the
  cassette's exact length, and every sampled recorded sequence aligns to it; the 3,312 reads
  carrying the junction still report it; rows with raw sequence in the acceptor-reference column
  went 607 -> 0; records whose `SEQ` disagreed with their `CIGAR` went 16,531 -> 92 of 162,712;
  "Unknown or Ambiguous Recombination Events" went 10,240 -> 0. On the bundled FHV short-read
  test data a 979nt deletion carrying a 13nt insert was previously discarded entirely and is now
  reported -- the only change to that dataset's results, and its SAM stays byte-identical. On
  SRR37793484 (HiFi) the 31 "Unknown or Ambiguous" events became 31 recombinations plus 32
  insertions, and invalid records went 4 -> 1.

  **Not fixed by this change**: the residual 92 records (and that 1 HiFi read) are a separate
  pre-existing defect in which two groups claim the same query bases, so the CIGAR over-claims
  the read and there is no soft clip left to trim against. samtools therefore still refuses the
  full SRR34426338 output, though at 0.06% of records rather than 10.2%.
  Also guarded: the single-interval BED writes name ONE reference but were handed the acceptor
  coordinate of cross-segment junctions, which indexes past the end of the coverage array in
  `BEDGraph_Plot()`. Previously unreachable because those events were being discarded; the guard
  the original author sketched in a comment on that line is now active, comparing bare reference
  names so same-segment strand flips keep their existing behaviour.
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
- **Round-2+ long-read rescue gap cost raised (`-O 2,4` -> `2,8`)** (`Minimap2_Module.py`, ONT/PB/
  HiFi only): a cheap long-gap-open cost let this preset "limp through" a real second-locus
  fragment via many small indels instead of soft-clipping it for proper rescue, silently absorbing
  true junctions into ordinary-looking noise. Verified directly against raw sequence on two
  independent misses (a ~194bp insertion-donor fragment, and a full ~1000bp inversion, both HiFi
  combo reads) -- raising the cost to 8 recovered a clean second-locus split in both. Re-run
  against the full benchmark suite: multi-event sensitivity improved substantially for inversions
  and insertions (e.g. HiFi `multi_inv` 0.231->0.585, ONT `multi_ins` 0.357->0.875), specificity
  improved across every multi-event category, with small mixed deltas elsewhere and zero effect on
  short-read results (that mode's own round-2+ preset is a different, already-strict value,
  untouched by this change -- see "Softclip Extraction and Re-alignment" above). Not a complete
  fix -- see "Separately, and not fixed by the above" below.
- **`--Aligner_Directory` now works with `--Aligner minimap2`**: every `build_minimap2_command()`
  branch hardcoded the bare string `'minimap2'`, so this flag (which already worked for
  bowtie/bowtie2/bwa) silently had no effect for minimap2. Added the same directory-prefix
  convention ViReMa.py itself uses. Verified end-to-end with `PATH` stripped of `minimap2` and the
  whole bowtie family, redirected entirely through `--Aligner_Directory`, producing identical
  output to a normal `PATH`-based run.
- **minimap2 now builds the bowtie index it still needs** (`ViReMa.py`): the virus-index-building
  step had no `minimap2` case, so it silently skipped building or even checking for an index.
  minimap2 itself doesn't need one, but `Compiler_Module.ExtractRefData()` calls `bowtie-inspect`
  on `Lib1` for every aligner, so a bowtie-format index was still required -- a fresh reference
  would align fine under `--Aligner minimap2` and then fail at the compile step. Added `minimap2`
  to the same build condition `bowtie` uses. Verified against a fresh reference with no
  pre-existing index.

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

**Separately, and not fixed by the above**: full-genome benchmark testing (the combo reference --
DEL+DUP+in-place INV+cross-locus INS all on one molecule, so a read crossing the inversion
typically crosses other real junctions too) found inversion detection incomplete, traced to two
distinct, confirmed mechanisms -- not one.

*Mechanism A -- the initial alignment absorbs the inversion as noise.* Confirmed on a synthetic
1000bp reverse-complemented block: a cheap-enough gap cost lets minimap2 plow straight through it
as a dense run of small-indel noise within one alignment record rather than soft-clipping it out.
With no softclip produced, the rescue mechanism never gets a chance to act -- there's nothing to
re-align, so this strand fix has nothing to preserve.

*Mechanism B -- the initial alignment finds it correctly, and later steps still discard it.*
Confirmed directly on a real simulated HiFi read (combo reference): minimap2's raw round-1 output
was three clean records -- a primary covering the DEL+DUP region, and two separate,
correctly-placed supplementary alignments covering the inversion and what follows it (right
strand, right position, low edit distance on each). Two things then threw this away: (i)
`rewrite_embedded_insertions_as_softclips()` found an unrelated embedded `I` earlier in the same
primary record (the ordinary way a tandem duplication's second copy gets encoded, nothing to do
with the inversion) and truncated the primary there, discarding its own already-correct context;
(ii) `merge_split_alignments()`'s segment-classification loop only recognizes a plain primary or a
round-2-rescue fragment (`_softclip_`-named) -- a genuine round-1 supplementary with neither trait
matches no branch and is silently dropped, which happens whenever the same read also needs any
rescue elsewhere (near-guaranteed for a read crossing multiple real junctions). The two correctly
found supplementaries vanish before compilation ever sees them, ending up absorbed into the same
kind of noise as Mechanism A even though the aligner never got it wrong.

The round-2+ gap-cost fix above (`-O 2,4` -> `2,8`) measurably helps fragments that *do* reach
round-2 rescue -- a real, separate win -- but doesn't touch either mechanism here, since both
discard the result before rescue (or round 1's own supplementaries) are ever consulted. Benchmark
evidence after that fix: combo `multi_inv` sensitivity is 0.585 (HiFi) and 0.729 (PacBio CLR) --
meaningfully above zero, so detection isn't uniformly broken -- but ONT's is essentially unchanged
at 0.303, and none of the three techs are close to DEL's ~0.9+ on the same reference. So while
strand information is no longer *discarded* once found (that fix remains valid and verified),
genuine inversion detection is still incomplete, for reasons beyond alignment-parameter tuning
alone. Neither mechanism above is fixed in code.
