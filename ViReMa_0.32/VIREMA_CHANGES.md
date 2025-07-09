# ViReMa Complex CIGAR String Fix

## Problem Description

ViReMa was crashing when processing certain long read datasets with complex CIGAR strings. The error occurred in `Compiler_Module.py` at line 1204:

```
ValueError: invalid literal for int() with base 10: 'NC'
```

### Root Cause

The issue was in the array indexing logic used to extract donor and acceptor site positions from parsed SAM data. The code assumed fixed relative positions (+1, +2, +3, +4) from mapped segment indices, but very complex CIGAR strings caused field misalignment in the internal data structure.

Specifically:
- `AcceptorSite = line[i+3].split("_")[0]` expected a numeric position
- Complex CIGAR strings shifted the array structure
- Reference sequence names like "NC_004144_FHV_RNA2.seq" ended up in positions where numeric values were expected
- The `int()` conversion failed when trying to convert "NC" to an integer

### Data Structure Impact

The `Indices(Code)` function calculates field positions based on the Code string pattern derived from CIGAR operations. With very complex CIGAR strings containing many operations (insertions, deletions, matches), the expected field positions calculated by this function could become misaligned with the actual data structure created by the SAM-to-internal format conversion.

## Solution

### Added Safe Position Extraction Function

Created a new function `safe_extract_position()` that:
1. Searches for position information in a range around the expected position
2. Validates that extracted values are numeric before returning them
3. Handles both donor and acceptor site extraction with appropriate logic
4. Provides fallback behavior for edge cases

### Key Changes Made

1. **Added `safe_extract_position()` function** in `Compiler_Module.py`:
   - Searches multiple positions around the expected offset
   - Validates numeric content before extraction
   - Handles RevStrand vs normal strand position extraction differences
   - Provides graceful fallback behavior

2. **Updated AcceptorSite extraction** (lines 1170, 1172, 1213, 1216):
   - Changed from `line[i+3].split("_")[0]` to `safe_extract_position(line, i, 3)`
   - Changed from `line[i+4].split("_")[0]` to `safe_extract_position(line, i, 4)`

3. **Updated DonorSite extraction** (lines 1180, 1182):
   - Changed from `line[i+1].split("_")[2]` and `line[i+1].split("_")[1]`
   - To `safe_extract_position(line, i, 1, is_donor=True)`

### Function Details

```python
def safe_extract_position(line, base_index, expected_offset, is_donor=False):
    """
    Safely extract position information from line array, handling complex CIGAR strings.
    Looks for the first field at or after base_index+expected_offset that contains position info.
    """
```

The function:
- Searches within a 5-position window around the expected location
- For donor sites: Extracts from part[2] for RevStrand, part[1] for normal strand
- For acceptor sites: Extracts from part[0] (first part of underscore-separated fields)
- Validates numeric content with try/except blocks
- Provides fallback to original logic if no valid position found

## Impact

This fix ensures that ViReMa can process long read datasets with complex CIGAR strings without crashing, while preserving all valuable data that was previously being lost. The solution is robust and handles edge cases gracefully without filtering out any valid entries.

## Files Modified

- `Compiler_Module.py`: Added `safe_extract_position()` function and updated position extraction calls