# Region-Aware Sampling Implementation

## Overview

This document describes the complete rewrite of Nova's region filtering system to address performance issues with the previous per-read checking approach. The new implementation pre-computes valid sampling windows based on BED file regions, eliminating the need to check every read against region filters.

## Problem Statement

The original implementation checked every selected read against inclusion/exclusion regions:
- **Performance**: O(n) region checks where n = number of reads
- **Inefficiency**: Millions of redundant overlap calculations
- **Unusable**: System becomes impractical for large-scale simulations

## Solution Architecture

### Core Concept
Instead of checking regions per-read, we:
1. Pre-compute valid genomic regions based on inclusion/exclusion BED files
2. Generate sampling windows only within valid regions
3. Select reads from these pre-validated windows

### Key Components

#### 1. RegionAwareWindowGenerator
A new class that handles all region-based window generation:

```python
class RegionAwareWindowGenerator:
    def __init__(self, bam_path: str):
        self.bam_path = bam_path
        self._valid_regions_cache = {}  # Cache PyRanges, not windows
        self._genome_regions = None     # Lazy-loaded
```

**Responsibilities:**
- Compute valid regions from inclusion/exclusion filters
- Cache expensive region operations (intersections/subtractions)
- Generate random windows within valid regions
- Ensure sufficient windows for requested read count

#### 2. Valid Region Computation

**Scenarios:**
1. **No filters**: Use entire genome
2. **Exclusion only**: Genome minus exclusion regions
3. **Inclusion only**: Only included regions
4. **Both**: Inclusion regions minus exclusions

**Caching Strategy:**
- Cache the resulting PyRanges objects (valid regions)
- Don't cache individual windows (respects random seed)
- Cache key: hash of BED file paths

#### 3. Window Generation Algorithm

```python
def generate_windows_for_regions(valid_regions, n_windows, window_size, seed):
    # Step 1: Calculate proportional allocation
    total_length = sum(region lengths)
    for region in valid_regions:
        proportion = region_length / total_length
        windows_for_region = int(proportion * n_windows)
    
    # Step 2: Generate windows with random positions
    random.seed(seed)  # Respect the seed
    for region, num_windows:
        for _ in range(num_windows):
            start = random.randint(region.start, region.end - window_size)
            windows.append((chr, start, start + window_size))
```

### Parameter Changes

**Rename**: `--max-reads-per-window` → `--reads-per-window`
- Clarifies this is the exact number, not a maximum
- Better reflects the 1:1 relationship with windows in current usage

### Edge Case Handling

1. **Insufficient Valid Regions**
   - Log warning with statistics
   - Generate as many windows as possible
   - Let downstream code handle shortfall

2. **Small Regions**
   - If region < window_size, use entire region as window
   - Adjust reads_per_window expectation accordingly

3. **Fragmented Regions**
   - Many small disjoint regions after exclusions
   - May need to adjust window_size dynamically

## Implementation Plan

### Phase 1: Core Infrastructure
1. Create `RegionAwareWindowGenerator` class
2. Implement valid region computation with caching
3. Add window generation algorithm

### Phase 2: Integration
1. Rename parameter throughout codebase
2. Rewrite `ReadSelector._lazy_select()`
3. Remove obsolete per-read checking code

### Phase 3: Polish
1. Add comprehensive logging
2. Performance metrics
3. Edge case handling

## Performance Implications

### Before
- Check every read: O(n × m) where n = reads, m = regions
- No caching of region operations
- Repeated PyRanges operations

### After
- Pre-compute once: O(m) for region operations
- O(w) to generate w windows
- No per-read checks needed

### Expected Improvements
- 100-1000x faster for large read counts
- Reduced memory usage (no repeated PyRanges objects)
- More predictable performance

## Testing Strategy

### Unit Tests
- Window generation with various region configurations
- Cache behavior with different seeds
- Edge cases (small regions, no valid regions)

### Integration Tests
- End-to-end with real BED files
- Verify reads come from correct regions
- Performance benchmarks

### Validation
- Compare output distribution with old implementation
- Ensure proper randomization
- Verify inclusion/exclusion logic

## Migration Notes

Since backwards compatibility is not a concern:
- Complete removal of old region checking code
- Clean API without legacy parameters
- Simplified ReadSelector logic

## Future Considerations

1. **Dynamic Window Sizing**
   - Adjust window size based on region characteristics
   - Optimize for fragmented regions

2. **Parallel Window Generation**
   - Generate windows for different chromosomes in parallel
   - Further performance improvements

3. **Advanced Caching**
   - Persist cache between runs
   - Share cache across processes

## Status Updates

- **2025-01-02**: ✅ **Implementation Complete**
  - `RegionAwareWindowGenerator` class implemented in `region_utils.py`
  - Parameter renamed from `--max-reads-per-window` to `--reads-per-window`
  - `ReadSelector._lazy_select()` completely rewritten to use region-aware windows
  - Removed obsolete per-read region checking code (`_passes_region_filters()`)
  - Added comprehensive test suite using real BAM and BED files
  - All 38 region-related tests passing
  - Backward compatibility maintained: no regions = genome-wide sampling

## Implementation Results

### Performance Improvements
- **Eliminated O(n) per-read region checks** where n = number of reads
- **Pre-computed valid regions** cached for reuse across variant groups
- **Window generation scales with region count**, not read count

### Key Features Delivered
1. **Region-aware window generation** that respects inclusion/exclusion regions
2. **Proportional sampling** within valid regions
3. **Caching system** for expensive region operations
4. **Comprehensive logging** for window generation statistics
5. **Robust edge case handling** (small regions, no valid regions, etc.)

### Test Coverage
- **14 unit tests** for `RegionAwareWindowGenerator`
- **9 integration tests** for `ReadSelector` region functionality  
- **15 existing tests** for `RegionFilter` (unchanged)
- **Real data testing** using test BAM file and subset BED files

### Files Modified
- `src/nova/region_utils.py`: Added `RegionAwareWindowGenerator` class
- `src/nova/read_selector.py`: Rewrote `_lazy_select()`, removed obsolete methods
- `src/nova/cli.py`: Renamed parameter `--max-reads-per-window` → `--reads-per-window`
- `tests/test_data/`: Added test BED files (`test_centromeres.bed`, `test_genomic_superdups.bed`)
- `tests/test_region_aware_window_generator.py`: New comprehensive test suite
- `tests/test_read_selector_regions.py`: Updated for new architecture

The region-aware sampling system is now ready for production use.