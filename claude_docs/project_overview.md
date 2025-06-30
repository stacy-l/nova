# Nova Project Overview

## Executive Summary

`nova` is a de novo variant insertion simulator designed to evaluate structural variant detection tools' ability to detect ultra-low frequency insertions in long-read sequencing data. The simulator generates random, simple repeat, and predefined insertion sequences and inserts them into BAM reads at ultra-low frequencies, creating ground truth datasets for benchmarking variant callers.

## Scientific Background

### The Challenge
Detecting ultra-low frequency structural variants (< 1% allele frequency) in long-read sequencing data presents unique challenges:
- Traditional variant callers assume moderate to high allele frequencies
- Low-frequency variants may appear in only 1-2 reads across an entire dataset
- Distinguishing true variants from sequencing artifacts becomes critical
- Clustering algorithms may incorrectly merge independent events

### Nova's Approach
Nova addresses these challenges by:
1. Inserting known sequences at controlled ultra-low frequencies
2. Tracking exact insertion locations and sequences for ground truth
3. Distributing insertions genome-wide to avoid clustering artifacts
4. Supporting multiple insertion types to test different detection scenarios

## Key Features

### Insertion Types
1. **Random Sequences**: Configurable length and GC content
2. **Simple Tandem Repeats**: CAG, CTG, and other repeat expansions
3. **Predefined Sequences**: Alu elements, LINE-1 segments, custom FASTA

### Simulation Control
- **Read Selection**: Filters by length, quality, and soft-clip ratio
- **Window Limiting**: Prevents genomic clustering with `--max-reads-per-window`
- **Proportional Distribution**: Allocates insertions across chromosomes by size
- **Reproducibility**: Seed-based random number generation

### Output Formats
- Modified BAM file with inserted sequences
- FASTA file of modified reads
- JSON registry of all insertions
- Detailed insertion records with genomic coordinates

## Use Cases

### Primary Applications
1. **Benchmarking SV Callers**: Test sensitivity at ultra-low frequencies
2. **Method Development**: Create ground truth for new detection algorithms
3. **Parameter Optimization**: Tune variant caller settings for rare variants
4. **Quality Assessment**: Evaluate false positive rates and clustering behavior

### Research Applications
- Somatic mutation detection in cancer
- Mosaic variant identification
- Mobile element insertion profiling
- Repeat expansion detection

## Project Status

### Current Capabilities
- ✅ Three insertion type categories
- ✅ Window-based clustering prevention
- ✅ Comprehensive output tracking
- ✅ Snakemake pipeline integration
- ✅ Analysis scripts for validation

### Recent Improvements
- Window limit parameter eliminates false positives
- Memory optimization for large-scale simulations
- Enhanced chromosome-proportional distribution
- Detailed false positive categorization

### Future Directions
- Deletion and duplication support
- Complex structural variant patterns
- Multi-sample simulation
- Integration with additional variant callers

## Related Documentation
- [Architecture Details](architecture.md)
- [Window Limit Experiments](window_limit_experiments.md)
- [Development Setup](development_setup.md)
- [Testing Philosophy](testing_philosophy.md)