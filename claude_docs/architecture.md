# Nova Architecture Guide

## Overview

Nova follows a modular architecture designed for extensibility and clear separation of concerns. The system processes BAM files, generates variant sequences, inserts them into reads, and produces comprehensive output tracking.

## Core Components

### 1. Read Selection Module (`read_selector.py`)

**Purpose**: Intelligently selects reads from BAM files for variant insertion.

**Key Classes**:
- `ReadSelector`: Main class handling read selection logic
- `ReadMetadata`: Data class storing read properties

**Key Methods**:
- `_select_reads_unified()`: Window-based selection preventing clustering
- `_get_chromosome_proportional_targets()`: Fair distribution across chromosomes
- `_passes_filters()`: Quality, length, and soft-clip filtering
- `_calculate_gc_content()`: GC content analysis for reads

**Critical Features**:
- **Window-based limiting**: Prevents genomic clustering (max_reads_per_window)
- **Proportional allocation**: Distributes reads by chromosome size
- **Quality filtering**: MAPQ, read length, soft-clip ratio thresholds
- **Dual strategies**: Window-limited (<500 reads) vs proportional (≥500 reads)

### 2. Variant Generation Module (`variant_generator.py`)

**Purpose**: Creates insertion sequences of various types.

**Key Classes**:
- `VariantGenerator`: Factory for creating variants
- `Variant`: Data class for variant properties

**Variant Types**:
1. **Random**: 
   - Configurable length and GC content
   - Uses weighted nucleotide selection
   
2. **Simple Repeats**:
   - Tandem repeats (CAG, CTG, etc.)
   - Configurable unit and count
   
3. **Predefined**:
   - Loads from FASTA files
   - Supports Alu, LINE-1, custom sequences

### 3. Variant Registry (`variant_registry.py`)

**Purpose**: Tracks all generated variants with unique identifiers.

**Key Features**:
- UUID-based variant identification
- Metadata storage (type, source, properties)
- JSON serialization for persistence
- Deduplication support

**Data Structure**:
```python
{
    "variant_id": "uuid",
    "sequence": "ATCG...",
    "type": "random|simple|predefined",
    "metadata": {
        "length": 300,
        "gc_content": 0.5,
        "source": "generated"
    }
}
```

### 4. Read Insertion Module (`read_inserter.py`)

**Purpose**: Inserts variants into selected reads at random positions.

**Key Methods**:
- `insert_variant()`: Core insertion logic
- `_find_optimal_position()`: Position selection algorithm
- `_update_quality_scores()`: Maintains quality score consistency

**Insertion Strategy**:
- Random position selection within read
- Quality score interpolation at junctions
- Preserves read structure integrity

### 5. CLI Interface (`cli.py`)

**Purpose**: Command-line interface using Click framework.

**Commands**:
- `simulate`: Main simulation command
- `validate-config`: Configuration file validation

**Key Parameters**:
- `--max-reads-per-window`: Clustering prevention (default=1)
- `--min-mapq`: Quality threshold (default=20)
- `--min-read-length`: Length filter (default=10000)
- `--max-soft-clip-ratio`: Soft-clip filter (default=0.1)

## Data Flow

```
1. Input BAM → ReadSelector
   ↓ (filtered reads with metadata)
2. Config JSON → VariantGenerator
   ↓ (generated variants)
3. VariantRegistry
   ↓ (tracked variants)
4. ReadInserter (combines reads + variants)
   ↓ (modified reads)
5. Output Files (BAM, FASTA, JSON)
```

## Configuration System

### Main Configuration (`variant_config.json`)
```json
{
    "insertions": {
        "random": {
            "sequences": [
                {"length": 250, "gc_content": 0.5, "count": 200}
            ]
        },
        "simple": {
            "sequences": [
                {"unit": "CAG", "units": 50, "count": 100}
            ]
        },
        "predefined": {
            "fasta_file": "predefined/alu_elements.fasta",
            "sequences": [
                {"seq_id": "AluYa5", "count": 50}
            ]
        }
    }
}
```

### Pipeline Configuration (`snakemake_config.yml`)
- Input/output paths
- Resource allocation
- Variant caller parameters
- Analysis settings

## Output Files

### 1. Modified Reads (`*.modified_reads.fasta`)
- FASTA format with full modified sequences
- Headers preserve original read names

### 2. Insertion Records (`*.insertions.json`)
```json
{
    "read_name": "read_001",
    "variant_id": "uuid",
    "insertion_position": 5000,
    "original_chr": "chr1",
    "original_pos": 1000000
}
```

### 3. Variant Registry (`*.registry.json`)
- Complete variant sequences and metadata
- Enables variant recovery from calls

### 4. Statistics (`*.statistics.json`)
- Read selection metrics
- Variant generation summary
- Processing performance data

## Extension Points

### Adding New Variant Types
1. Extend `VariantGenerator` with new generation method
2. Update configuration schema
3. Add corresponding tests

### Custom Read Filters
1. Add filter method to `ReadSelector`
2. Expose via CLI parameter
3. Update documentation

### Output Format Extensions
1. Add writer method to output module
2. Integrate into main workflow
3. Update file manifest

## Performance Considerations

### Memory Usage
- Streaming BAM processing (no full load)
- Batch processing for large selections
- Efficient variant storage

### Optimization Strategies
- Indexed BAM access for random sampling
- Preprocessed chromosome lengths
- Cached GC calculations

### Scalability
- Supports millions of reads
- Handles genome-wide distribution
- Parallel-friendly design

## Error Handling

### Validation Layers
1. Configuration validation at startup
2. BAM file integrity checks
3. Read filter warnings
4. Output verification

### Recovery Mechanisms
- Partial output on failure
- Detailed error logging
- State preservation for debugging

## Integration with Analysis Pipeline

### Snakemake Workflow
1. Nova simulation
2. Read alignment (minimap2)
3. Variant calling (Sniffles2)
4. Result analysis (custom scripts)

### Analysis Scripts
- `analyze_vcf_results.py`: Detection validation
- `analyze_false_positives.py`: FP categorization
- `compare_nova_results.py`: Cross-run comparison

## Best Practices

### Module Independence
- Clear interfaces between modules
- Minimal coupling
- Testable components

### Data Validation
- Type hints throughout
- Runtime validation
- Comprehensive error messages

### Logging Strategy
- Module-level loggers
- Configurable verbosity
- Structured log format

## Future Architecture Considerations

### Planned Enhancements
1. Deletion/duplication support
2. Complex variant patterns
3. Multi-sample simulation
4. Real-time progress tracking

### Modular Extensions
- Variant effect prediction
- Read technology simulation
- Quality score modeling
- Coverage-aware selection

## Related Code References

### GenomicInterval Implementation
For VCF annotation with genomic intervals, see the GenomicInterval and BEDProcessor classes in:
`/Users/stacy/sudmant/vesper/src/vesper/`

This implementation provides more comprehensive genomic interval operations including proximal and overlapping BED annotations, which may be useful for future Nova enhancements requiring complex genomic coordinate operations beyond the current region-targeting features.