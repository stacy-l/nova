# Nova Analysis Scripts

This directory contains analysis scripts for processing and evaluating nova simulation results.

## Scripts Overview

### 1. `analyze_vcf_results.py`
**Purpose**: Analyzes variant detection results from VCF files to calculate performance metrics.

**Key Features**:
- Matches detected variants with ground truth insertions
- Calculates true positive rates and precision
- Categorizes variants by insertion type
- Identifies false positives and false negatives

**Usage**:
```bash
python scripts/analyze_vcf_results.py \
    --vcf results/nova_simulation.vcf.gz \
    --insertions results/nova_insertions.json \
    --registry results/nova_registry.json \
    --output results/nova_variant_analysis
```

**Output Files**:
- `*_analysis.json`: Detailed analysis results
- `*_analysis_tabular.csv`: CSV format for easy analysis

### 2. `analyze_false_positives.py`
**Purpose**: Deep analysis of false positive variants to understand their origins.

**Key Features**:
- Categorizes false positives by read composition
- Detects genomic clustering patterns
- Analyzes CIGAR string complexity
- Provides detailed breakdowns by category

**Usage**:
```bash
python scripts/analyze_false_positives.py \
    results/nova_variant_analysis_tabular.csv \
    results/false_positives
```

**Categories Identified**:
- **multiple_nova_only**: Multiple nova reads clustered
- **mixed_nova_minority**: Nova reads < 50% of support
- **mixed_nova_majority**: Nova reads ≥ 50% of support

### 3. `compare_nova_results.py`
**Purpose**: Compares results between different nova simulation runs.

**Key Features**:
- Matches variants between runs by position and type
- Identifies common and unique variants
- Calculates concordance metrics
- Supports flexible matching tolerances

**Usage**:
```bash
python scripts/compare_nova_results.py \
    run1/output_dir \
    run2/output_dir \
    --output comparison_results \
    --pos-tolerance 50 \
    --size-tolerance 0.2
```

**Output**:
- Matched variant pairs
- Run-specific variants
- Concordance statistics

## Common Workflows

### 1. Full Analysis Pipeline
```bash
# Run nova simulation
nova simulate input.bam config.json -o results/

# Align reads (if using Snakemake, this is automated)
minimap2 -ax map-hifi ref.fa results/nova_modified_reads.fasta | \
    samtools sort -o results/nova_aligned.bam

# Call variants
sniffles -i results/nova_aligned.bam -v results/variants.vcf.gz

# Analyze results
python scripts/analyze_vcf_results.py \
    --vcf results/variants.vcf.gz \
    --insertions results/nova_insertions.json \
    --registry results/nova_registry.json \
    --output results/analysis

# Investigate false positives
python scripts/analyze_false_positives.py \
    results/analysis_tabular.csv \
    results/fp_analysis
```

### 2. Parameter Comparison
```bash
# Run with different window limits
nova simulate input.bam config.json -o window1/ --max-reads-per-window 1
nova simulate input.bam config.json -o window3/ --max-reads-per-window 3

# Compare results
python scripts/compare_nova_results.py window1/ window3/ -o comparison/
```

## Output Format Details

### Variant Analysis JSON
```json
{
    "summary": {
        "total_variants_detected": 1086,
        "total_insertions_simulated": 1000,
        "true_positives": 979,
        "false_positives": 107,
        "false_negatives": 21
    },
    "by_type": {
        "random": {"detected": 385, "expected": 400},
        "simple": {"detected": 195, "expected": 200},
        "AluYa5": {"detected": 98, "expected": 100}
    }
}
```

### False Positive Analysis JSON
```json
{
    "metadata": {
        "total_false_positives": 107,
        "analysis_date": "2024-01-15"
    },
    "category_breakdown": {
        "multiple_nova_only": 0,
        "mixed_nova_minority": 45,
        "mixed_nova_majority": 62
    },
    "pattern_statistics": {
        "variants_with_genomic_clustering": 0,
        "avg_support_reads": 5.2
    }
}
```

## Script Parameters

### Common Parameters
- `--threads`: Number of parallel threads (default: 1)
- `--verbose`: Enable detailed logging
- `--output-prefix`: Custom prefix for output files

### Matching Tolerances
- `--pos-tolerance`: Maximum position difference in bp (default: 50)
- `--size-tolerance`: Maximum size difference as fraction (default: 0.2)
- `--min-reciprocal-overlap`: For complex variant matching (default: 0.8)

## Performance Considerations

### Large VCF Files
- Scripts use streaming parsing for memory efficiency
- Consider filtering VCF first for very large files
- Use `--threads` for parallel processing where available

### Many Simulation Runs
- Comparison script can handle multiple runs
- Output becomes large with many pairwise comparisons
- Consider hierarchical comparison strategies

## Extending Scripts

### Adding New Metrics
1. Modify the analysis logic in the script
2. Update output schema
3. Add command-line parameter if needed
4. Update this documentation

### Custom Filters
Scripts support custom filtering via modification:
```python
# Example: Add custom quality filter
if variant.qual < args.min_qual:
    continue
```

## Troubleshooting

### Common Issues

1. **"File not found" errors**
   - Verify paths to input files
   - Check file extensions (.gz vs uncompressed)

2. **Memory errors with large files**
   - Use streaming mode if available
   - Process in chunks
   - Increase system memory

3. **No variants matched**
   - Check matching tolerances
   - Verify coordinate systems match
   - Ensure variant types are compatible

### Debug Mode
Add `--verbose` to any script for detailed logging:
```bash
python scripts/analyze_vcf_results.py --verbose ...
```

## Integration with Snakemake

These scripts are integrated into the Snakemake pipeline:
```yaml
rule analyze_variants:
    input:
        vcf = "{sample}.vcf.gz",
        insertions = "{sample}_insertions.json",
        registry = "{sample}_registry.json"
    output:
        "{sample}_analysis.json"
    shell:
        "python scripts/analyze_vcf_results.py ..."
```

## Future Enhancements

- Real-time analysis during simulation
- Interactive visualization outputs
- Machine learning for FP prediction
- Multi-sample analysis support