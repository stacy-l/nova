# Test Data Directory

This directory contains test files for `nova`'s test suite.

## Provided Files

The `tests/test_data` directory contains several sample configuration files (`test_config_*.json`) for test use.

## Required Files for Full Testing

For complete integration testing, you'll need to provide the following files:

1. **test_reads.bam** - A sample BAM file containing long reads (PacBio)
    - At least 100 mapped reads
    - Reads should be between 10-20kb in length
    - Should contain mapped reads with MAPQ scores ≥ 20
    - Include proper CIGAR strings for soft-clip ratio calculation
   
2. **test_reads.bam.bai** - BAM index file (created with `samtools index test_reads.bam`)

Store these files at:
- `tests/test_data/test_reads.bam`
- `tests/test_data/test_reads.bam.bai`

## Running Tests

```bash
# Run all tests (most use mocks, don't require BAM files)
python -m pytest tests/ -v

# Run only unit tests (no BAM files required)
python -m pytest tests/test_variant_registry.py tests/test_variant_generator.py tests/test_read_inserter.py

# Run integration tests (may require BAM files)
python -m pytest tests/test_integration.py

# Run with coverage
python -m pytest tests/ --cov=nova
``