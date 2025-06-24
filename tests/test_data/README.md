# Test Data Directory

This directory contains test files for Nova's test suite.

## Provided Files

1. **test_sequences.fasta** - Sample FASTA file containing Alu and LINE elements for predefined insertion testing
2. **test_config.json** - Sample configuration file for integration testing

## Required Files for Full Testing

For complete integration testing, you'll need to provide the following files:

### BAM Files
Nova requires BAM files for read selection testing. You'll need to provide:

1. **test_reads.bam** - A sample BAM file containing long reads (PacBio/Oxford Nanopore)
   - Reads should be between 10-20kb in length
   - Should contain mapped reads with MAPQ scores ≥ 20
   - Include proper CIGAR strings for soft-clip ratio calculation
   
2. **test_reads.bam.bai** - BAM index file (created with `samtools index test_reads.bam`)

### How to Create Test BAM Files

If you don't have suitable BAM files, you can create minimal test files:

```bash
# Option 1: Use existing BAM files from your data
cp /path/to/your/long_reads.bam tests/test_data/test_reads.bam
samtools index tests/test_data/test_reads.bam

# Option 2: Create a minimal test BAM with samtools
# (This would require reference genome and read data)
```

### Minimal BAM Requirements for Testing

The test BAM file should contain:
- At least 100 mapped reads
- Read lengths between 10,000-20,000 bp
- Reads with MAPQ ≥ 20
- Proper reference chromosome names (e.g., chr1, chr2, etc.)
- Minimal soft clipping (< 10% of read length)

### Alternative: Mock Testing

Most tests use mocked pysam objects to avoid requiring actual BAM files. Only integration tests that specifically test BAM file reading will need real BAM files.

## Test File Locations

When running tests that require these files, place them as:
- `tests/test_data/test_reads.bam`
- `tests/test_data/test_reads.bam.bai`

## Running Tests

### Setup Test Environment

Tests use pytest and require a Python virtual environment. You can set this up using `uv` or `pip`:

```bash
# Using uv (recommended)
uv venv test-env
source test-env/bin/activate
uv pip install -r requirements.txt
uv pip install -e .

# Or using pip
python -m venv test-env
source test-env/bin/activate
pip install -r requirements.txt
pip install -e .
```

### Running Tests

```bash
# Activate your virtual environment first
source test-env/bin/activate

# Run all tests (most use mocks, don't require BAM files)
python -m pytest tests/ -v

# Run only unit tests (no BAM files required)
python -m pytest tests/test_variant_registry.py tests/test_variant_generator.py tests/test_read_inserter.py

# Run integration tests (may require BAM files)
python -m pytest tests/test_integration.py

# Run with coverage
python -m pytest tests/ --cov=nova
```