# Nova: De Novo Variant Insertion Simulator

A Python simulation framework to evaluate structural variant detection tools' ability to detect ultra-low frequency de novo insertions.

## Overview

Nova simulates de novo variant insertions in long-read sequencing data to test the sensitivity and specificity of structural variant detection tools like `vesper`. The simulator can generate random, simple repeat, and predefined insertion sequences and insert them into reads at ultra-low frequencies typical of de novo structural variants.

## Installation

```bash
conda env create --file environment.yml
conda activate nova
uv pip install -e .
```

## Usage

### Basic Simulation

```bash
nova simulate input.bam config.json -o output_dir
```

### Configuration File

Create a JSON configuration file specifying the types and numbers of variants to simulate:

```json
{
  "random": {
    "n": 50,
    "length": 100,
    "gc_content": 0.5
  },
  "simple": {
    "n": 50,
    "repeat": "CAG",
    "units": 40
  },
  "predefined": {
    "Alu": {
      "fasta": "dfam_AluY_homininae.fasta",
      "spec": {
        "AluYa5": 30,
        "AluYb8": 20
      }
    }
  }
}
```

### Command Line Options

- `--output-dir, -o`: Output directory for results (default: current directory)
- `--output-prefix, -p`: Prefix for output files (default: 'nova_sim')
- `--min-mapq`: Minimum mapping quality (default: 20)
- `--max-soft-clip-ratio`: Maximum soft clipping ratio (default: 0.1)
- `--min-read-length`: Minimum read length (default: 10000)
- `--max-read-length`: Maximum read length (default: 20000)
- `--min-distance-from-ends`: Minimum distance from read ends for insertion (default: 1000)
- `--random-seed`: Random seed for reproducibility

### Utility Commands

#### Validate Configuration
```bash
nova validate-config config.json
```

## Output Files

Nova generates several output files:

1. **`*.insertions.json`**: Detailed insertion records with positions and metadata
2. **`*.modified_reads.fasta`**: FASTA file with modified reads containing insertions
3. **`*.registry.json`**: Registry of all generated insertion sequences

## Features

### Read Selection
- Filters reads by length and soft-clip ratio
- Future implementation:
    - Stratifies reads by genomic complexity (repeats, segmental duplications, unique regions)
    - Supports inclusion/exclusion of specific genomic regions

### Insertion Generation
- **Random insertions**: Sequences with specified length and GC content
- **Simple repeats**: Tandem repeats with specified unit and count
- **Predefined sequences**: From FASTA files (e.g., Alu elements)

### Read Modification
- Random positioning within reads (avoiding ends)
- Preserves original read metadata
- Generates comprehensive insertion records

## Module Architecture

- `read_selector.py`: BAM file processing and read selection
- `variant_generator.py`: Insertion sequence generation
- `variant_registry.py`: Sequence management with unique IDs
- `read_inserter.py`: Read modification and insertion
- `cli.py`: Command-line interface

## Testing

### Setup Test Environment

Tests use pytest and require a Python virtual environment:

```bash
# Using uv (recommended)
uv venv nova-test-env
source nova-test-env/bin/activate
uv pip install -r requirements.txt
uv pip install -e .

# Or using pip
python -m venv nova-test-env
source nova-test-env/bin/activate
pip install -r requirements.txt
pip install -e .
```

### Run the Test Suite

```bash
# Activate your virtual environment first
source nova-test-env/bin/activate

# Run all tests
python -m pytest tests/ -v

# Run with coverage
python -m pytest tests/ --cov=nova

# Run specific test modules
python -m pytest tests/test_variant_generator.py
```

Most tests use mocks and don't require real BAM files. See `tests/test_data/README.md` for information about providing BAM files for full integration testing.

## Example Workflow

1. Prepare your BAM file and optional BED annotation files
2. Use `bedtools`, `samtools`, or the Snakemake pipeline to subsample reads
3. Create a configuration JSON file specifying desired variants
4. Run `nova`:
   ```bash
   nova simulate reads.bam config.json -o output_dir
   ```
5. Analyze the output files to evaluate your variant detection tool

## Requirements

- Python 3.8+
- pysam
- pybedtools
- numpy
- biopython
- click

## License

MIT License