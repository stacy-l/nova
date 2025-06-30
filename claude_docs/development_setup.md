# Development Environment Setup

## Overview

This guide provides detailed instructions for setting up a development environment for nova. We support two installation methods: conda (recommended) and pip with virtual environments.

## Prerequisites

- Python 3.8 or higher
- Git for cloning the repository
- Conda (Miniforge/Mambaforge recommended) OR pip/venv
- ~2GB disk space for dependencies

## Primary Installation Method (Conda)

### Step 1: Clone the Repository
```bash
git clone https://github.com/yourusername/nova.git
cd nova
```

### Step 2: Create Conda Environment
```bash
conda env create --file environment.yml
```

This creates a `nova` environment with all bioinformatics dependencies including:
- pysam for BAM manipulation
- biopython for sequence handling
- samtools, bedtools, minimap2 binaries
- All Python package dependencies

### Step 3: Activate Environment
```bash
conda activate nova
```

### Step 4: Install Nova in Development Mode
```bash
# Using uv (faster, recommended if available)
uv pip install -e .

# Or using standard pip
pip install -e .
```

### Step 5: Verify Installation
```bash
# Check Python path (should show nova conda env)
which python
# Expected: /Users/[username]/miniforge3/envs/nova/bin/python

# Test nova CLI
nova --help
# Should show: simulate, validate-config commands
```

## Alternative Installation (pip/venv)

### Step 1: Create Virtual Environment
```bash
cd nova
python -m venv nova-test-env
```

### Step 2: Activate Virtual Environment
```bash
# On macOS/Linux
source nova-test-env/bin/activate

# On Windows
nova-test-env\Scripts\activate
```

### Step 3: Install Dependencies
```bash
pip install -r requirements.txt
pip install -e .
```

⚠️ **Note**: This method requires system installations of samtools, bedtools, and minimap2.

### Step 4: Verify Installation
```bash
nova --help
```

## Environment Management

### For Development Work
Always ensure the correct environment is active:
```bash
# For conda users
conda activate nova

# For venv users
source nova-test-env/bin/activate
```

### For Testing
```bash
# Activate test environment
source nova-test-env/bin/activate

# Run tests
python -m pytest tests/ -v
```

### Environment Variables
Nova supports these optional environment variables:
- `NOVA_LOG_LEVEL`: Set logging verbosity (DEBUG, INFO, WARNING, ERROR)
- `NOVA_TEMP_DIR`: Override temporary file location

## Troubleshooting

### Common Issues

#### 1. "nova: command not found"
- Ensure environment is activated
- Verify installation with `pip show nova`
- Try `python -m nova.cli` as alternative

#### 2. Import errors for pysam
- Conda: Reinstall with `conda install pysam`
- Pip: May need system dependencies: `apt-get install libbz2-dev liblzma-dev`

#### 3. BAM index errors
- Ensure BAM files have accompanying .bai index files
- Create with: `samtools index your_file.bam`

#### 4. Memory issues with large BAMs
- Use `--max-reads` to limit selection
- Ensure sufficient system RAM (8GB+ recommended)

### Platform-Specific Notes

#### macOS
- Xcode Command Line Tools required: `xcode-select --install`
- Apple Silicon (M1/M2): Use Miniforge for native ARM64 support

#### Linux
- Ensure development headers installed:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install build-essential libbz2-dev liblzma-dev

  # CentOS/RHEL
  sudo yum groupinstall "Development Tools"
  sudo yum install bzip2-devel xz-devel
  ```

#### Windows
- Use WSL2 for best compatibility
- Native Windows requires Visual Studio Build Tools

## Development Workflow

### 1. Make Code Changes
Edit source files in `src/nova/`

### 2. Test Changes
```bash
# Run specific test
python -m pytest tests/test_variant_generator.py -v

# Run with coverage
python -m pytest tests/ --cov=nova
```

### 3. Lint and Format
```bash
# If configured
black src/nova/
flake8 src/nova/
```

### 4. Build Documentation
```bash
# If using Sphinx
cd docs
make html
```

## IDE Configuration

### VS Code
Recommended extensions:
- Python
- Pylance
- Python Test Explorer

Settings (`.vscode/settings.json`):
```json
{
    "python.defaultInterpreterPath": "~/miniforge3/envs/nova/bin/python",
    "python.testing.pytestEnabled": true,
    "python.testing.pytestArgs": ["tests"]
}
```

### PyCharm
1. Set Project Interpreter to nova conda environment
2. Mark `src` as Sources Root
3. Configure pytest as test runner

## Next Steps

- Read the [Architecture Guide](architecture.md) to understand code structure
- Review [Testing Philosophy](testing_philosophy.md) for test practices
- Check [Project Overview](project_overview.md) for scientific context
- See [CLAUDE.md](../CLAUDE.md) for quick command reference