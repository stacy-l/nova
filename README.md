# `nova`: *de novo* insertion simulator

## Overview

`nova` simulates de novo variant insertions in long-read sequencing data. 

The simulator can generate random, simple repeat, and predefined insertion sequences and insert them into reads at ultra-low frequencies typical of de novo structural variants.

## Installation

### `conda` (recommended)
We recommend using `mamba` or `conda` for (almost) one-step installation. 

```bash
# we recommend using mamba for faster environment resolution
mamba env create -f environment.yml --subdir osx-64
conda activate nova
conda config --env --set subdir osx-64 

uv pip install -e. # switch to conda later, but uv recognizes and respects activate conda env
```

### Standalone install
`nova` can be installed using `uv` or `pip`. We recommend installing it into a venv managed by `uv`, like so:

```bash
uv venv nova-venv
source nova-venv/bin/activate
uv pip install -e. # atm this is dev only, switch to nova later if up on pypi
```

The requirements for `nova` are detailed in both `setup.py` and `requirements.txt`. However, there are several dependencies that aren't available via PyPi, hence why we recommend `conda`/`mamba`. 

If you want to avoid that path entirely, you'll need to install additional dependencies to use the `snakemake` pipeline:
```
bedtools==2.31.1
samtools==1.22
bcftools==1.22
minimap2==2.30
sniffles==2.6.1
truvari==5.3.0
snakemake-minimal==9.6.2
```

## Quick Start

### Read simulation
`nova` has one main program (`nova simulate`), which samples reads from an input BAM file and modifies them with insertion variants according to a configuration JSON file.

```bash
nova simulate input.bam config.json \
--random-seed 783 \
-o output_dir \
-p nova
```

For more options and info, read the [simulate](https://github.com/stacy-l/nova/docs/simulate) and [configuration](https://github.com/stacy-l/nova/docs/configuration.md) docs.

### Analysis

`nova` is a variant simulator and doesn't perform variant calling and simulation analysis on its own. There are a number of accompanying scripts that you can use to assess your simulations, described in the [analysis docs](https://github.com/stacy-l/nova/docs/analysis.md) and run by default in the analysis pipeline.


## Analysis pipeline

We recommend using the `environment.yml` file to install all requirements to run the full pipeline.

1. Prepare your BAM file and variant configuration JSON file specifying desired variants
2. Prepare the config file (template in `examples/snakemake_config.yml`).
    - **Note**: You will need to self-supply the path to the reference genome and repeat file (for `sniffles2`).
3. Dry run the pipeline, specifying the path to your edited config file:
```bash
snakemake --use-conda --configfile {edited_config.yml} -np
```
4. Once the dry run is working, run the pipeline with resource-appropriate core allocation:
```bash
snakemake --use-conda --configfile {edited_config.yml} --cores {cores}
```

## License

(tbd?)