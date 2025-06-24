"""
Nova: De novo variant insertion simulator - Snakemake pipeline
"""

import json
from pathlib import Path

# Configuration
configfile: "snakemake_config.yml"

OUTPUT_DIR = config["output_dir"]
OUTPUT_PREFIX = config["output_prefix"]

rule all:
    input:
        # Core simulation outputs
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_insertions.json",
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_modified_reads.fasta",
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_registry.json",
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_statistics.json"

rule simulate_variant_reads:
    input:
        bam=config["bam_file"],
        variant_config=config["variant_config"]
    output:
        registry=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_registry.json",
        insertions=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_insertions.json",
        read_fasta = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_modified_reads.fasta",
        statistics=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_statistics.json"
    params:
        prefix=config["output_prefix"],
        outdir=config["output_dir"],
        seed=config["random_seed"],
        max_soft_clip=config["max_soft_clip_ratio"],
        min_length=config["min_read_length"],
        max_length=config["max_read_length"],
        min_distance_from_ends=config["min_distance_from_ends"],
    conda:
        "environment.yml"
    shell:
        """
        nova simulate \
        {input.bam} \
        {input.variant_config} \
        --random-seed {params.seed} \
        --max-soft-clip-ratio {params.max_soft_clip} \
        --min-read-length {params.min_length} \
        --max-read-length {params.max_length} \
        --min-distance-from-ends {params.min_distance_from_ends} \
        -p {params.prefix} \
        -o {params.outdir}
        """