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
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.insertions.json",
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.modified_reads.fasta",
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.registry.json",

rule subsample_and_filter_reads:
    input:
        bam=config["bam_file"],
        bai=config["bam_file"] + ".bai"
    output:
        bam=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.subsampled_reads.bam"
    params:
        fraction=config["subsample_fraction"],
        seed=config["random_seed"],
        exclude_flags="0x900"  # secondary (0x100) + supplementary (0x800)
    conda:
        "environment.yml"
    params:
        min_mapq=config["min_mapq"],
    threads: 8
    shell:
        """
        samtools view \
        -@ {threads} \
        --subsample {params.fraction} \
        --subsample-seed {params.seed} \
        --exclude-flags {params.exclude_flags} \
        -q {params.min_mapq} \
        {input.bam} -bo {output.bam}
        """

rule index_subsampled_reads:
    input:
        bam=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.subsampled_reads.bam"
    output:
        bai=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.subsampled_reads.bam.bai"
    conda:
        "environment.yml"
    threads: 4
    shell:
        "samtools index \
        -@ {threads} \
        {input.bam} {output.bai}"

rule simulate_variant_reads:
    input:
        bam=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.subsampled_reads.bam",
        variant_config=config["variant_config"]
    output:
        registry=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.registry.json",
        insertions=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.insertions.json",
        read_fasta = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}.modified_reads.fasta"
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
        --input {input.bam} \
        --config {input.variant_config} \
        --seed {params.seed} \
        --max-soft-clip {params.max_soft_clip} \
        --min-length {params.min_length} \
        --max-length {params.max_length} \
        --min-distance-from-ends {params.min_distance_from_ends} \
        -p {params.prefix} \
        -o {params.outdir}
        """