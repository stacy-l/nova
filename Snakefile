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
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_statistics.json",
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_simulation.jl",
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_base.jl"

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

rule map_reads:
    input:
        read_fasta = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_modified_reads.fasta"
    output: 
        bam = temp(f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_modified_reads.temp.bam")
    params:
        refgenome = config['reference_genome'],
        readgroup = f"@RG\\tID:nova_sim_{config['output_prefix']}",
        minQ = config['min_mapq']
    conda: 
        "environment.yml"
    threads: 4
    shell: 
        """
        minimap2 {params.refgenome} {input.read_fasta} -t {threads} -ax map-hifi -Y -y -L --eqx --cs --MD -R '{params.readgroup}' | samtools view -q {params.minQ} -bT {params.refgenome} -o {output.bam}
        """

rule sort_reads:
    input:
        bam = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_modified_reads.temp.bam"
    output:
        bam = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_modified_reads.bam"
    shell:
        "samtools sort -@ {threads} {input.bam} -o {output.bam}"

rule index_bam:
    input:
        bam = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_{{fn}}.bam"
    output:
        bai = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_{{fn}}.bam.bai"
    conda: 
        "environment.yml"
    threads: 4
    shell:
        "samtools index -@ {threads} {input.bam} {output.bai}"

rule merge_sim_base_bams:
    input:
        sim_bam = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_modified_reads.bam",
        base_bam = config["bam_file"]
    output:
        merged_bam = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_simulation.bam"
    conda: 
        "environment.yml"
    threads: 4
    shell:
        "samtools merge -@ {threads} {output.merged_bam} {input.sim_bam} {input.base_bam}"

rule call_variants_base:
    input:
        bam = config['bam_file'],
        index = config['bam_file'] + '.bai'
    output:
        vcf=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_base.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_base.vcf.gz.tbi"
    conda: 
        "environment.yml"
    threads:
        4
    params:
        refgenome = config['reference_genome'],
        repeats = config['repeats'],
        minsupport = config['minsupport'],
        mapq = config['mapq'],
        mosaic_af_min = config['mosaic_af_min'],
        mosaic_af_max = config['mosaic_af_max'],
        mosaic_qc_strand = config['mosaic_qc_strand']
    log:
        f"{OUTPUT_DIR}/logs/{OUTPUT_PREFIX}_base.qc_all.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --reference {params.refgenome} \
        --tandem-repeats {params.repeats} \
        --threads {threads} \
        --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min {params.mosaic_af_min} \
        --mosaic-af-max {params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} \
        --no-qc &> {log}
        """

use rule call_variants_base as call_variants_sim with:
    input:
        bam = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_simulation.bam"
    output:
        vcf=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_simulation.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_simulation.vcf.gz.tbi"
    log:
        f"{OUTPUT_DIR}/logs/{OUTPUT_PREFIX}_simulation.qc_all.log"

rule simulation_benchmark:
    input:
        query = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_simulation.vcf.gz",
        query_index = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_simulation.vcf.gz.tbi",
        benchmark = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_base.vcf.gz",
        benchmark_index = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_base.vcf.gz.tbi"
    output:
        expand(f"{OUTPUT_DIR}/{{outfiles}}",
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: 
        "environment.yml"
    threads: 1
    params:
        refgenome = config['reference_genome'],
        outdir=config["output_dir"],
    shell:
        """
        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification?)
        truvari bench \
        -f {params.refgenome} \
        -b {input.base} \
        -c {input.query} \
        -o {params.outdir}/bench \
        -r 1000 \
        --dup-to-ins
        """

rule vcf2df:
    input:
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_{{fn}}.vcf.gz"
    output:
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_{{fn}}.jl"
    conda: 
        "environment.yml"
    threads: 1
    shell:
        """
        truvari vcf2df --info --format {input} {output}
        """