"""Snakemake file that runs the analysis.

Written by Jesse Bloom."""


import os


configfile: 'config.yaml'


# suffixes on FASTQ files extracted from SRA files
SRA_FQ_SUFFIXES = ['', '_1', '_2']


rule all:
    input:
        expand("results/SRA_files/{name}{sra_fq_suffix}.fastq.gz",
               sra_fq_suffix=SRA_FQ_SUFFIXES,
               name=config['SRA_accessions']),
        expand("results/SRA_bams/{name}.fastq.gz",
               name=config['SRA_bams']),

rule get_sra_file:
    """Get `*.sra` files from SRA."""
    output: sra_file=protected("results/SRA_files/{name}.sra")
    params: path=lambda wc: config['SRA_accessions'][wc.name]
    conda: 'environment.yml'
    shell: "wget {params.path} -O {output.sra_file}"

rule sra_to_fastq:
    """Convert `*.sra` files to `*.fastq` files."""
    input: sra_file=rules.get_sra_file.output.sra_file
    output:
        fastqs=temp(expand("results/SRA_files/{name}{sra_fq_suffix}.fastq",
                           sra_fq_suffix=SRA_FQ_SUFFIXES,
                           allow_missing=True))
    params: outdir=lambda wc, output: os.path.dirname(output[0])
    conda: 'environment.yml'
    shell: "fasterq-dump -O {params.outdir} {input.sra_file}"

rule gzip_fastq:
    """Compress a fastq file."""
    input: fastq="{base}.fastq"
    output: fastq_gz="{base}.fastq.gz"
    conda: 'environment.yml'
    shell: "gzip {input.fastq}"

rule get_SRA_bam:
    """Get `*.bam` files from SRA."""
    output: sra_bam=protected("results/SRA_bams/{name}.bam")
    params: path=lambda wc: config['SRA_bams'][wc.name]
    conda: 'environment.yml'
    shell: "wget {params.path} -O {output.sra_bam}"

rule bam_to_fastq:
    """Convert `*.bam` files to `*.fastq`."""
    input: bam="results/SRA_bams/{name}.bam"
    output: fastq="results/SRA_bams/{name}.fastq.gz"
    conda: 'environment.yml'
    shell: "samtools bam2fq {input.bam} | gzip > {output.fastq}"
