"""Snakemake file that runs the analysis.

Written by Jesse Bloom."""


import glob
import os


configfile: 'config.yaml'


# suffixes on FASTQ files extracted from SRA files
SRA_FQ_SUFFIXES = {'': 'single',
                   '_1': 'read_1',
                   '_2': 'read_2',
                   }


rule all:
    input:
        expand("results/aggregated_variants/{name}.csv",
               name=config['SRA_bams'])

rule get_genbank_fasta:
    """Get FASTA from Genbank."""
    output: fasta="results/genbank/{genbank}.fasta"
    conda: 'environment.yml'
    shell: "efetch -format fasta -db nuccore -id {wildcards.genbank} > {output.fasta}"

rule minimap2_genome:
    """Build ``minimap2`` reference genome."""
    input: fasta=f"results/genbank/{config['refgenome']}.fasta"
    output: mmi=f"results/genbank/{config['refgenome']}.mmi"
    conda: 'environment.yml'
    shell: "minimap2 -d {output.mmi} {input.fasta}"

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

checkpoint split_fastq_by_run:
    """Split a FASTQ into different files for each Illumina runs."""
    input: fastq=rules.bam_to_fastq.output.fastq
    output: subdir=directory("results/SRA_bams/{name}_by_run")
    conda: 'environment.yml'
    script: 'scripts/split_fastq_by_run.py'

def list_runs(wildcards):
    """Get a list of all Illumina runs for a BAM for a given name."""
    subdir = checkpoints.split_fastq_by_run.get(**wildcards).output.subdir
    return [os.path.basename(f).replace('.fastq.gz', '')
            for f in glob.glob(f"{subdir}/*.fastq.gz")]

rule align_fastq:
    """Align FASTQ file using ``minimap2``."""
    input:
        fastq="results/{path}.fastq.gz",
        mmi=rules.minimap2_genome.output.mmi
    output:
        sam=temp("results/alignments/{path}.sam"),
        unsorted_bam=temp("results/alignments/{path}.bam"),
        bam="results/alignments/{path}_sorted.bam",
    conda: 'environment.yml'
    shell:
        """
        minimap2 -a {input.mmi} {input.fastq} > {output.sam}
        samtools view -b -F 4 -o {output.unsorted_bam} {output.sam}
        samtools sort -o {output.bam} {output.unsorted_bam}
        """

rule ivar:
    """Call variants and consensus using ``ivar``."""
    input:
        bam=rules.align_fastq.output.bam,
        ref=f"results/genbank/{config['refgenome']}.fasta",
    output:
        pileup=temp("results/ivar/{path}.pileup"),
        consensus="results/ivar/{path}.fa",
        variants="results/ivar/{path}.tsv",
    params:
        minq=config['ivar_minq'],
        depth=config['ivar_mindepth'],
        variants_minfreq=config['ivar_variants_minfreq'],
        prefix="results/ivar/{path}"
    conda: 'environment.yml'
    shell:
        """
        samtools mpileup -aa -A -d 0 -B -Q {params.minq} {input.bam} > {output.pileup}
        cat {output.pileup} | ivar consensus \
            -p {params.prefix} \
            -t 0.5 \
            -q {params.minq} \
            -m {params.depth}
        cat {output.pileup} | ivar variants \
            -p {params.prefix} \
            -t {params.variants_minfreq} \
            -q {params.minq} \
            -m {params.depth} \
            -r {input.ref}
        """

rule aggregate_variants_by_sample:
    """Aggregate all the variants for a sample."""
    input:
        bams=lambda wc: [f"results/ivar/SRA_bams/{wc.name}_by_run/{run}.tsv"
                         for run in list_runs(wc)],
        sras=lambda wc: [f"results/ivar/SRA_files/{wc.name}{sra_fq_suffix}.tsv"
                         for sra_fq_suffix in SRA_FQ_SUFFIXES],
    output: csv="results/aggregated_variants/{name}.csv"
    params:
        bam_runs=list_runs,
        sras=SRA_FQ_SUFFIXES
    conda: 'environment.yml'
    script: 'scripts/aggregate_variants_by_sample.py'
