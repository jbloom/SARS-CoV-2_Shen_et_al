"""Snakemake file that runs the analysis.

Written by Jesse Bloom."""


configfile: 'config.yaml'


rule all:
    input:
        expand("results/SRA_files/{name}.sra",
               name=config['SRA_accessions']),
        expand("results/SRA_bams/{name}.bam",
               name=config['SRA_bams']),

rule get_SRA_file:
    output: sra_file=protected("results/SRA_files/{name}.sra")
    params: accession=lambda wc: config['SRA_accessions'][wc.name]
    conda: 'environment.yml'
    shell:
        """
        echo "For {wildcards.name} getting {params.accession} to {output.sra_file}"
        prefetch {params.accession} -o {output.sra_file}
        """

rule get_SRA_bam:
    output: sra_bam=protected("results/SRA_bams/{name}.bam")
    params: path=lambda wc: config['SRA_bams'][wc.name]
    conda: 'environment.yml'
    shell:
        """
        echo "For {wildcards.name} getting {params.path} to {output.sra_bam}"
        wget {params.path} -O {output.sra_bam}
        """
