"""Implements ``snakemake`` rule `aggregate_variants_by_sample`."""


import pandas as pd


# get variables from snakemake
bam_tsvs = snakemake.input.bams
sra_tsvs = snakemake.input.sras
output_csv = snakemake.output.csv
bam_runs = snakemake.params.bam_runs
sra_suffixes = snakemake.params.sras
name = snakemake.wildcards.name

print(f"Analyzing variants for {name}")

dfs = []
assert len(bam_tsvs) == len(bam_runs)
print(f"Reading from {len(bam_tsvs)} BAM-derived runs")
for tsv, bam_run in zip(bam_tsvs, bam_runs):
    print(f"Reading {tsv} for BAM-derived run {bam_run}")
    dfs.append(
        pd.read_csv(tsv, sep='\t')
        .assign(source='BAM',
                run=bam_run,
                )
        )

assert len(sra_tsvs) == len(sra_suffixes)
print(f"Reading from {len(sra_tsvs)} SRA-derived files")
for tsv, sra_suffix in zip(sra_tsvs, sra_suffixes):
    print(f"Reading {tsv} for SRA-derived file suffixed {sra_suffix}")
    dfs.append(
        pd.read_csv(tsv, sep='\t')
        .assign(source='SRA',
                run=sra_suffixes[sra_suffix],
                )
        )

assert all(dfs[0].columns.tolist() == df.columns.tolist()
           for df in dfs)
df = pd.concat(dfs, ignore_index=True).assign(name=name)
assert len(df) == len(df.drop_duplicates())

print(f"Writing results to {output_csv}")
df.to_csv(output_csv, index=False)
