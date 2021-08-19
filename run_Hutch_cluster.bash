#!/bin/bash

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

printf "Running snakemake...\n"
snakemake \
    -j 20 \
    --cluster "sbatch" \
    --latency-wait 60 \
    --use-conda \
    --keep-going \
    --scheduler greedy \
    --rerun-incomplete
printf "Run of snakemake complete.\n"
