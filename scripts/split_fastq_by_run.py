"""Implements ``snakemake`` rule `split_fastq_by_machine`.

Reads entire FASTQ into memory, so not a good idea to run
on really large files.

"""


import collections
import gzip
import os
import re

import pysam


# get variables from `snakemake`
fastq = snakemake.input.fastq
subdir = snakemake.output.subdir

os.makedirs(subdir, exist_ok=True)

# regex for Illumina Casava 1.8 header:
# https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
regex = re.compile('(?P<run>[\w\-]+:\d+:[\w\-]+:\d+):\d+:\d+:\d+')

print(f"Reading reads from {fastq}")
with pysam.FastxFile(fastq) as f_in:
    reads = collections.defaultdict(list)
    for read in f_in:
        m = regex.match(read.name)
        if not m:
            raise ValueError(f"cannot match {read.name}")
        run = m.group('run')
        reads[run].append(read)

print(f"Found reads from {len(reads)} Illumina runs.")
for run, run_reads in reads.items():
    outfile = os.path.join(subdir, f"{run}.fastq.gz")
    print(f"Writing {len(run_reads)} reads for {run} to {outfile}")
    with gzip.open(outfile, 'wt') as f_out:
        for read in run_reads:
            f_out.write(str(read))
            f_out.write('\n')
