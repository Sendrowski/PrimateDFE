"""
Joint call counts with metadata.
"""

import pandas as pd

try:
    testing = False
    metadata = snakemake.input.metadata
    call_counts = snakemake.input.call_counts
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    metadata = "resources/samples.csv"
    call_counts = "results/vcf/primates.subset.1000000.call_counts.csv"
    out = "scratch/joined_call_counts.csv"

metadata = pd.read_csv(metadata)
call_counts = pd.read_csv(call_counts, header=None, names=["ID", "call_count"], sep=" ")

joined = call_counts.merge(metadata, on="ID")

joined.to_csv(out, index=False)
