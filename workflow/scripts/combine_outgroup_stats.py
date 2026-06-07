"""
Combine per-population outgroup stats into one CSV.
"""
import numpy as np
import pandas as pd

try:
    in_files = snakemake.input
    populations = list(snakemake.params.populations)
    out_file = snakemake.output[0]
except NameError:
    # minimal local test
    in_files = [
        "results/vcf/Gorilla_gorilla/Gorilla_beringei.outgroup.stats.csv",
        "results/vcf/Homo_sapiens/Homo_sapiens.outgroup.stats.csv"
    ]
    populations = ["Gorilla_gorilla", "Homo_sapiens"]
    out_file = "scratch/combined.csv"

stats = pd.DataFrame(columns=[
    'population',
    'match',
    'mismatch',
    'no_data',
    'mismatch_fraction'
])
for i, p in enumerate(populations):
    df = pd.read_csv(in_files[i], index_col=0)  # df index = outgroup name
    row = {
        'population': p,
        'match': np.array([r.match for r in df.itertuples()]),
        'mismatch': np.array([r.mismatch for r in df.itertuples()]),
        'no_data': np.array([r.no_data for r in df.itertuples()]),
        'mismatch_fraction': np.array([r.mismatch_fraction for r in df.itertuples()]),
    }

    stats.loc[len(stats)] = row

stats.to_csv(out_file, index=False)

