"""
Merge Castellano neutral (4-fold) and selected (0-fold) SFS into one file.
"""

import pandas as pd
import fastdfe as fd

try:
    neutral = snakemake.input["neutral"]
    selected = snakemake.input["selected"]
    out = snakemake.output[0]
except NameError:
    # testing
    neutral = "resources/Castellano/human_4fold_all_sfs.txt"
    selected = "resources/Castellano/human_0fold_all_sfs.txt"
    out = "scratch/sfs.csv"

# read input files
neutral_sfs = pd.read_csv(neutral, sep=",", header=None).iloc[0].tolist()
selected_sfs = pd.read_csv(selected, sep=",", header=None).iloc[0].tolist()

sfs = fd.Spectra(dict(
    neutral=neutral_sfs,
    selected=selected_sfs
))

# write to file
sfs.to_file(out)
