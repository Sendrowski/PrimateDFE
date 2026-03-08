"""
Plot SFS
"""

import fastdfe as fd
import pandas as pd

try:
    testing = False
    file_neut = snakemake.input.neutral
    file_sel = snakemake.input.selected
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    file_neut = "resources/Castellano/human_4fold_all_sfs.txt"
    file_sel = "resources/Castellano/human_0fold_all_sfs.txt"
    out = f"scratch/sfs.png"

spectra = fd.Spectra(dict(
    neutral=pd.read_csv(file_neut, header=None).iloc[0],
    selected=pd.read_csv(file_sel, header=None).iloc[0]
))

spectra.plot(show=testing, file=out)
