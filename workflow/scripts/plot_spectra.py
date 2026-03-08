"""
Plot SFS
"""

import fastdfe as fd

try:
    testing = False
    file = snakemake.input[0]
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    file = "results/sfs/Homo_sapiens/Homo_sapiens/sfs.8.csv"
    out = f"scratch/sfs.png"

fd.Spectra.from_file(file).plot(show=testing, file=out)
