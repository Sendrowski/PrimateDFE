"""
Combine multiple SFS files into a single SFS file.
"""

import fastdfe as fd
import numpy as np
from collections import defaultdict

try:
    testing = False
    files = snakemake.input
    labels = snakemake.params.labels
    categories = snakemake.params.get("categories", ['neutral', 'selected'])
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    files = [
        "results/sfs/Homo_sapiens/Gorilla_gorilla_gorilla/sfs.unfolded.8.csv",
        "results/sfs/Homo_sapiens/Pan_paniscus/sfs.unfolded.8.csv",
        "results/sfs/Homo_sapiens/Pan_troglodytes/sfs.unfolded.8.csv"
    ]
    labels = ['Gorilla_gorilla_gorilla', 'Pan_paniscus', 'Pan_troglodytes']
    categories = ['neutral', 'selected']
    out = f"scratch/combined_sfs.csv"

spectra = defaultdict(dict)
for label, file in zip(labels, files):
    sfs = fd.Spectra.from_file(file)
    for c in categories:
        spectra[c + '.' + label] = sfs[c]

spectra = fd.Spectra(spectra)

spectra.to_file(out)