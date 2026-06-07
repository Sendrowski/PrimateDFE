"""
Compute the site frequency spectrum (SFS) for a VCF file.
"""

import fastdfe as fd

try:
    testing = False
    sfs_file = snakemake.input[0]
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    sfs_file = "results/sfs/Homo_sapiens/Pongo_abelii/sfs.unfolded.20.csv"
    out = "scratch/nn_ns.txt"

sfs = fd.Spectra.from_file(sfs_file)
r = sfs['selected'].data[0] / sfs['neutral'].data[0]

with open(out, 'w') as h:
    h.write(str(r))
