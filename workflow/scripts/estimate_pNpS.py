"""
Estimate pN/pS from SFS.
"""

import fastdfe as fd
import numpy as np

try:
    testing = False
    sfs_file = snakemake.input[0]
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    sfs_file = "results/sfs/Homo_sapiens/Homo_sapiens/sfs.unfolded.8.csv"
    out = "scratch/Ne.txt"

sfs = fd.Spectra.from_file(sfs_file)


def get_pi(sfs: fd.Spectrum) -> float:
    """
    Nucleotide diversity from an unfolded SFS.
    """
    i = np.arange(sfs.n + 1)

    w = 2 * i * (sfs.n - i) / (sfs.n * (sfs.n - 1))

    return np.sum(w * sfs.data) / sfs.n_sites


pS = get_pi(sfs['neutral'])
pN = get_pi(sfs['selected'])

pnps = pN / pS

with open(out, "w") as h:
    h.write(str(pnps))
