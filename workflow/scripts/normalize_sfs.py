"""
Normalize SFS
"""

import fastdfe as fd
import numpy as np

try:
    testing = False
    file = snakemake.input[0]
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    file = "results/sfs/hg38/chr/chr22.degeneracy/ingroup_pan_troglodytes/outgroups_ref_gorilla_gorilla_gorilla.ref_pongo_abelii/sfs.11.csv"
    out = f"scratch/folded_sfs.csv"

spectra = fd.Spectra.from_file(file).normalize().to_file(out)