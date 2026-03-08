"""
Remove monomorphic sites from site frequency spectra.
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
    out = f"scratch/sfs_no_monomorphic.csv"

spectra = fd.Spectra.from_file(file)

# remove monomorphic sites
spectra.data.iloc[0] = 0
spectra.data.iloc[spectra.n] = 0

spectra.to_file(out)