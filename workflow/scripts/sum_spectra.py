"""
Sum up multiple SFS.
"""

import fastdfe as fd
import numpy as np

try:
    testing = False
    files = snakemake.input
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    files = [
        "results/sfs/hg38/chr/chr22.degeneracy/ingroup_pan_troglodytes/outgroups_ref_gorilla_gorilla_gorilla.ref_pongo_abelii/sfs.11.csv",
        "results/sfs/hg38/chr/chr22.degeneracy/ingroup_pan_troglodytes/outgroups_ref_gorilla_gorilla_gorilla.ref_pongo_abelii/sfs.11.csv",
        "results/sfs/hg38/chr/chr22.degeneracy/ingroup_pan_troglodytes/outgroups_ref_gorilla_gorilla_gorilla.ref_pongo_abelii/sfs.11.csv",
    ]
    out = f"scratch/summed_sfs.csv"

fd.Spectrum(np.sum([fd.Spectrum.from_file(file).data for file in files], axis=0)).to_file(out)
