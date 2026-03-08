"""
Compute the unfolded site frequency spectrum (SFS) for a VCF file.
"""

import fastdfe as fd
import numpy as np

try:
    testing = False
    vcf = snakemake.input[0]
    n = snakemake.params["n"]
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf = "results/vcf/hg38/ingroup/Gorilla_gorilla_gorilla/chr/chr22.vcf.gz"
    n = 1000
    out = f"scratch/sfs.{n}.csv"

p = fd.Parser(
    vcf=vcf,
    n=n,
    skip_non_polarized=False
)

sfs = p.parse()

if len(sfs.types) == 0:
    sfs = fd.Spectra(dict(all=np.zeros(n + 1)))

sfs.to_file(out)

pass