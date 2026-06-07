"""
Compute the site frequency spectrum (SFS) for a VCF file.
"""

import fastdfe as fd
import numpy as np

try:
    testing = False
    vcf = snakemake.input.vcf
    n = snakemake.params["n"]
    ingroups = snakemake.params["ingroups"]
    max_sites = snakemake.params.get("max_sites", np.inf)
    skip_non_polarized = snakemake.params.get("skip_non_polarized", True)
    n_target_sites = snakemake.params.get("n_target_sites", 0)
    fold = snakemake.params.get("fold", False)
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf = "results/vcf/Homo_sapiens/Pongo_abelii.degeneracy.polarized.vcf.gz"
    n = 10
    ingroups = ['PD_0262', 'PD_0263', 'SAMN01920542', 'SAMN01920543', 'SAMN01920544', 'SAMN01920545', 'SAMN01920546']
    max_sites = 100
    skip_non_polarized = False
    n_target_sites = 10 ** 6
    fold = False
    out = f"scratch/chr22.sfs.{n}.vcf.gz"

p = fd.Parser(
    vcf=vcf,
    n=n,
    include_samples=ingroups,
    max_sites=max_sites,
    skip_non_polarized=skip_non_polarized
)

sfs = p.parse()

if fold:
    sfs = sfs.fold()

# we only have one category, so simply take remaining sites to be monomorphic
sfs.data.loc[0, 'all'] = n_target_sites - sfs.all.n_polymorphic

sfs.to_file(out)
