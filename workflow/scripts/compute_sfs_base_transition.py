"""
Compute the site frequency spectrum (SFS) for a VCF file.
"""

import fastdfe as fd
import numpy as np
import pandas as pd

try:
    testing = False
    vcf = snakemake.input.vcf
    gff = snakemake.input.gff
    fasta = snakemake.input.fasta
    n = snakemake.params["n"]
    ingroups = snakemake.params["ingroups"]
    n_target_sites = snakemake.params["n_target_sites"]
    skip_non_polarized = snakemake.params.get("skip_non_polarized", True)
    fold = snakemake.params.get("fold", False)
    filter_bgc = snakemake.params.get("filter_bgc", False)
    n_samples = snakemake.params.get("n_samples", 1000000)
    max_sites = snakemake.params.get("max_sites", np.inf)
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf = "results/vcf/Homo_sapiens/Homo_sapiens.degeneracy.polarized.vcf.gz"
    gff = "resources/gff/Homo_sapiens.renamed.gff.gz"
    fasta = "resources/ref/Homo_sapiens.fasta"
    n = 8
    samples = pd.read_csv(f"resources/metadata/Homo_individuals.csv", sep="\t")
    ingroups = samples[samples["BAM_FOLDER"].str.contains('Homo_sapiens')].GVCF_ID.tolist()
    n_target_sites = 2 * 10 ** 7
    skip_non_polarized = False
    fold = False
    filter_bgc = False
    n_samples = 1000
    max_sites = 30000
    out = "scratch/sfs.png"

p = fd.Parser(
    vcf=vcf,
    n=n,
    gff=gff,
    fasta=fasta,
    max_sites=max_sites,
    include_samples=ingroups,
    target_site_counter=fd.TargetSiteCounter(
        n_target_sites=n_target_sites,
        n_samples=n_samples
    ),
    stratifications=[
        fd.DegeneracyStratification(),
        fd.BaseTransitionStratification()
    ],
    annotations=[
        fd.DegeneracyAnnotation()
    ],
    filtrations=[fd.BiasedGCConversionFiltration()] if filter_bgc else [],
    polarize_probabilistically=True,
    skip_non_polarized=skip_non_polarized
)

sfs = p.parse()

if fold:
    sfs = sfs.fold()

sfs.to_file(out)

print(sfs.data)

if testing:
    sfs.plot()

sfs_bgc_con = fd.Spectra(dict(
    neutral=sfs['neutral.A>T'] + sfs['neutral.T>A'] + sfs['neutral.C>G'] + sfs['neutral.G>C'],
    selected=sfs['selected.A>T'] + sfs['selected.T>A'] + sfs['selected.C>G'] + sfs['selected.G>C']
))

sfs_bgc_div = fd.Spectra(dict(
    neutral=sfs['neutral.A>G'] + sfs['neutral.A>C'] + sfs['neutral.T>G'] + sfs['neutral.T>C'] + sfs['neutral.C>A'] + sfs['neutral.C>T'] + sfs['neutral.G>A'] + sfs['neutral.G>T'],
    selected=sfs['selected.A>G'] + sfs['selected.A>C'] + sfs['selected.T>G'] + sfs['selected.T>C'] + sfs['selected.C>A'] + sfs['selected.C>T'] + sfs['selected.G>A'] + sfs['selected.G>T']
))

sfs_all = fd.Spectra(dict(
    neutral=sfs['neutral.*'].all,
    selected=sfs['selected.*'].all
))

np.testing.assert_array_almost_equal(sfs_bgc_con.data + sfs_bgc_div.data, sfs_all.data)

pass