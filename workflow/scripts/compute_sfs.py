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
    vcf = "results/vcf/Macaca_nemestrina/Macaca_nigra.degeneracy.polarized.vcf.gz"
    gff = "resources/gff/Macaca_nemestrina.gff.gz"
    fasta = "resources/ref/Macaca_nemestrina.fasta"
    n = 8
    samples = pd.read_csv(f"resources/metadata/Macaca_individuals.csv", sep="\t")
    ingroups = samples[samples["BAM_FOLDER"].str.contains('Macaca_nemestrina')].GVCF_ID.tolist()
    n_target_sites = 2 * 10 ** 7
    skip_non_polarized = False
    fold = False
    filter_bgc = True
    n_samples = 1000
    max_sites = 10000
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
        fd.DegeneracyStratification()
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

pass