"""
Determines the outgroup coverage of a VCF file.
"""
from collections import defaultdict

import cyvcf2
import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    testing = False
    vcf_in = snakemake.input.vcf
    ingroups = snakemake.params.ingroups
    outgroups = snakemake.params.outgroups
    max_sites = snakemake.params.get('max_sites', np.inf)
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf_in = "results/vcf/Homo_sapiens.vcf.gz"
    # vcf_in = "results/vcf/hg38/autosomes/exons.french_samples.norm.trim.french_snps.snps.pass.vcf.gz"
    #samples = pd.read_csv("resources/samples.csv")
    ingroups = pd.read_csv('resources/PDGP_french.txt').iloc[:, 0].tolist()
    outgroups = [
        # 'ref_pongo_abelii',
        #'ref_pan_troglodytes'
        'Gorilla_gorilla_gorilla'
        #'ref_gorilla_gorilla_gorilla',
    ]
    max_sites = np.inf
    out = "scratch/outgroup_coverage.csv"

reader = cyvcf2.VCF(vcf_in)
ingroup_mask = np.isin(reader.samples, ingroups)
outgroup_mask = np.isin(reader.samples, outgroups)

# n_sites = fd.io_handlers.count_sites(vcf_in)

# for each site for which ingroup is called, check if outgroup is called
coverage = defaultdict(int)
divergence = defaultdict(int)
base_counts = defaultdict(int)
n_sites = 0
for variant in tqdm(reader):

    if n_sites >= max_sites:
        break

    n_sites += 1
    try:
        genotypes = np.array(variant.genotypes)
    except ValueError:
        # inhomogenous shape for sex chromosomes
        continue

    # check if ingroup is called
    if np.sum(genotypes[ingroup_mask][:, [0, 1]] != -1) > 0:

        n_called_outgroups = int(np.sum(genotypes[outgroup_mask][:, [0, 1]] != -1))
        coverage[n_called_outgroups] += 1

        if n_called_outgroups > 0:
            bases = np.array([variant.REF] + variant.ALT)
            ingroup_genotypes = np.unique(genotypes[ingroup_mask][:, [0, 1]])
            outgroup_genotypes = np.unique(genotypes[outgroup_mask][:, [0, 1]])
            ingroup_bases = bases[ingroup_genotypes[ingroup_genotypes != -1]]
            outgroup_bases = bases[outgroup_genotypes[outgroup_genotypes != -1]]

            # make sure bases are correct
            # assert len(outgroup_bases) == 1
            # assert outgroup_bases[0] in "".join(variant.gt_bases[outgroup_mask])
            # assert np.all([base in "".join(variant.gt_bases[ingroup_mask]) for base in ingroup_bases])

            # add 1 if outgroup is in ingroup, 0 otherwise
            divergence[int(outgroup_bases[0].upper() in ingroup_bases)] += 1
            base_counts[(tuple(str(s) for s in ingroup_bases), tuple(str(s) for s in outgroup_bases))] += 1

            pass

stats = pd.DataFrame([{
    'divergence_ratio': divergence[0] / sum(divergence.values()),  # ingroup/outgroup divergence
    'relative_coverage': coverage[2] / sum(coverage.values()),  # coverage of outgroup wrt ingroup
    'fraction_polymorphism': sum(coverage.values()) / n_sites,  # ingroup polymorphisms vs all polymorphisms
    'n_sites': n_sites,  # total number of sites
    'n_called_sites_ingroup': sum(coverage.values()),  # number of sites where ingroup is called
    'n_sites_outgroup_covered': coverage[2],  # number of sites where outgroups and ingroup are called
    'coverage': dict(coverage),  # coverage
    'divergence': dict(divergence),  # divergence
    'base_counts': str(dict(base_counts))  # base counts
}])

stats.to_csv(out)

pass
