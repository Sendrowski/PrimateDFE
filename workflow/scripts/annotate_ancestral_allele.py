"""
Annotate polarization of variants in a VCF file.
"""
import sys

import fastdfe as fd
import jsonpickle
import numpy as np
import pandas as pd

try:
    testing = False
    vcf_in = snakemake.input.vcf
    fasta = snakemake.input.fasta
    outgroups = snakemake.params.outgroups
    ingroups = snakemake.params.ingroups
    max_sites = snakemake.params.get("max_sites", np.inf)
    n_target_sites = snakemake.params.n_target_sites
    n_ingroups = snakemake.params.get("n_ingroups", 9)
    parallelize = False
    vcf_out = snakemake.output.vcf
    ann_out = snakemake.output.ann
    stats_out = snakemake.output.stats
    probs_out = snakemake.output.probs
except NameError:
    # testing
    testing = True
    vcf_in = "results/vcf/Papio_anubis/Papio_papio.degeneracy.existing_outgroups.vcf.gz"
    fasta = "resources/ref/Papio_anubis.fasta"
    outgroups = ['Papio_hamadryas']
    samples = pd.read_csv(f"resources/metadata/Papio_individuals.csv", sep="\t")
    ingroups = samples[samples["BAM_FOLDER"].str.contains('Papio_papio')].GVCF_ID.tolist()
    max_sites = 10000
    n_target_sites = 100000
    n_ingroups = 9
    parallelize = False
    vcf_out = 'scratch/Pongo_abelii.polarized.vcf.gz'
    ann_out = 'scratch/Pongo_abelii.polarized.ann.json'
    stats_out = 'scratch/Pongo_abelii.polarized.stats.txt'
    probs_out = 'scratch/Pongo_abelii.polarized.probs.png'

# initialize ancestral annotation
anc = fd.MaximumLikelihoodAncestralAnnotation(
    outgroups=outgroups,
    ingroups=ingroups,
    n_ingroups=n_ingroups,
    max_sites=max_sites,
    n_target_sites=n_target_sites,
    parallelize=parallelize,
    # prior=fd.AdaptivePolarizationPrior(),
    model=fd.K2SubstitutionModel(
        bounds={'K': (1e-8, 10), 'k': (0.1, 10)},
        fix_transition_transversion_ratio=True
    )
)

# initialize annotator
ann = fd.Annotator(
    vcf=vcf_in,
    fasta=fasta,
    output=vcf_out,
    annotations=[anc]
)

# run annotation
ann.annotate()

# save probabilities plot
anc.prior.plot(file=probs_out)

"""
no_minor_base = anc.configs.minor_base == -1
n_monomorphic = int(np.round(anc.configs[no_minor_base].multiplicity.sum()))
n_polymorphic = anc.n_sites - n_monomorphic
div = np.zeros(anc.n_outgroups)
for i in range(anc.n_outgroups):
    is_mismatch = anc.configs.outgroup_bases.map(lambda x: x[i]) != anc.configs.major_base
    div[i] = (is_mismatch.astype(int) * anc.configs.multiplicity)[~no_minor_base].sum() / n_polymorphic
"""

stats = pd.DataFrame({
    'n_sites': [anc.n_sites],
    'n_annotated': [anc.n_annotated],
    'n_mismatches': [len(anc.mismatches)],
    'outgroup_divergence': [str(anc.get_outgroup_divergence())],
    'params_mle': [str(anc.params_mle)],
})
print(stats.to_string(), file=sys.stderr)

# save stats
stats.to_csv(stats_out)

# save annotation
with open(ann_out, 'w') as f:
    anc._reader = None
    anc._handler = None
    anc.configs = anc.configs.to_json()
    f.write(jsonpickle.encode(anc, indent=4))

pass
