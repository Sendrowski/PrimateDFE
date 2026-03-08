import fastdfe as fd
import numpy as np
import pandas as pd

from utils import Parametrizer

try:
    species_metadata_file = snakemake.input.species_metadata
    n_individuals_file = snakemake.input.n_individuals
    ne_file = snakemake.input.ne
    ns_file = snakemake.input.ns
    spectra_file = snakemake.input.spectra
    ancestral_stats_file = snakemake.input.ancestral_stats
    outgroup_stats_file = snakemake.input.outgroup_stats
    dfe_types = snakemake.params.dfe_types
    inference_files = dict(zip(dfe_types, snakemake.input.dfes))
    populations = snakemake.params.populations
    outgroups = snakemake.params.outgroups
    references = snakemake.params.references
    species = snakemake.params.species
    out = snakemake.output[0]
except NameError:
    species_metadata_file = 'resources/Kuderna/species.csv'
    n_individuals_file = 'results/stats/sample_names/counts/population_counts.csv'
    ne_file = 'results/stats/Ne/comp/original_ref/all/8.csv'
    ns_file = 'results/stats/comp/ns/ingroup/population/ns.csv'
    spectra_file = 'results/sfs/comp/original_ref/all/sfs.unfolded.8.csv'
    ancestral_stats_file = 'results/stats/comp/great_apes/ancestral_stats.csv'
    outgroup_stats_file = 'results/stats/comp/catarrhini/outgroup_stats.8.csv'
    inference_files = dict(
        dfe_discrete_8_deleterious_noeps="results/dfe/great_apes/joint/original_ref/S1+S2+S3/unfolded.8.discrete.del.noeps.json"
    )
    populations = [
        'Homo_sapiens',
        'Pan_troglodytes',
        'Pan_paniscus'
    ]
    outgroups = dict(
        Homo_sapiens=['Pan_troglodytes', 'Gorilla_gorilla', 'Pongo_abelii'],
        Pan_troglodytes=['Homo_sapiens', 'Gorilla_gorilla', 'Pongo_abelii'],
        Pan_paniscus=['Homo_sapiens', 'Gorilla_gorilla', 'Pongo_abelii'],
    )
    references = dict(
        Homo_sapiens='Homo_sapiens',
        Pan_troglodytes='Pan_troglodytes',
        Pan_paniscus='Pan_paniscus',
    )
    species = dict(
        Homo_sapiens='Homo_sapiens',
        Pan_troglodytes='Pan_troglodytes',
        Pan_paniscus='Pan_paniscus',
    )
    out = "out.csv"


def is_increasing(s: str):
    """
    Whether given array encoded as string is increasing.
    """
    return np.all(np.diff(np.fromstring(s.strip("[]"), sep=" ")) >= 0)


df = pd.DataFrame({'population': populations})
df['outgroup_species'] = [outgroups[species[p]] for p in populations]
df['n_outgroups'] = [len(outgroups[species[p]]) for p in populations]
df['reference'] = [references[p] for p in populations]
df['species'] = [species[p] for p in populations]

# load metadata and merge
metadata = pd.read_csv(species_metadata_file)
# prefix metadata columns
metadata = metadata.add_prefix('species_metadata.')
df = df.merge(metadata, left_on='species', right_on='species_metadata.SPECIES_BINOMIAL', how='left')

# load n_individuals and join
n_individuals = pd.read_csv(n_individuals_file)
n_individuals = n_individuals.drop_duplicates(subset='population', keep='last')
df = df.merge(n_individuals[['population', 'n_individuals']], on='population', how='left')

# load Ne and join
ne = pd.read_csv(ne_file)
ne.rename(columns={'x': 'ne_estimate', 'label': 'population'}, inplace=True)
ne = ne.drop_duplicates(subset='population', keep='last')
df = df.merge(ne[['population', 'ne_estimate']], on='population', how='left')

# load N/S and join
ns = pd.read_csv(ns_file)
ns = ns.drop_duplicates(subset='population', keep='last')
df = df.merge(ns[['population', 'ns']], on='population', how='left')

# load ancestral stats and join
ancestral = pd.read_csv(ancestral_stats_file, index_col=1)
ancestral = ancestral.drop_duplicates(subset='population', keep='last')
ancestral['fraction_mismatches'] = ancestral['n_mismatches'] / ancestral['n_annotated']

ancestral['increasing_outgroup_divergence'] = ancestral['outgroup_divergence'].apply(is_increasing)

df = df.merge(
    ancestral.add_prefix('ancestral_annotation.'),
    left_on='population',
    right_on='ancestral_annotation.population',
    how='left'
).drop(columns=['ancestral_annotation.population'])

# load outgroup stats and join
outgroup = pd.read_csv(outgroup_stats_file, index_col=1)
outgroup = outgroup.drop_duplicates(subset='population', keep='last')
outgroup['increasing_mismatch_fraction'] = outgroup['mismatch_fraction'].apply(is_increasing)

df = df.merge(
    outgroup.add_prefix('outgroup_stats.'),
    left_on='population',
    right_on='outgroup_stats.population',
    how='left'
).drop(columns=['outgroup_stats.population'])

# load SFS and join
spectra = fd.Spectra.from_file(spectra_file)
neutral = spectra['neutral.*'].merge_groups(level=1)
selected = spectra['selected.*'].merge_groups(level=1)
df['sfs_neutral_8'] = [neutral[p].to_list() for p in populations]
df['sfs_selected_8'] = [selected[p].to_list() for p in populations]
df['theta'] = [neutral[p].theta for p in populations]

for key, file in inference_files.items():
    inference = fd.JointInference.from_file(file)
    df[key + '.params'] = [inference.marginal_inferences[p].bootstraps.mean().to_dict() for p in populations]
    df[key + '.params_percentile_5'] = [
        inference.marginal_inferences[p].bootstraps.quantile(0.05).to_dict() for p in populations]
    df[key + '.params_percentile_95'] = [
        inference.marginal_inferences[p].bootstraps.quantile(0.95).to_dict() for p in populations]
    df[key + '.dfe_discretized'] = [inference.marginal_inferences[p].get_discretized()[0] for p in populations]
    df[key + '.dfe_discretized_percentile_5'] = [
        inference.marginal_inferences[p].get_discretized()[0] -
        inference.marginal_inferences[p].get_discretized()[1][0] for p in populations]
    df[key + '.dfe_discretized_percentile_95'] = [
        inference.marginal_inferences[p].get_discretized()[0] +
        inference.marginal_inferences[p].get_discretized()[1][1] for p in populations]
    df[key + '.log_likelihood'] = [inference.marginal_inferences[p].likelihood for p in populations]
    df[key + '.L2_residual'] = [inference.marginal_inferences[p].L2_residual for p in populations]
    df[key + '.S_d'] = [np.mean(Parametrizer.get_S_d(inference.marginal_inferences[p].get_dfe())) for p in populations]

df.to_csv(out, index=False)
