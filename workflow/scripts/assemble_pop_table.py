import fastdfe as fd
import numpy as np
import pandas as pd

from parametrizer import Parametrizer
from populations import Populations

try:
    species_metadata_file = snakemake.input.species_metadata
    n_individuals_file = snakemake.input.n_individuals
    ne_file = snakemake.input.ne
    pnps_file = snakemake.input.ns
    spectra_file = snakemake.input.spectra
    ancestral_stats_file = snakemake.input.ancestral_stats
    outgroup_stats_file = snakemake.input.outgroup_stats
    dfe_types = snakemake.params.dfe_types
    dfe_files = dict(zip(dfe_types, snakemake.input.dfes))
    pops = snakemake.params.populations
    outgroups = snakemake.params.outgroups
    references = snakemake.params.references
    species = snakemake.params.species
    out = snakemake.output[0]
except NameError:
    species_metadata_file = 'resources/Kuderna/species.csv'
    n_individuals_file = 'results/stats/sample_names/counts/population_counts.csv'
    ne_file = 'results/stats/Ne/comp/original_ref/all/8.csv'
    pnps_file = 'results/stats/pNpS/comp/original_ref/catarrhini/8.csv'
    spectra_file = 'results/sfs/comp/original_ref/all/sfs.unfolded.8.csv'
    ancestral_stats_file = 'results/stats/comp/catarrhini/ancestral_stats.csv'
    outgroup_stats_file = 'results/stats/comp/catarrhini/outgroup_stats.8.csv'
    dfe_files = dict(
        dfe_discrete_8_deleterious="results/dfe/catarrhini/dfe.unfolded.8.gamma.del.noeps.csv"
    )
    pops = Populations.get_pops(8, 'catarrhini')
    references = {p: p for p in pops} # not accessible here
    species = {p: Populations.to_species(p) for p in pops}
    outgroups = {s: Populations.ingroup_outgroups[s] for s in species.values()}
    out = "scratch/pop_table.csv"

df = pd.DataFrame({'population': pops})
df['outgroup_species'] = [outgroups[species[p]] for p in pops]
df['n_outgroups'] = [len(outgroups[species[p]]) for p in pops]
df['reference'] = [references[p] for p in pops]
df['species'] = [species[p] for p in pops]

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
pnps = pd.read_csv(pnps_file)
pnps.rename(columns={'x': 'pnps', 'label': 'population'}, inplace=True)
#pnps = pnps.drop_duplicates(subset='population', keep='last')
df = df.merge(pnps[['population', 'pnps']], on='population', how='left')

# load ancestral stats and join
ancestral = pd.read_csv(ancestral_stats_file, index_col=1)
ancestral = ancestral.drop_duplicates(subset='population', keep='last')
ancestral['fraction_mismatches'] = ancestral['n_mismatches'] / ancestral['n_annotated']

df = df.merge(
    ancestral.add_prefix('ancestral_annotation.'),
    left_on='population',
    right_on='ancestral_annotation.population',
    how='left'
).drop(columns=['ancestral_annotation.population'])

# load outgroup stats and join
outgroup = pd.read_csv(outgroup_stats_file, index_col=1)
outgroup = outgroup.drop_duplicates(subset='population', keep='last')

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
df['sfs_neutral_8'] = [neutral[p].to_list() for p in pops]
df['sfs_selected_8'] = [selected[p].to_list() for p in pops]
df['theta'] = [neutral[p].theta for p in pops]
bins = np.array([-np.inf, -100, -10, -1, 0, 1, np.inf])

for key, file in dfe_files.items():
    df_dfes = pd.read_csv(file)
    df_dfes = df_dfes[df_dfes['population'].isin(pops)]
    dfes = {row.population: fd.DFE.from_json(row.json) for _, row in df_dfes.iterrows()}
    params = list(dfes[pops[0]].params.keys())
    df[key + '.params'] = [dfes[p].bootstraps[params].mean().to_dict() for p in pops]
    df[key + '.params_percentile_5'] = [
        dfes[p].bootstraps[params].quantile(0.05).to_dict() for p in pops]
    df[key + '.params_percentile_95'] = [
        dfes[p].bootstraps[params].quantile(0.95).to_dict() for p in pops]
    df[key + '.dfe_discretized'] = [dfes[p].discretize(bins)[0] for p in pops]
    df[key + '.dfe_discretized_percentile_5'] = [
        dfes[p].discretize(bins)[0] -
        dfes[p].discretize(bins)[1][0] for p in pops]
    df[key + '.dfe_discretized_percentile_95'] = [
        dfes[p].discretize(bins)[0] +
        dfes[p].discretize(bins)[1][1] for p in pops]

df.to_csv(out, index=False)
