import pandas as pd
from populations import Populations

df = pd.read_csv("Bjarke/apes_in_prado.csv")

in_prado = df[df['in_prado']]

# used populations have at least 4 individuals in Prado dataset
prado_pop_sizes = in_prado.groupby('species_full').size().reset_index(name='count')

new = pd.read_csv("results/stats/sample_names/all_samples.csv")

grouped = (
    new.groupby('sample')['population']
       .agg(list)
       .reset_index()
)

# find samples that map to multiple populations
# only alternative or nested populations contain the same samples
duplicates = grouped[grouped['population'].str.len() > 1]

new_pop_sizes = (
    new.groupby('population')
       .size()
       .reset_index(name='count')
)

# get all individuals from prado that are also in new samples
in_prado_and_new = in_prado[in_prado['ID'].isin(new['sample'])]

assert len(in_prado_and_new) == len(in_prado)

# double check individually for each population
for pop in Populations.names_castellano:
    assert set(in_prado[in_prado['species_full'].str.contains(pop)].ID).issubset(
        set(new[new['population'].str.contains(pop)]['sample'])
    )

pass