"""
Extract regression variables from DFE and Ne estimates.
"""
import fastdfe as fd
import numpy as np
import pandas as pd

from parametrizer import Parametrizer
from populations import Populations

try:
    testing = False
    dfe_file = snakemake.input.dfe
    ne_file = snakemake.input.ne
    populations = snakemake.params.populations
    out = snakemake.output[0]
except NameError:
    testing = True
    model = 'discrete'
    kind = 'full'
    dfe_file = f"results/dfe/catarrhini/dfe.unfolded.8.{model}.{kind}.noeps.csv"
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    kuderna_file = "resources/Kuderna/species.csv"
    populations = Populations.get_pops(8, 'catarrhini')
    out = f"scratch/reg_vars.{model}.{kind}.csv"


dfe_df = pd.read_csv(dfe_file)
dfes = {row.population: fd.DFE.from_json(row.json) for _, row in dfe_df.iterrows()}

ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

# add after loading ne_df
gen_df = (
    pd.read_csv(kuderna_file)[["SPECIES_BINOMIAL", "GENERATION_LENGTH"]]
    .rename(columns={"SPECIES_BINOMIAL": "species", "GENERATION_LENGTH": "generation_time"})
)

gen_map = dict(
    gen_df.assign(species=gen_df["species"])
    .set_index("species")["generation_time"]
)

rows = []
for pop in populations:
    species = Populations.to_species(pop)
    rows.append({
        'population': pop,
        'Ne': ne_dict[pop],
        'generation_time': gen_map.get(species, np.nan),
        'range_inf_-10': Parametrizer.get_S_range(dfes[pop], -np.inf, -10).mean(),
        'range_-10_-1': Parametrizer.get_S_range(dfes[pop], -10, -1).mean(),
        'range_-1_0': Parametrizer.get_S_range(dfes[pop], -1, 0).mean(),
        'range_0_inf': Parametrizer.get_S_range(dfes[pop], 0, np.inf).mean(),
    } | dfes[pop].params)

reg_df = pd.DataFrame(rows)
reg_df.to_csv(out, index=False)

pass
