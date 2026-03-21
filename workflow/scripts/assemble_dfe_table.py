import fastdfe as fd
import numpy as np
import pandas as pd

from parametrizer import Parametrizer
from populations import Populations

try:
    dfe_file = snakemake.input[0]
    pops = snakemake.params.populations
    out = snakemake.output[0]
except NameError:
    dfe_file = "results/dfe/catarrhini/dfe.unfolded.8.gamma.del.noeps.csv"
    pops = Populations.get_pops(8, 'catarrhini')
    out = "scratch/dfe_discrete_8_deleterious.csv"

df = pd.DataFrame({'population': pops})

bins = np.array([-np.inf, -100, -10, -1, 0, 1, np.inf])

df_dfes = pd.read_csv(dfe_file)
df_dfes = df_dfes[df_dfes['population'].isin(pops)]
dfes = {row.population: fd.DFE.from_json(row.json) for _, row in df_dfes.iterrows()}
params = list(dfes[pops[0]].params.keys())
df['params'] = [dfes[p].bootstraps[params].mean().to_dict() for p in pops]
df['params_ci_5'] = [
    dfes[p].bootstraps[params].quantile(0.05).to_dict() for p in pops]
df['params_ci_95'] = [
    dfes[p].bootstraps[params].quantile(0.95).to_dict() for p in pops]
df['discretized'] = [dfes[p].discretize(bins)[0] for p in pops]
df['discretized_ci_5'] = [
    dfes[p].discretize(bins)[0] -
    dfes[p].discretize(bins)[1][0] for p in pops]
df['discretized_ci_95'] = [
    dfes[p].discretize(bins)[0] +
    dfes[p].discretize(bins)[1][1] for p in pops]

df.to_csv(out, index=False)
