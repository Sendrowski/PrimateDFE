"""
Plot marginal DFE statistics against Ne for multiple populations.
"""

import fastdfe as fd
import pandas as pd

from populations import Populations
from visualization import DFEvsNePlotter

try:
    testing = False
    dfe_file = snakemake.input.dfe
    ne_file = snakemake.input.ne
    populations = snakemake.params.populations
    labels = snakemake.params.get("labels", populations)
    legend = snakemake.params.legend
    out = snakemake.output[0]
except NameError:
    testing = True
    dfe_file = "results/dfe/catarrhini/dfe.unfolded.8.discrete.full.noeps.csv"
    # dfe_file = "results/dfe/catarrhini/dfe.unfolded.8.gamma.fixed_S_b=1.noeps.csv"
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    populations = pd.read_csv(dfe_file)["population"].tolist()
    labels = [Populations.get_group_from_pop(p) for p in populations]
    legend = True
    out = "scratch/dfe_Ne.png"

# -------------------- load --------------------

df = pd.read_csv(dfe_file)
dfes = {row.population: fd.DFE.from_json(row.json) for _, row in df.iterrows()}

ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

stat_list = ["S_b", "p_b", "range_S_-inf_-10", "range_S_0_1", 'range_S_-1_0']
plotter = DFEvsNePlotter(
    dfes=dfes,
    ne_dict=ne_dict,
    populations=populations,
    stat_list=stat_list,
    labels=labels,
)

fig, axes = plotter.plot(
    file=out,
    figsize=(8, 7),
    show=testing,
    show_legend=legend,
)
