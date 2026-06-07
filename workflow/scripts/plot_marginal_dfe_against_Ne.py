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
    tree_file = snakemake.input.tree
    populations = snakemake.params.populations
    labels = snakemake.params.get("labels", populations)
    scaled = snakemake.params.get("scaled", False)
    legend = snakemake.params.legend
    stat_list = snakemake.params.stat_list
    style = snakemake.params.style
    out = snakemake.output[0]
except NameError:
    testing = True
    dfe_file = "results/dfe/catarrhini/dfe.unfolded.8.gamma.full.noeps.csv"
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    populations = Populations.get_pops(8, 'catarrhini')
    labels = [Populations.get_group_from_pop(p) for p in populations]
    scaled = True
    legend = True
    stat_list = ["range_S_-inf_-10", "range_S_-10_-1", "range_S_-1_0", "range_S_0_0.1", "range_S_0.1_inf"]
    style = "dfe"
    out = "scratch/dfe_Ne.png"

df = pd.read_csv(dfe_file)
dfes = {row.population: fd.DFE.from_json(row.json) for _, row in df.iterrows()}

ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

plotter = DFEvsNePlotter(
    dfes=dfes,
    ne_dict=ne_dict,
    tree_file=tree_file,
    populations=populations,
    stat_list=stat_list,
    labels=labels,
    reg_type="phylo",
)

plotter.plot(
    file=out,
    show=testing,
    show_legend=legend,
    legend_n_cols=6,
    style=style,
)
