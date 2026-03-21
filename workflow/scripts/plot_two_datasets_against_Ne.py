"""
Plot marginal DFE statistics against Ne for two datasets on shared axes.
"""
import fastdfe as fd
import pandas as pd

from populations import Populations
from visualization import DFEvsNePlotter

try:
    testing = False
    dfe_file_a = snakemake.input.dfe_a
    dfe_file_b = snakemake.input.dfe_b
    ne_file = snakemake.input.ne
    tree_file = snakemake.input.tree
    populations = snakemake.params.populations
    labels = snakemake.params.get("labels", populations)
    scaled = snakemake.params.get("scaled", False)
    legend = snakemake.params.legend
    out = snakemake.output[0]
except NameError:
    testing = True
    dfe_file_a = "results/dfe/catarrhini/dfe.unfolded.8.gamma.del.noeps.csv"
    dfe_file_b = "results/dfe/catarrhini/dfe.unfolded.filter_bgc.8.gamma.del.noeps.csv"
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    populations = Populations.get_pops(8, 'catarrhini')
    labels = [Populations.get_group_from_pop(p) for p in populations]
    dataset_labels = ["unfiltered", "GC-conservative"]
    scaled = True
    legend = True
    out = "scratch/dfe_Ne_both.png"

df_a = pd.read_csv(dfe_file_a)
dfes_a = {row.population: fd.DFE.from_json(row.json) for _, row in df_a.iterrows()}

df_b = pd.read_csv(dfe_file_b)
dfes_b = {row.population: fd.DFE.from_json(row.json) for _, row in df_b.iterrows()}

ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

if scaled:
    stat_list = [
        "range_S_-inf_-10",
        "range_S_-10_-1",
        "range_S_-1_0",
        #"range_S_0_inf",
        #"alpha"
    ]
else:
    stat_list = [
        "range_s_-inf_-1e-3",
        "range_s_-1e-3_-1e-5",
        "range_s_-1e-5_0",
        "range_s_0_inf",
    ]

plotter = DFEvsNePlotter(
    dfes=dfes_a,
    ne_dict=ne_dict,
    tree_file=tree_file,
    populations=populations,
    stat_list=stat_list,
    labels=labels,
    reg_type="linear",
)

plotter.plot_two_datasets_stacked(
    datasets={
        dataset_labels[0]: (dfes_a, ne_dict, labels, "C0"),
        dataset_labels[1]: (dfes_b, ne_dict, labels, "C1"),
    },
    file=out,
    show=testing,
)