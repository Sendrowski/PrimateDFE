"""
Plot alpha vs Ne for up to three datasets
using plot_alpha_three_datasets_stacked.
"""

import fastdfe as fd
import pandas as pd

from populations import Populations
from visualization import DFEvsNePlotter


try:
    testing = False
    dfe_file_a = snakemake.input.dfe_a
    dfe_file_b = snakemake.input.dfe_b
    dfe_file_c = snakemake.input.get("dfe_c", None)

    ne_file = snakemake.input.ne
    tree_file = snakemake.input.tree

    populations = snakemake.params.populations
    labels = snakemake.params.get("labels", populations)
    dataset_labels = snakemake.params.dataset_labels

    out = snakemake.output[0]

except NameError:
    testing = True

    dfe_file_a = "results/dfe/catarrhini/dfe.unfolded.8.gamma.full.noeps.csv"
    dfe_file_b = "results/dfe/catarrhini/dfe.folded.8.gamma.full.noeps.csv"
    dfe_file_c = "results/dfe/catarrhini/dfe.unfolded.8.discrete.full.noeps.csv"

    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"

    populations = Populations.get_pops(8, "catarrhini")
    labels = [Populations.get_group_from_pop(p) for p in populations]

    dataset_labels = ["gamma unfolded", "gamma folded", "discrete unfolded"]
    out = "scratch/alpha_Ne_stacked.png"


# ------------------------------------------------------------
# Load DFEs
# ------------------------------------------------------------

df_a = pd.read_csv(dfe_file_a)
dfes_a = {row.population: fd.DFE.from_json(row.json) for _, row in df_a.iterrows()}

df_b = pd.read_csv(dfe_file_b)
dfes_b = {row.population: fd.DFE.from_json(row.json) for _, row in df_b.iterrows()}

datasets = {
    dataset_labels[0]: (dfes_a, None, labels, "C0"),
    dataset_labels[1]: (dfes_b, None, labels, "C1"),
}

if dfe_file_c is not None:
    df_c = pd.read_csv(dfe_file_c)
    dfes_c = {row.population: fd.DFE.from_json(row.json) for _, row in df_c.iterrows()}
    datasets[dataset_labels[2]] = (dfes_c, None, labels, "C2")


# ------------------------------------------------------------
# Load Ne
# ------------------------------------------------------------

ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

# inject ne_dict into dataset tuples
datasets = {
    name: (dfes, ne_dict, labels, color)
    for name, (dfes, _, labels, color) in datasets.items()
}


# ============================================================
# Alpha stacked plot (NEW METHOD ONLY)
# ============================================================

plotter = DFEvsNePlotter(
    dfes=dfes_a,  # dummy init, swapped internally
    ne_dict=ne_dict,
    tree_file=tree_file,
    populations=populations,
    stat_list=["alpha"],
    labels=labels,
    reg_type="linear",
)

plotter.plot_alpha_three_datasets_stacked(
    datasets=datasets,
    file=out,
    show=testing,
)