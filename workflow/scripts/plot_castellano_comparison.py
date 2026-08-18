"""
DFE statistics against Ne for the estimates of Castellano et al. and for the present dataset.

One column per dataset and one row per statistic. Each dataset is drawn with the Ne estimates
of its own study.
"""
import fastdfe as fd
import pandas as pd
from fastdfe import DFE, GammaExpParametrization

from visualization import DFEvsNePlotter

try:
    testing = False
    dfe_file = snakemake.input.dfe
    ne_castellano = snakemake.input.ne_castellano
    json_file = snakemake.input.json
    ne_new = snakemake.input.ne_new
    tree_file = snakemake.input.tree
    populations = snakemake.params.populations
    out = snakemake.output[0]
except NameError:
    testing = True
    dfe_file = "resources/Castellano/dfe/3s.csv"
    ne_castellano = "resources/Castellano/Ne/all.csv"
    json_file = "results/dfe/Castellano/joint/original_ref/b/unfolded.8.gamma.del.noeps.json"
    ne_new = "results/stats/Ne/comp/original_ref/Castellano/8.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    populations = pd.read_csv(ne_castellano)["label"].tolist()
    out = "scratch/castellano_comparison.png"


def make_dfe_from_ci_row(row) -> DFE:
    """Build a DFE whose bootstraps are the reported lower, central and upper estimates,
    so the published confidence intervals propagate to every derived statistic.
    """
    boot = pd.DataFrame({
        "S_d": [row.S_d_low, row.S_d, row.S_d_high],
        "b": [row.b_low, row.b, row.b_high],
        "S_b": [0.0, 0.0, 0.0],
        "p_b": [0.0, 0.0, 0.0],
    })

    params = dict(S_d=row.S_d, b=row.b, S_b=0.0, p_b=0.0)

    return DFE(params=params, model=GammaExpParametrization(), bootstraps=boot)


dfe_df = pd.read_csv(dfe_file).set_index("species")
dfes_castellano = {pop: make_dfe_from_ci_row(dfe_df.loc[pop]) for pop in populations}

inferences = fd.JointInference.from_file(json_file).joint_inferences
dfes_new = {pop: inferences[pop].get_dfe() for pop in populations}

ne_dicts = []
for f in (ne_castellano, ne_new):
    df = pd.read_csv(f)
    ne_dicts.append(dict(zip(df["label"], df["x"])))

stat_list = [
    "S_d",
    "s_d",
    "range_S_-inf_-10",
    "range_S_-10_-1",
    "range_S_-1_0",
]

plotter = DFEvsNePlotter(
    dfes=dfes_castellano,
    ne_dict=ne_dicts[0],
    populations=populations,
    stat_list=stat_list,
    labels=populations,
    label_color={p: f"C{idx % 10}" for idx, p in enumerate(populations)},
    tree_file=tree_file,
    reg_type="phylo"
)

plotter.label_order = populations

plotter.plot_two_datasets_columns(
    datasets={
        "Castellano et al.": (dfes_castellano, ne_dicts[0], populations),
        "This study": (dfes_new, ne_dicts[1], populations),
    },
    legend_n_cols=3,
    file=out,
    show=testing,
)
