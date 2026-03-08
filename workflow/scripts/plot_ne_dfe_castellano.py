"""
Plot DFE statistic against Ne for multiple populations.
"""

import pandas as pd
from fastdfe import DFE, GammaExpParametrization
from visualization import DFEvsNePlotter

try:
    testing = False
    dfe_file = snakemake.input.dfe
    ne_file = snakemake.input.ne
    populations = snakemake.params.populations
    tree_file = snakemake.input.tree
    labels = snakemake.params.get("labels", populations)
    out = snakemake.output[0]
except NameError:
    testing = True
    dfe_file = "resources/Castellano/dfe/3s_gc.csv"
    ne_file = "resources/Castellano/Ne/all.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    populations = pd.read_csv(ne_file)["label"].tolist()
    labels = populations
    out = "scratch/dfe_Ne.png"


def make_dfe_from_ci_row(row):
    """
    Build a fastDFE.DFE whose bootstrap DFEs are exactly:
    (low, mid, high) from Castellano confidence intervals.
    This propagates CIs deterministically to all derived stats.
    """
    Sd_mid  = row.S_d
    Sd_low  = row.S_d_low
    Sd_high = row.S_d_high

    b_mid  = row.b
    b_low  = row.b_low
    b_high = row.b_high

    boot = pd.DataFrame({
        "S_d": [Sd_low,  Sd_mid,  Sd_high],
        "b":   [b_low,   b_mid,   b_high],
        "S_b": [0.0,     0.0,     0.0],
        "p_b": [0.0,     0.0,     0.0],
    })

    params = dict(
        S_d = Sd_mid,
        b   = b_mid,
        S_b = 0.0,
        p_b = 0.0,
    )

    return DFE(params=params, model=GammaExpParametrization(), bootstraps=boot)


dfe_df = pd.read_csv(dfe_file).set_index("species")
ne_df  = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

dfes = {
    pop: make_dfe_from_ci_row(dfe_df.loc[pop])
    for pop in populations
}

stat_list = [
    "S_d",
    "s_d",
    "range_S_-inf_-10",
    "range_S_-10_-1",
    "range_S_-1_0",
]

plotter = DFEvsNePlotter(
    dfes=dfes,
    ne_dict=ne_dict,
    populations=populations,
    stat_list=stat_list,
    labels=labels,
    label_color={p: f"C{idx%10}" for idx, p in enumerate(populations)},
    tree_file=tree_file,
    reg_type="phylo"
)

plotter.plot(file=out, show=testing)