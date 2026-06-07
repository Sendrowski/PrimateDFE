"""
Plot joint DFE statistics against Ne for multiple populations.
"""
from collections import defaultdict

import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

from parametrizer import Parametrizer
from populations import Populations
from visualization import DFEvsNePlotter

try:
    testing = False
    json_file = snakemake.input.json
    ne_file = snakemake.input.ne
    tree_file = snakemake.input.tree
    populations = snakemake.params.populations
    labels = snakemake.params.get("labels", populations)
    legend = snakemake.params.legend
    out = snakemake.output[0]
except NameError:
    testing = True
    #json_file = "results/dfe/Castellano/joint/original_ref/b/unfolded.8.gamma.del.noeps.json"
    json_file = "scratch/serialized.joint.json"
    #ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    ne_file = "resources/Castellano/Ne/all.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    populations = fd.JointInference.from_file(json_file).types
    labels = populations
    legend = True
    out = "scratch/dfe_Ne.png"

inf = fd.JointInference.from_file(json_file)
ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

inferences = inf.joint_inferences
dfes = {pop: inferences[pop].get_dfe() for pop in populations}

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
    reg_type='phylo'
)

fig = plotter.plot(
    file=out,
    show=testing,
    show_legend=legend,
    style='default'
)
