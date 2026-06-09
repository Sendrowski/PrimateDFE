"""
Plot marginal DFE statistics against Ne for two (or four) datasets on shared axes.

The optional dfe_c/dfe_d inputs enable a four-way overlay (used for the dominance
x parametrization figure): parametrization is encoded by colour+marker
(gamma: black/circle, discrete: grey/triangle) and dominance by line style
(additive: solid, recessive: dashed).
"""
import fastdfe as fd
import pandas as pd

from populations import Populations
from visualization import DFEvsNePlotter

try:
    testing = False
    stat_list = snakemake.params.stat_list
    dfe_files = [snakemake.input.dfe_a, snakemake.input.dfe_b]
    for extra in ("dfe_c", "dfe_d"):
        f = getattr(snakemake.input, extra, None)
        if f:
            dfe_files.append(f)
    ne_file = snakemake.input.ne
    tree_file = snakemake.input.tree
    populations = snakemake.params.populations
    labels = snakemake.params.get("labels", populations)
    dataset_labels = snakemake.params.dataset_labels
    scaled = snakemake.params.get("scaled", False)
    legend = snakemake.params.legend
    color_by = snakemake.params.get("color_by", "dataset")
    out = snakemake.output[0]
except NameError:
    testing = True
    base = "results/dfe/catarrhini/dfe.unfolded.8.{}.noeps.csv"
    dfe_files = [base.format("gamma.del"), base.format("gamma.recessive"),
                 base.format("discrete.del"), base.format("discrete.recessive")]
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    populations = Populations.get_pops(8, 'catarrhini')
    labels = [Populations.get_group_from_pop(p) for p in populations]
    dataset_labels = ["gamma additive", "gamma recessive", "discrete additive", "discrete recessive"]
    scaled = True
    legend = True
    color_by = "taxon"
    stat_list = ["range_S_-inf_-10", "range_S_-10_-1", "range_S_-1_0", "omega"]
    out = "scratch/dfe_Ne_both.png"

dfes_list = []
for f in dfe_files:
    df = pd.read_csv(f)
    dfes_list.append({row.population: fd.DFE.from_json(row.json) for _, row in df.iterrows()})

ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

plotter = DFEvsNePlotter(
    dfes=dfes_list[0],
    ne_dict=ne_dict,
    tree_file=tree_file,
    populations=populations,
    stat_list=stat_list,
    labels=labels,
    reg_type="phylo",
)

n_sets = len(dfes_list)
if color_by == "taxon":
    if n_sets == 4:
        # parametrization -> colour + marker; dominance -> line style
        line_colors = ["black", "black", "0.25", "0.25"]
        point_markers = ["o", "o", "^", "^"]
        line_styles = ["-", "--", "-", "--"]
    else:
        # colour by taxonomic group, distinguish datasets by marker + line style
        line_colors = ["black", "0.25"]
        point_markers = ["o", "^"]
        line_styles = ["-", "--"]
else:
    line_colors = ["C0", "C1"]
    point_markers = ["o", "o"]
    line_styles = ["--", "--"]

# draw the GammaExpParametrization regression lines in a softened red so they
# stand out against the (darker grey) discrete parametrization
line_colors = ["#ff6b6b" if lab.lower().startswith("gamma") else c
               for lab, c in zip(dataset_labels, line_colors)]

datasets, markers, linestyles = {}, {}, {}
for i, (lab, dfes) in enumerate(zip(dataset_labels, dfes_list)):
    datasets[lab] = (dfes, ne_dict, labels, line_colors[i])
    markers[lab] = point_markers[i]
    linestyles[lab] = line_styles[i]

plotter.plot_two_datasets_stacked(
    datasets=datasets,
    color_by=color_by,
    markers=markers,
    linestyles=linestyles,
    file=out,
    show=testing,
)
