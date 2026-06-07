"""
Additive-vs-recessive scatter showing that assuming additive fitness barely
changes the inferred DFE.

For each species and each parametrization (gamma, discrete) a point compares a
DFE statistic inferred under additive ($h=0.5$) inference (x-axis) against the
value under a partially recessive $h(S)$ (y-axis). Points cluster along the
$y=x$ line with only a small, systematic shift toward stronger selection;
parametrization is encoded by marker shape and taxonomic group by colour. One
panel per DFE statistic.
"""
import fastdfe as fd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from parametrizer import Parametrizer
from populations import Populations
from visualization import DFEvsNePlotter

try:
    testing = False
    stat_list = snakemake.params.stat_list
    add_files = {"gamma": snakemake.input.dfe_add_gamma,
                 "discrete": snakemake.input.dfe_add_disc}
    rec_files = {"gamma": snakemake.input.dfe_rec_gamma,
                 "discrete": snakemake.input.dfe_rec_disc}
    populations = snakemake.params.populations
    labels = snakemake.params.get("labels", populations)
    out = snakemake.output[0]
except NameError:
    testing = True
    base = "results/dfe/catarrhini/dfe.unfolded.8.{}.noeps.csv"
    add_files = {"gamma": base.format("gamma.del"),
                 "discrete": base.format("discrete.del")}
    rec_files = {"gamma": base.format("gamma.recessive"),
                 "discrete": base.format("discrete.recessive")}
    populations = Populations.get_pops(8, "catarrhini")
    labels = [Populations.get_group_from_pop(p) for p in populations]
    stat_list = ["range_S_-inf_-10", "range_S_-10_-1", "range_S_-1_0", "omega"]
    out = "scratch/dominance_scatter.png"

label_color = Populations.get_label_color_map(labels)
pop_label = dict(zip(populations, labels))
markers = {"gamma": "o", "discrete": "^"}


def load(path):
    df = pd.read_csv(path)
    return {row.population: fd.DFE.from_json(row.json) for _, row in df.iterrows()}


def stat_boot(dfe, stat):
    if stat == "omega":
        return Parametrizer.get_omega(dfe)
    if stat.startswith("range_S_"):
        lo, hi = DFEvsNePlotter._parse_range(stat)
        return Parametrizer.get_S_range(dfe, lo, hi)
    raise ValueError(stat)


add = {p: load(f) for p, f in add_files.items()}
rec = {p: load(f) for p, f in rec_files.items()}

ncol = 2
nrow = int(np.ceil(len(stat_list) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(7.0, 7.0), dpi=300)
axes = np.atleast_1d(axes).ravel()

for ax, stat in zip(axes, stat_list):
    xs, ys = [], []
    for param in ("gamma", "discrete"):
        for pop in populations:
            if pop not in add[param] or pop not in rec[param]:
                continue
            xm = np.median(stat_boot(add[param][pop], stat))
            ym = np.median(stat_boot(rec[param][pop], stat))
            ax.scatter(xm, ym, marker=markers[param], color=label_color[pop_label[pop]],
                       s=26, alpha=0.85, edgecolors="none", zorder=2)
            xs.append(xm)
            ys.append(ym)

    lo, hi = min(xs + ys), max(xs + ys)
    pad = 0.05 * (hi - lo)
    lim = (lo - pad, hi + pad)
    ax.plot(lim, lim, ls="--", color="0.5", lw=1.0, zorder=0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(DFEvsNePlotter.make_label(stat), fontsize=12)
    ax.tick_params(labelsize=8)

for ax in axes[len(stat_list):]:
    ax.set_visible(False)

fig.subplots_adjust(left=0.10, right=0.98, top=0.95, bottom=0.15,
                    hspace=0.30, wspace=0.25)
fig.supxlabel("additive ($h=0.5$) estimate", y=0.095, fontsize=12)
fig.supylabel("recessive estimate", x=0.015, fontsize=12)

order = sorted(set(labels), key=Populations.get_label_rank)
clade_handles = [Line2D([], [], marker="o", ls="none", color=label_color[l],
                        label=l.replace("_", " ")) for l in order]
param_handles = [Line2D([], [], marker=markers[p], ls="none", color="0.35",
                        label=p) for p in ("gamma", "discrete")]
fig.legend(handles=clade_handles + param_handles, loc="lower center",
           ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))

fig.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
else:
    plt.close(fig)
