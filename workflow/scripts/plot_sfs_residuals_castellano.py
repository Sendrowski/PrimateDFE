"""
Per-bin difference between the spectra of this study and those of Castellano et al.

Each spectrum is normalized over its polymorphic bins before differencing, so the two datasets
are compared on a common footing despite differing in the number of SNPs and in the reference
genome they were mapped to. Positive values indicate a larger proportion in this study.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

try:
    testing = False
    ours = list(snakemake.input.ours)
    theirs_neutral = list(snakemake.input.theirs_neutral)
    theirs_selected = list(snakemake.input.theirs_selected)
    labels = list(snakemake.params.labels)
    out = snakemake.output[0]
except NameError:
    testing = True
    import sys

    sys.path.insert(0, "workflow/scripts")
    from populations import Populations

    labels = Populations.names_castellano
    english = Populations.names_castellano_english
    refs = {p: ("Gorilla_gorilla" if p.startswith("Gorilla") else
                "Pan_troglodytes" if p.startswith("Pan_troglodytes") else p) for p in labels}
    ours = [f"results/sfs/{refs[p]}/{p}/sfs.unfolded.8.csv" for p in labels]
    theirs_neutral = [f"resources/Castellano/sfs/{e}_4fold_all_sfs.txt" for e in english]
    theirs_selected = [f"resources/Castellano/sfs/{e}_0fold_all_sfs.txt" for e in english]
    out = "scratch/castellano_residuals.png"

N_COLS = 3


def proportions(counts: np.ndarray) -> np.ndarray:
    """Polymorphic bins, scaled to sum to one."""
    poly = counts[1:-1]
    return poly / poly.sum()


def axes_center(fig, axes) -> float:
    """Horizontal centre of the plotting area, excluding the y-label gutter."""
    boxes = [ax.get_position() for ax in np.ravel(axes) if ax.get_visible()]
    return (min(b.x0 for b in boxes) + max(b.x1 for b in boxes)) / 2


n_rows = int(np.ceil(len(labels) / N_COLS))
fig, axes = plt.subplots(n_rows, N_COLS, figsize=(3.9 * N_COLS, 2.6 * n_rows), dpi=300,
                         sharex=True, sharey=True)
axes = np.atleast_2d(axes)

for ax, label, f_ours, f_neut, f_sel in zip(axes.flat, labels, ours, theirs_neutral, theirs_selected):
    mine = pd.read_csv(f_ours)
    diffs = {}
    for cls, f_theirs in [("neutral", f_neut), ("selected", f_sel)]:
        theirs = np.array([float(x) for x in open(f_theirs).read().strip().split(",")])
        diffs[cls] = proportions(mine[cls].to_numpy()) - proportions(theirs)

    x = np.arange(1, len(diffs["neutral"]) + 1)
    width = 0.4
    ax.bar(x - width / 2, diffs["neutral"], width, label="neutral", color="tab:blue")
    ax.bar(x + width / 2, diffs["selected"], width, label="selected", color="tab:orange")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_title(label.replace("_", " "), fontsize=15, style="italic")
    ax.set_xticks(x)
    ax.tick_params(labelsize=8)

for ax in axes.flat[len(labels):]:
    ax.set_visible(False)

fig.supylabel("difference in proportion", fontsize=12)

fig.tight_layout(rect=(0, 0.055, 1, 1))

center = axes_center(fig, axes)
fig.supxlabel("derived allele count", fontsize=12, x=center, y=0.055)

handles, legend_labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc="lower center", ncol=2, frameon=True, fontsize=11,
           bbox_to_anchor=(center, -0.02))
fig.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
