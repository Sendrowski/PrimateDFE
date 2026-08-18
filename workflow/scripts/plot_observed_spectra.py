"""
Plot the observed neutral and selected SFS for every population as a grid of small multiples,
using fastDFE's own spectrum plotting so that bar style and colors match the other SFS figures.
Each spectrum is scaled to proportions over the polymorphic bins, and panels follow the
same order as the phylogeny figure.
"""

import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from populations import Populations

try:
    testing = False
    sfs_files = snakemake.input.sfs
    order_file = snakemake.input.order
    pops = snakemake.params.pops
    out = snakemake.output[0]
except NameError:
    testing = True


    def get_ref_from_pop(pop: str) -> str:
        s = "_".join(pop.split("_")[:2])
        df = pd.read_csv(f"resources/metadata/{s.split('_')[0]}_individuals.csv", sep="\t")
        return df[df["BAM_FOLDER"].str.contains(s)].iloc[0].REFERENCE_FOLDER.replace("_ssp", "")


    pops = Populations.get_pops(8, "catarrhini")
    sfs_files = [f"results/sfs/{get_ref_from_pop(p)}/{p}/sfs.unfolded.8.csv" for p in pops]
    order_file = "results/tables/phylogeny/catarrhini/species_order.8.csv"
    out = "scratch/observed_spectra.png"

N_COLS = 3

# -------------------- load --------------------

def normalize(spectra: 'fd.Spectra') -> 'fd.Spectra':
    """
    Scale each spectrum to proportions over the polymorphic bins, so that the neutral and selected
    shapes are directly comparable despite the excess of 0-fold over 4-fold target sites.
    """
    df = spectra.data.copy()
    df.iloc[0] = 0
    df.iloc[-1] = 0
    return fd.Spectra.from_dataframe(df / df.sum())


spectra = {pop: normalize(fd.Spectra.from_file(path)) for pop, path in zip(pops, sfs_files)}

# lay the panels out in the same order as the phylogeny figure, subspecies grouped by species
rank = pd.read_csv(order_file).set_index("species")["rank"]
order = sorted(spectra, key=lambda p: (rank.get("_".join(p.split("_")[:2]), len(rank)), p))

# -------------------- plot --------------------

n_rows = int(np.ceil(len(order) / N_COLS))
fig, axes = plt.subplots(n_rows, N_COLS, figsize=(3.6 * N_COLS, 1.08 * n_rows), dpi=300,
                         sharex=True, sharey=True)
axes = np.atleast_2d(axes)

handles, labels = [], []
for ax, pop in zip(axes.flat, order):
    spectra[pop].plot(ax=ax, show=False, title=pop.replace("_", " "),
                      kwargs_legend=dict(prop=dict(size=6)))
    ax.set_title(pop.replace("_", " "), fontsize=12, style="italic")
    ax.set_xlabel("")
    ax.tick_params(labelsize=9)
    ax.yaxis.get_offset_text().set_fontsize(6)
    if not handles:
        handles, labels = ax.get_legend_handles_labels()
    if ax.get_legend() is not None:
        ax.get_legend().remove()

for ax in axes.flat[len(order):]:
    ax.set_visible(False)

# the final row is partly empty, so restore x tick labels on the lowest visible panel per column
for col in range(N_COLS):
    visible = [axes[row, col] for row in range(n_rows) if axes[row, col].get_visible()]
    if visible:
        visible[-1].tick_params(labelbottom=True)

fig.supxlabel("derived allele count", fontsize=11, y=0.03)
fig.legend(handles, labels, loc="lower center", ncol=2, frameon=True, fontsize=11,
           bbox_to_anchor=(0.5, -0.015))
fig.tight_layout(rect=(0, 0.012, 1, 1))
fig.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
