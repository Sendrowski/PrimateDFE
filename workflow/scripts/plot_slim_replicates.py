"""
Simulated spectra per replicate and their departure from the \\texttt{fastDFE} expectation.

The upper row shows the \\texttt{SLiM} and \\texttt{fastDFE} spectra as counts, the lower row
their difference after normalizing each over its polymorphic bins, so the residual measures a
difference in shape rather than in the number of segregating sites.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

try:
    testing = False
    slim = list(snakemake.input.slim)
    fastdfe = list(snakemake.input.fastdfe)
    labels = list(snakemake.params.labels)
    out = snakemake.output[0]
except NameError:
    testing = True
    slim = [f"resources/slim_replicates/slim/n_replicate={r}/sfs.csv" for r in (1, 2, 3)]
    fastdfe = [f"resources/slim_replicates/sfs_slim/n_replicate={r}/sfs.csv" for r in (1, 2, 3)]
    labels = [f"replicate {r}" for r in (1, 2, 3)]
    out = "scratch/slim_replicates.png"


def counts(path: str) -> np.ndarray:
    """Polymorphic bins of the selected spectrum.

    The \\texttt{SLiM} spectra carry separate neutral and selected columns, the \\texttt{fastDFE}
    spectra a single column holding the selected spectrum.
    """
    df = pd.read_csv(path)
    values = df["selected"].to_numpy() if "selected" in df else df.iloc[:, 0].to_numpy()
    return values[1:-1]


def axes_center(fig, axes) -> float:
    """Horizontal centre of the plotting area, excluding the y-label gutter."""
    boxes = [ax.get_position() for ax in np.ravel(axes) if ax.get_visible()]
    return (min(b.x0 for b in boxes) + max(b.x1 for b in boxes)) / 2


fig, axes = plt.subplots(2, len(labels), figsize=(2.8 * len(labels), 4.0), dpi=300,
                         sharex=True, sharey="row")

for col, (label, f_slim, f_fd) in enumerate(zip(labels, slim, fastdfe)):
    a, b = counts(f_slim), counts(f_fd)
    x = np.arange(1, len(a) + 1)
    width = 0.4

    top = axes[0, col]
    top.bar(x - width / 2, a, width, label="SLiM", color="tab:blue")
    top.bar(x + width / 2, b, width, label="fastDFE", color="tab:orange")
    top.set_title(label, fontsize=14)

    bottom = axes[1, col]
    bottom.bar(x, a / a.sum() - b / b.sum(), color="tab:blue")
    bottom.axhline(0, color="black", lw=0.8)

    for ax in (top, bottom):
        ax.set_xticks(x[::2])
        ax.tick_params(labelsize=9)
    if col == 0:
        top.set_ylabel("SFS count", fontsize=12)
        bottom.set_ylabel("residual", fontsize=12)

fig.tight_layout(rect=(0, 0.10, 1, 1))

center = axes_center(fig, axes)
fig.supxlabel("derived allele count", fontsize=12, x=center, y=0.078)

handles, legend_labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc="lower center", ncol=2, frameon=True, fontsize=11,
           bbox_to_anchor=(center, -0.035))
fig.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
