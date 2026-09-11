"""
Simulated \\textit{selected} spectra under the default settings, lower recombination rates and
strongly purifying selection. Each panel contrasts the \\texttt{SLiM} simulation with the
\\texttt{fastDFE} prediction and carries the simulated DFE parameters as a subtitle.
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

try:
    testing = False
    slim = list(snakemake.input.slim)
    fastdfe = list(snakemake.input.fastdfe)
    labels = list(snakemake.params.labels)
    subtitles = list(snakemake.params.subtitles)
    out = snakemake.output[0]
except NameError:
    testing = True
    scenarios = ["default", "r_1e-8", "r_1e-9", "s_d_100"]
    slim = [f"resources/slim_recombination/{d}/sfs.csv" for d in scenarios]
    fastdfe = [f"resources/slim_recombination/{d}/sfs.fastdfe.csv" for d in scenarios]
    labels = ["default ($r=10^{-7}$)", "normal ($r=10^{-8}$)", "reduced ($r=10^{-9}$)",
              "strongly purifying ($s_d=100$)"]
    subtitles = ["$s_d=0.3$, $b=0.3$, $p_b=0$, $s_b=10^{-3}$"] * 3 + ["$s_d=100$, $b=0.3$, $p_b=0$, $s_b=10^{-3}$"]
    out = "scratch/spectra_recombination.pdf"

N_COLS = 2
ROW_HEIGHT = 1.8
LEGEND_HEIGHT = 0.3


def counts(path: str) -> np.ndarray:
    """Polymorphic bins of the selected spectrum.

    The \\texttt{SLiM} spectra carry separate neutral and selected columns, the \\texttt{fastDFE}
    spectra a single column holding the selected spectrum.
    """
    df = pd.read_csv(path)
    values = df["selected"].to_numpy() if "selected" in df else df.iloc[:, 0].to_numpy()
    return values[1:-1]


n_rows = int(np.ceil(len(labels) / N_COLS))
height = ROW_HEIGHT * n_rows + LEGEND_HEIGHT
fig, axes = plt.subplots(n_rows, N_COLS, figsize=(7.5, height), dpi=300, sharex=True)
axes = np.atleast_2d(axes)

for ax, label, subtitle, f_slim, f_fd in zip(axes.flat, labels, subtitles, slim, fastdfe):
    a, b = counts(f_slim), counts(f_fd)
    x = np.arange(1, len(a) + 1)
    width = 0.4
    ax.bar(x - width / 2, a, width, label="SLiM", color="tab:blue")
    ax.bar(x + width / 2, b, width, label="fastDFE", color="tab:orange")
    ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)

    ax.set_title(label, fontsize=12.4, pad=24)
    ax.text(0.5, 1.05, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=10.4)
    ax.set_xticks(x[::2])
    ax.tick_params(labelsize=9)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_fontsize(9)

for ax in axes.flat[len(labels):]:
    ax.set_visible(False)

fig.tight_layout(h_pad=1.5, w_pad=1.5, rect=(0, LEGEND_HEIGHT / height, 1, 1))

handles, legend_labels = axes.flat[0].get_legend_handles_labels()
center = (axes[0, 0].get_position().x0 + axes[0, -1].get_position().x1) / 2
fig.legend(handles, legend_labels, loc="lower center", ncol=len(handles), frameon=True, fontsize=11,
           bbox_to_anchor=(center, 0))
fig.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
