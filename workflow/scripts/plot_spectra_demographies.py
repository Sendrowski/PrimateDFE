"""
Simulated \\textit{selected} spectra under each demographic scenario.

Each panel contrasts the \\texttt{SLiM} simulation with the \\texttt{fastDFE} prediction obtained
using nuisance parameters. The constant-size spectrum is repeated in every panel as a reference,
rescaled to the total of that panel so that it marks the expected shape rather than a count.
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
    scenarios = ["constant", "expansion_4", "reduction_4", "bottleneck_20",
                 "substructure_0.0001", "dominance_0.2"]
    slim = [f"resources/slim_demographies/{d}/sfs.csv" for d in scenarios]
    fastdfe = [f"resources/slim_demographies/{d}/sfs.fastdfe.csv" for d in scenarios]
    labels = ["constant", "expansion", "reduction", "bottleneck", "substructure", "recessiveness"]
    out = "scratch/spectra_demographies.png"

N_COLS = 2


def counts(path: str) -> np.ndarray:
    """Polymorphic bins of the selected spectrum.

    The \\texttt{SLiM} spectra carry separate neutral and selected columns, the \\texttt{fastDFE}
    spectra a single column holding the selected spectrum.
    """
    df = pd.read_csv(path)
    values = df["selected"].to_numpy() if "selected" in df else df.iloc[:, 0].to_numpy()
    return values[1:-1]


constant = counts(slim[0])
reference_shape = constant / constant.sum()

n_rows = int(np.ceil(len(labels) / N_COLS))
fig, axes = plt.subplots(n_rows, N_COLS, figsize=(1.1 * 2048 / 300, 1.1 * 1317 / 300), dpi=300,
                         sharex=True)
axes = np.atleast_2d(axes)

for ax, label, f_slim, f_fd in zip(axes.flat, labels, slim, fastdfe):
    a, b = counts(f_slim), counts(f_fd)
    x = np.arange(1, len(a) + 1)
    width = 0.4
    ax.bar(x - width / 2, a, width, label="SLiM", color="tab:blue")
    ax.bar(x + width / 2, b, width, label="fastDFE", color="tab:orange")

    # constant-size shape at this panel's scale, padded so the end steps span a full bin
    scaled = reference_shape * a.sum()
    xs = np.concatenate([[x[0] - 1], x, [x[-1] + 1]])
    ys = np.concatenate([[scaled[0]], scaled, [scaled[-1]]])
    ax.step(xs, ys, where="mid", color="0.25", lw=1.3, alpha=0.5, label="constant size")
    ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)

    ax.set_title(label, fontsize=15)
    ax.set_xticks(x[::2])
    ax.tick_params(labelsize=9)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_fontsize(9)
    ax.legend(fontsize=8)

for ax in axes.flat[len(labels):]:
    ax.set_visible(False)

fig.tight_layout()
fig.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
