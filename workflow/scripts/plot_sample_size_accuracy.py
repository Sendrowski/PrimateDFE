"""
Inferred DFEs across SFS sample sizes under the constant population size scenario.

Black dashes mark the simulated (ground-truth) fraction in each interval. Ground truth is
converted to population-scaled coefficients with the effective population size estimated from
each simulation's neutral spectrum ($N_e = \\hat{\\theta} / 4\\mu$), which is below the census
size because of background selection.
"""
import json
import re

import fastdfe as fd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.container import BarContainer
from matplotlib.lines import Line2D

try:
    testing = False
    files = list(snakemake.input.inferences)
    sample_sizes = list(snakemake.params.sample_sizes)
    out = snakemake.output[0]
except NameError:
    testing = True
    base = ("resources/slim/n_replicate=1/n_chunks=100/g=1e4/L=1e7/mu=1e-8/r=1e-7/N=1e3/"
            "s_b={s_b}/b={b}/s_d={s_d}/p_b={p_b}/n={n}/constant/unfolded/inference.json")
    configs = [("1e-3", "0.3", "3e-1", "0.00"), ("1e-3", "0.1", "3e-2", "0.00"),
               ("1e-2", "0.1", "3e-1", "0.01"), ("1e-3", "0.3", "3e-2", "0.05")]
    sample_sizes = [8, 20, 100]
    files = [base.format(s_b=s_b, b=b, s_d=s_d, p_b=p_b, n=n)
             for s_b, b, s_d, p_b in configs for n in sample_sizes]
    out = "scratch/sample_size_accuracy.png"

INTERVALS = np.array([-np.inf, -100, -10, -1, 1, np.inf])


def simulated_params(path: str) -> dict:
    """Simulation parameters encoded in the result path."""
    return {k: float(re.search(rf"/{k}=([0-9.eE+-]+)", path).group(1))
            for k in ["s_d", "s_b", "b", "p_b", "N", "mu"]}


n_rows = len(files) // len(sample_sizes)
def axes_center(fig, axes) -> float:
    """Horizontal centre of the plotting area, excluding the y-label gutter."""
    boxes = [ax.get_position() for ax in np.ravel(axes) if ax.get_visible()]
    return (min(b.x0 for b in boxes) + max(b.x1 for b in boxes)) / 2


fig, axes = plt.subplots(n_rows, len(sample_sizes), figsize=(13.5, 9), dpi=150,
                         sharex=True, sharey="row")

for i, f in enumerate(files):
    row, col = divmod(i, len(sample_sizes))
    ax = axes[row, col]

    record = json.load(open(f))
    p = simulated_params(f)
    N_e = record["theta_neutral"] / (4 * p["mu"])

    truth = fd.DFE(dict(S_d=-4 * N_e * p["s_d"], b=p["b"], p_b=p["p_b"],
                        S_b=4 * N_e * p["s_b"], eps=0))
    values_true, _ = truth.discretize(INTERVALS, confidence_intervals=False)

    q = record["params_mle"]
    fd.DFE(q, bootstraps=pd.DataFrame(record["bootstraps"])).plot(
        intervals=INTERVALS,
        ax=ax,
        show=False,
        title=(f"$S_d$={q['S_d']:.0f}, $b$={q['b']:.2f}, "
               f"$p_b$={q['p_b']:.2f}, $S_b$={q['S_b']:.1f}"),
    )
    ax.set_title(ax.get_title(), fontsize=14)

    # ground truth as a thin line across each bar, so the bars and intervals are unchanged
    bars = [b for c in ax.containers if isinstance(c, BarContainer) for b in c]
    for bar, value in zip(bars, values_true):
        x = bar.get_x()
        ax.hlines(value, x + 0.06 * bar.get_width(), x + 0.94 * bar.get_width(),
                  color="tab:orange", lw=3.2, linestyles=(0, (1, 0.7)), zorder=6)

    if row < n_rows - 1:
        ax.set_xlabel("")
    else:
        ax.set_xlabel(ax.get_xlabel(), fontsize=13)
    if col > 0:
        ax.set_ylabel("")
    else:
        ax.set_ylabel(ax.get_ylabel(), fontsize=13)
    if row == 0:
        ax.text(0.5, 1.30, f"$n$ = {sample_sizes[col]}", transform=ax.transAxes,
                ha="center", fontsize=24)


fig.tight_layout()

center = axes_center(fig, axes)
fig.legend(handles=[Line2D([], [], color="tab:orange", lw=3.2, ls=(0, (1, 0.7)), label="SLiM (ground truth)")],
           loc="lower center", ncol=1, fontsize=13, frameon=True, bbox_to_anchor=(center, -0.055))
fig.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
