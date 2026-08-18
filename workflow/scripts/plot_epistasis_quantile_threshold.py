"""
Tail-robust quantile test for the sign of epistasis on the deleterious DFE,
across the four DFE-model variants ({gamma,discrete} x {del,full}).

Instead of the mean deleterious effect S_d (a tail-dominated, parametrization-
dependent quantity whose Ne-slope flips sign across variants), we regress the
deleterious-mass quantile |S_q| on N_e. Under no epistasis a fixed quantile
rescales as S_q = 4 N_e s_q, so the slope of log|S_q| on log N_e is 1 (dashed
reference); deviation from 1 measures epistasis. The reported p-value tests
H0: slope = 1, obtained by regressing the unscaled quantile log|s_q|
(= log|S_q| - log(4 N_e)) on log N_e.

The quantile is robust to the parametrization because it fixes a cumulative-
probability level (a self-normalising, observable position for every species),
unlike an iso-threshold whose scaled evaluation point drifts into the poorly-
constrained strong tail at large N_e. Each DFE-model variant occupies its own
panel, with colour encoding the quantile level f.
"""
import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator

from parametrizer import Parametrizer
from populations import Populations
from regression import Regression

# variant key -> (legend label, colour). Muted Seaborn palette matching the
# other manuscript figures; gamma = warm, discrete = cool.
VARIANTS = {
    "gamma_del": ("gamma (del)", "#C44E52"),    # muted red
    "gamma_full": ("gamma (full)", "#DD8452"),  # muted orange
    "discrete_del": ("discrete (del)", "#4C72B0"),    # muted blue
    "discrete_full": ("discrete (full)", "#8172B2"),  # muted purple
}

def log_label_pow(x, pos):
    """Format log-scaled x-axis labels as powers of 10 (matching Figure 2)."""
    if x <= 0:
        return ""
    k = int(np.floor(np.log10(x)))
    n = x / 10 ** k
    if np.isclose(n, 1.0):
        return rf"$10^{{{k}}}$"
    return rf"${int(round(n))}\times10^{{{k}}}$"


try:
    testing = False
    dfe_files = {
        "gamma_del": snakemake.input.gamma_del,
        "gamma_full": snakemake.input.gamma_full,
        "discrete_del": snakemake.input.discrete_del,
        "discrete_full": snakemake.input.discrete_full,
    }
    ne_file = snakemake.input.ne
    tree_file = snakemake.input.tree
    populations = snakemake.params.populations
    quantiles = snakemake.params.get("quantiles", [0.4, 0.5, 0.6, 0.7])
    out = snakemake.output.plot
    table_out = snakemake.output.table
except NameError:
    testing = True
    base = "results/dfe/catarrhini/dfe.unfolded.8.{}.noeps.csv"
    dfe_files = {
        "gamma_del": base.format("gamma.del"),
        "gamma_full": base.format("gamma.full"),
        "discrete_del": base.format("discrete.del"),
        "discrete_full": base.format("discrete.full"),
    }
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    populations = Populations.get_pops(8, "catarrhini")
    quantiles = [0.4, 0.5, 0.6]
    out = "scratch/epistasis_quantile.png"
    table_out = "scratch/epistasis_quantile.tex"

ne_df = pd.read_csv(ne_file)
ne = dict(zip(ne_df["label"], ne_df["x"]))

dfes = {}
for key, path in dfe_files.items():
    df = pd.read_csv(path)
    dfes[key] = {row.population: fd.DFE.from_json(row.json) for _, row in df.iterrows()}


def summarize(boot: dict) -> dict:
    """Median and 2.5/97.5 percentile half-widths per population."""
    out = {}
    for pop, vals in boot.items():
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue
        m = np.median(vals)
        out[pop] = (m, m - np.percentile(vals, 2.5), np.percentile(vals, 97.5) - m)
    return out


def fit_line(stat: dict, transform):
    """PGLS fit of transform(median, Ne) on log(Ne). Returns (x_ln_ne, intercept, slope, p, r)."""
    pops = list(stat)
    x = np.array([np.log(ne[p]) for p in pops])
    y = np.array([transform(stat[p][0], ne[p]) for p in pops])
    ok = np.isfinite(x) & np.isfinite(y)
    pops = list(np.array(pops)[ok])
    x, y = x[ok], y[ok]
    intercept, slope, p, r = Regression.regress(
        "phylo", x=x, y=y, pops=np.array(pops), tree_file=tree_file
    )
    return x, intercept, slope, p, r


# pre-compute all quantile bootstrap distributions per variant (one CDF pass each)
stats = {}  # (key, f) -> summarized dict
for key in VARIANTS:
    pops = [p for p in populations if p in dfes[key] and p in ne]
    for p in pops:
        qd = Parametrizer.get_S_quantiles(dfes[key][p], quantiles)
        for f in quantiles:
            stats.setdefault((key, f), {})[p] = np.abs(qd[f])
stats = {k: summarize(v) for k, v in stats.items()}

# line style encodes the quantile level f; colour encodes the model variant
fstyles = ["-", "--", ":", "-."]
fstyle = {f: fstyles[i % len(fstyles)] for i, f in enumerate(quantiles)}

# fit each (variant, f) once: scaled slope b1 for the drawn line, plus the
# H0: slope=1 test (slope, r, p) via the unscaled quantile
# log|s_q| = log|S_q| - log(4 Ne).
fits = {}
for f in quantiles:
    for key in VARIANTS:
        stat = stats[(key, f)]
        if not stat:
            continue
        x, b0, b1, _, _ = fit_line(stat, lambda m, n: np.log(m))
        _, _, _, p_null, r_dev = fit_line(stat, lambda m, n: np.log(m) - np.log(4 * n))
        fits[(key, f)] = dict(xmin=x.min(), xmax=x.max(), b0=b0, b1=b1, p=p_null, r=r_dev)

# one panel per DFE-model variant, so colour is free to encode the quantile level f
QCOLORS = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]
qcolor = {f: QCOLORS[i % len(QCOLORS)] for i, f in enumerate(quantiles)}

fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.3), dpi=300, sharex=True, sharey=True)

# log scales are set before any data is drawn, so the shared autoscale sees
# every panel on the final scale
for ax in axes.flat:
    ax.set_xscale("log")
    ax.set_yscale("log")

xref = np.array([np.log(min(ne[p] for p in ne)), np.log(max(ne[p] for p in ne))])

for ax, (key, (label, _)) in zip(axes.flat, VARIANTS.items()):
    for f in quantiles:
        if (key, f) not in fits:
            continue
        stat, fit = stats[(key, f)], fits[(key, f)]
        for p in stat:
            m, lo, hi = stat[p]
            ax.errorbar(ne[p], m, yerr=[[max(lo, 0)], [max(hi, 0)]], fmt="o",
                        ms=3, color=qcolor[f], alpha=0.25, capsize=0, lw=0,
                        markeredgewidth=0)
        xx = np.linspace(fit["xmin"], fit["xmax"], 200)
        ax.plot(np.exp(xx), np.exp(fit["b0"] + fit["b1"] * xx), "-",
                color=qcolor[f], lw=2.0, solid_capstyle="round")
    ax.set_title(label, fontsize=14)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax.xaxis.set_major_formatter(FuncFormatter(log_label_pow))
    ax.grid(True, which="major", ls=":", lw=0.5, color="0.8", alpha=0.7)
    ax.set_axisbelow(True)

    # r and p test H0: slope = 1, i.e. no epistasis
    panel_handles = [
        Line2D([], [], color=qcolor[f], lw=1.8,
               label=fr"$f={f:g}$: r={fits[(key, f)]['r']:.2f}, p={fits[(key, f)]['p']:.2g}")
        for f in quantiles if (key, f) in fits
    ]
    ax.legend(handles=panel_handles, fontsize=8, loc="lower right",
              frameon=True, framealpha=0.85, borderpad=0.4)

# bound both axes by the plotted values, error bars included
vals = [(ne[pop], m, max(m - lo, 0), m + hi)
        for key_f in fits for pop, (m, lo, hi) in stats[key_f].items()]
xs = [v[0] for v in vals]
ys = [y for v in vals for y in (v[1], v[2], v[3]) if y > 0]
axes[0, 0].set_xlim(min(xs) * 0.85, max(xs) * 1.18)
axes[0, 0].set_ylim(min(ys) * 0.5, max(ys) * 2.0)

# slope-1 (no-epistasis) reference in every panel, anchored at the geometric
# centre of the shared y-range so it sits amid the fits
ymin, ymax = axes[0, 0].get_ylim()
anchor = np.log(np.sqrt(ymin * ymax)) - np.mean(xref)
for ax in axes.flat:
    ax.plot(np.exp(xref), np.exp(anchor + 1.0 * xref), color="0.25",
            lw=3.0, ls=(0, (7, 4)), solid_capstyle="round", zorder=1.5, alpha=0.9)

fig.supxlabel(r"$N_e$", fontsize=13)
fig.supylabel(r"$|S_q|$ (deleterious quantile)", fontsize=12)

handles = [Line2D([], [], color="0.25", lw=3.0, ls=(0, (7, 4)),
                  label="no epistasis (slope $=1$)")]

fig.tight_layout()
fig.legend(handles=handles, fontsize=9, ncol=len(handles), loc="lower center",
           frameon=True, bbox_to_anchor=(0.5, -0.045))

fig.savefig(out, bbox_inches="tight")


# ---- companion LaTeX table: PGLS slope with r and p for H0: slope = 1 ----
# Rendered as a plain `tabular` fragment for \input into the manuscript. The
# four DFE-model variants are the columns; rows are grouped by quantile level
# f, with slope, r and p as sub-rows within each group.
def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "---"
    if p >= 0.01:
        return f"{p:.2f}"
    exp = int(f"{p:e}".split('e')[1])
    mant = p / 10 ** exp
    return f"${mant:.1f}\\!\\times\\!10^{{{exp}}}$"


# variant labels used as column headers (spelled out)
tex_labels = {
    "gamma_del": r"gamma del",
    "gamma_full": r"gamma full",
    "discrete_del": r"discrete del",
    "discrete_full": r"discrete full",
}

header = " & ".join([r"$f$", ""] + [tex_labels[k] for k in VARIANTS])
lines = [
    r"\begin{tabular}{ll" + "r" * len(VARIANTS) + r"}",
    r"    \toprule",
    "    " + header + r" \\",
    r"    \midrule",
]
stat_rows = [
    ("slope", lambda fit: fr"${fit['b1']:.2f}$"),
    (r"$r$", lambda fit: fr"${fit['r']:.2f}$"),
    (r"$p$", lambda fit: fmt_p(fit['p'])),
]
for i, f in enumerate(quantiles):
    if i > 0:
        lines.append(r"    \midrule")
    for j, (stat_name, render) in enumerate(stat_rows):
        flabel = fr"${f:g}$" if j == 0 else ""
        cells = [render(fits[(key, f)]) if (key, f) in fits else "---"
                 for key in VARIANTS]
        lines.append("    " + " & ".join([flabel, stat_name] + cells) + r" \\")
lines += [r"    \bottomrule", r"\end{tabular}", ""]

with open(table_out, "w") as fh:
    fh.write("\n".join(lines))

if testing:
    plt.show()
