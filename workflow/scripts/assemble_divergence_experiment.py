"""
Assemble the catarrhini divergence experiment: per-population alpha with and
without divergence, alpha-vs-Ne regressions (OFF vs ON), and a consistency check
of DFE-based alpha against the MK-style ``get_alpha_divergence``.

Deliverables:
* a summary table CSV (one row per population);
* an alpha-vs-Ne figure (OFF vs ON side by side);
* a small stats JSON answering the two questions (does divergence narrow the
  alpha CI; does an alpha-Ne signal appear / strengthen).
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

try:
    off_files = list(snakemake.input.off)
    on_files = list(snakemake.input.on)
    ne_file = snakemake.input.ne
    pops = list(snakemake.params.populations)
    out_table = snakemake.output.table
    out_plot = snakemake.output.plot
    out_stats = snakemake.output.stats
except NameError:
    import glob
    off_files = sorted(glob.glob("results/dfe/divergence/*/*/summary.div_off.8.json"))
    on_files = sorted(glob.glob("results/dfe/divergence/*/*/summary.div_on.8.json"))
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    pops = None
    out_table = "scratch/divergence/experiment_table.csv"
    out_plot = "scratch/divergence/alpha_vs_Ne.png"
    out_stats = "scratch/divergence/experiment_stats.json"


def pop_of(path):
    # results/dfe/divergence/{ref}/{population}/summary.div_*.{n}.json
    return path.split("/")[-2]


def load(files):
    out = {}
    for f in files:
        with open(f) as fh:
            out[pop_of(f)] = json.load(fh)
    return out


off = load(off_files)
on = load(on_files)
ne = pd.read_csv(ne_file).set_index("label")["x"].to_dict()

rows = []
for p in sorted(set(off) | set(on)):
    o, n = off.get(p, {}), on.get(p, {})
    rows.append(dict(
        population=p,
        Ne=ne.get(p, np.nan),
        alpha_off=o.get("alpha", np.nan),
        alpha_off_ci_width=o.get("alpha_ci_width", np.nan),
        alpha_off_ci_lo=o.get("alpha_ci_lower", np.nan),
        alpha_off_ci_hi=o.get("alpha_ci_upper", np.nan),
        alpha_on=n.get("alpha", np.nan),
        alpha_on_ci_width=n.get("alpha_ci_width", np.nan),
        alpha_on_ci_lo=n.get("alpha_ci_lower", np.nan),
        alpha_on_ci_hi=n.get("alpha_ci_upper", np.nan),
        alpha_dfe_based=n.get("alpha_dfe_based", np.nan),
        alpha_divergence_based=n.get("alpha_divergence_based", np.nan),
        n_div_sel=n.get("n_div_sel", np.nan),
        n_div_neut=n.get("n_div_neut", np.nan),
    ))
df = pd.DataFrame(rows)
df.to_csv(out_table, index=False)


def regress(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return dict(n=int(len(x)), slope=None, intercept=None, pearson_r=None,
                   pearson_p=None, spearman_r=None, spearman_p=None)
    lr = stats.linregress(x, y)
    sr = stats.spearmanr(x, y)
    return dict(n=int(len(x)), slope=float(lr.slope), intercept=float(lr.intercept),
                pearson_r=float(lr.rvalue), pearson_p=float(lr.pvalue),
                spearman_r=float(sr.correlation), spearman_p=float(sr.pvalue))


reg_off = regress(df.Ne.to_numpy(float), df.alpha_off.to_numpy(float))
reg_on = regress(df.Ne.to_numpy(float), df.alpha_on.to_numpy(float))

mean_w_off = float(np.nanmean(df.alpha_off_ci_width))
mean_w_on = float(np.nanmean(df.alpha_on_ci_width))
med_w_off = float(np.nanmedian(df.alpha_off_ci_width))
med_w_on = float(np.nanmedian(df.alpha_on_ci_width))

result = dict(
    n_populations=int(len(df)),
    ci_width_off_mean=mean_w_off, ci_width_on_mean=mean_w_on,
    ci_width_off_median=med_w_off, ci_width_on_median=med_w_on,
    ci_width_ratio_median=(med_w_off / med_w_on) if med_w_on else None,
    alpha_vs_Ne_off=reg_off,
    alpha_vs_Ne_on=reg_on,
    dfe_vs_divergence_alpha=regress(df.alpha_dfe_based.to_numpy(float),
                                    df.alpha_divergence_based.to_numpy(float)),
)
with open(out_stats, "w") as fh:
    json.dump(result, fh, indent=2)

# ---- plot ----
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
for ax, mode, reg in [(axes[0], "off", reg_off), (axes[1], "on", reg_on)]:
    a = df[f"alpha_{mode}"].to_numpy(float)
    lo = df[f"alpha_{mode}_ci_lo"].to_numpy(float)
    hi = df[f"alpha_{mode}_ci_hi"].to_numpy(float)
    x = df.Ne.to_numpy(float)
    yerr = np.vstack([a - lo, hi - a])
    ax.errorbar(x, a, yerr=yerr, fmt="o", ms=4, capsize=2, alpha=0.7, lw=0.8)
    if reg["slope"] is not None:
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 50)
        ax.plot(xs, reg["intercept"] + reg["slope"] * xs, "r-", lw=1.5)
        ax.set_title(f"divergence {mode.upper()}  "
                     f"(r={reg['pearson_r']:.2f}, p={reg['pearson_p']:.3f}, "
                     f"slope={reg['slope']:.2e})")
    ax.set_xlabel("Ne")
    ax.axhline(0, color="grey", lw=0.6, ls="--")
axes[0].set_ylabel(r"$\alpha$")
fig.suptitle("Catarrhini: alpha vs Ne, polymorphism-only (OFF) vs +divergence (ON)")
fig.tight_layout()
fig.savefig(out_plot, dpi=150)
print(json.dumps(result, indent=2))
