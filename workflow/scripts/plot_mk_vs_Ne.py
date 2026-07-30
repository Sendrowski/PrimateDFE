"""
Species-wise McDonald-Kreitman alpha (= 1 - NI) against Ne, as a scatter.

Same pooled-to-species MK counts as plot_mk_test.py; Ne is averaged over the
subspecies of each species. OLS fit + Pearson/Spearman reported.
"""
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import fastdfe as fd
from populations import Populations

warnings.filterwarnings("ignore")

try:
    sfs_file = snakemake.input.sfs
    ne_file = snakemake.input.ne
    out = snakemake.output[0]
except NameError:
    sfs_file = "results/divergence/sfs/comp/catarrhini/sfs.8.csv"
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    out = "results/graphs/divergence/catarrhini/mk_vs_Ne.8.png"


def to_species(pop):
    return pop.rsplit("_", 1)[0] if pop.count("_") >= 2 else pop


spectra = fd.Spectra.from_file(sfs_file)
pops = sorted({t.split(".", 1)[1] for t in spectra.types if t.startswith("neutral.")})
ne_pop = dict(pd.read_csv(ne_file)[["label", "x"]].values)

rec = {}
for p in pops:
    neut, sel = spectra[f"neutral.{p}"], spectra[f"selected.{p}"]
    sp = to_species(p)
    r = rec.setdefault(sp, dict(Pn=0., Ps=0., Dn=0., Ds=0., ne=[],
                                group=Populations.get_group_from_pop(p)))
    r["Pn"] += float(np.sum(sel.polymorphic)); r["Dn"] += float(sel.n_div)
    r["Ps"] += float(np.sum(neut.polymorphic)); r["Ds"] += float(neut.n_div)
    if p in ne_pop and np.isfinite(ne_pop[p]):
        r["ne"].append(ne_pop[p])

rows = []
for sp, r in rec.items():
    if not r["ne"]:
        continue
    Pn, Ps, Dn, Ds = r["Pn"], r["Ps"], r["Dn"], r["Ds"]
    NI = (Pn / Ps) / (Dn / Ds)
    se = np.sqrt(1 / Pn + 1 / Ps + 1 / Dn + 1 / Ds)
    rows.append(dict(species=sp, group=r["group"], Ne=float(np.mean(r["ne"])),
                     alpha=1 - NI, alpha_lo=1 - np.exp(np.log(NI) + 1.96 * se),
                     alpha_hi=1 - np.exp(np.log(NI) - 1.96 * se)))

df = pd.DataFrame(rows)
x = np.log10(df["Ne"].to_numpy())
y = df["alpha"].to_numpy()
lr = stats.linregress(x, y)
sr = stats.spearmanr(x, y)

fig, ax = plt.subplots(figsize=(8, 6))
for g in df["group"].unique():
    m = df["group"] == g
    ax.errorbar(df["Ne"][m], df["alpha"][m],
                yerr=[df["alpha"][m] - df["alpha_lo"][m], df["alpha_hi"][m] - df["alpha"][m]],
                fmt="o", ms=6, capsize=2, lw=0.7, color=Populations.get_color(g), label=g)
xs = np.linspace(x.min(), x.max(), 50)
ax.plot(10 ** xs, lr.intercept + lr.slope * xs, "k--", lw=1.3)
ax.set_xscale("log")
ax.axhline(0, color="grey", lw=0.6, ls=":")
ax.set_xlabel(r"$N_e$")
ax.set_ylabel(r"MK $\alpha = 1 - \mathrm{NI}$")
ax.set_title(f"Species-wise MK alpha vs Ne (catarrhini, n={len(df)})\n"
             f"Pearson r={lr.rvalue:.2f}, p={lr.pvalue:.1e}; "
             f"Spearman r={sr.correlation:.2f}, p={sr.pvalue:.1e}")
ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(out, dpi=180, bbox_inches="tight")
print("wrote", out)
print(df[["species", "group", "Ne", "alpha"]].round(4).sort_values("Ne").to_string(index=False))
print(f"\nPearson r={lr.rvalue:.3f} p={lr.pvalue:.2e} | Spearman r={sr.correlation:.3f} p={sr.pvalue:.2e}")
