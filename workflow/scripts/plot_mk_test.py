"""
Species-wise McDonald-Kreitman test from the divergence-carrying SFS.

For each population the 2x2 MK table is
    [[Pn, Dn],        Pn = nonsyn (0-fold) polymorphisms,  Dn = 0-fold substitutions
     [Ps, Ds]]        Ps = syn (4-fold) polymorphisms,     Ds = 4-fold substitutions
read directly from the spectra (polymorphic bins summed; last bin = divergence).
Subspecies are pooled to species (counts summed). alpha = 1 - NI with
NI = (Pn/Ps)/(Dn/Ds) (the table odds ratio); CI and p via the normal
approximation to ln(NI) (SE = sqrt(1/Pn+1/Ps+1/Dn+1/Ds)).
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
    out = snakemake.output[0]
except NameError:
    sfs_file = "results/divergence/sfs/comp/catarrhini/sfs.8.csv"
    out = "results/graphs/divergence/catarrhini/mk_test.8.png"


def to_species(pop):
    # drop subspecies suffix when there are >=2 underscores (mirrors the PGLS aggregation)
    return pop.rsplit("_", 1)[0] if pop.count("_") >= 2 else pop


spectra = fd.Spectra.from_file(sfs_file)
pops = sorted({t.split(".", 1)[1] for t in spectra.types if t.startswith("neutral.")})

# per-population MK counts
rec = {}
for p in pops:
    neut, sel = spectra[f"neutral.{p}"], spectra[f"selected.{p}"]
    sp = to_species(p)
    Ps, Ds = float(np.sum(neut.polymorphic)), float(neut.n_div)
    Pn, Dn = float(np.sum(sel.polymorphic)), float(sel.n_div)
    r = rec.setdefault(sp, dict(Pn=0., Ps=0., Dn=0., Ds=0., pops=0,
                                group=Populations.get_group_from_pop(p)))
    r["Pn"] += Pn; r["Ps"] += Ps; r["Dn"] += Dn; r["Ds"] += Ds; r["pops"] += 1

rows = []
for sp, r in rec.items():
    Pn, Ps, Dn, Ds = r["Pn"], r["Ps"], r["Dn"], r["Ds"]
    NI = (Pn / Ps) / (Dn / Ds)
    alpha = 1 - NI
    se = np.sqrt(1 / Pn + 1 / Ps + 1 / Dn + 1 / Ds)   # SE of ln(NI)
    z = np.log(NI) / se
    p = 2 * stats.norm.sf(abs(z))
    ni_lo, ni_hi = np.exp(np.log(NI) - 1.96 * se), np.exp(np.log(NI) + 1.96 * se)
    rows.append(dict(species=sp, group=r["group"], n_sub=r["pops"],
                     Pn=Pn, Ps=Ps, Dn=Dn, Ds=Ds, NI=NI, alpha=alpha,
                     alpha_lo=1 - ni_hi, alpha_hi=1 - ni_lo, p=p))

df = pd.DataFrame(rows)
df["rank"] = df["group"].map(Populations.get_label_rank)
df = df.sort_values(["rank", "species"], ascending=[False, False]).reset_index(drop=True)


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


y = np.arange(len(df))
colors = [Populations.get_color(g) for g in df["group"]]
fig, ax = plt.subplots(figsize=(9, 0.34 * len(df) + 1.5))
xerr = np.vstack([df["alpha"] - df["alpha_lo"], df["alpha_hi"] - df["alpha"]])
ax.barh(y, df["alpha"], xerr=xerr, color=colors, edgecolor="0.3", lw=0.3,
        error_kw=dict(lw=0.7, capsize=2))
ax.axvline(0, color="grey", lw=0.7, ls="--")
ax.set_yticks(y); ax.set_yticklabels(df["species"], fontsize=7)
ax.set_ylim(-0.6, len(df) - 0.4)
ax.set_xlabel(r"MK $\alpha = 1 - \mathrm{NI}$")
xmax = df["alpha_hi"].max()
for yi, (a, hi, pv) in enumerate(zip(df["alpha"], df["alpha_hi"], df["p"])):
    ax.text(max(hi, a) + 0.01, yi, stars(pv), va="center", fontsize=6.5, color="0.3")

seen = {}
for g, c in zip(df["group"], colors):
    seen.setdefault(g, c)
ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in seen.values()],
          list(seen.keys()), fontsize=7, title="taxon group", loc="lower right", frameon=False)
ax.set_title("Species-wise McDonald-Kreitman test (catarrhini)\n"
             "pooled subspecies; * p<0.05 ** p<0.01 *** p<0.001 (NI≠1)", fontsize=10)
fig.tight_layout()
fig.savefig(out, dpi=180, bbox_inches="tight")
print("wrote", out)
print(df[["species", "n_sub", "alpha", "NI", "p"]].round(4).to_string(index=False))
