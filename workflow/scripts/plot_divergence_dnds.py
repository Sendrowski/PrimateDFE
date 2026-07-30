"""
Subspecies-level dN/dS and target-site ratio from the divergence counts.

For each population (subspecies), from the divergence summary:
  dN = n_div_sel / n_sites_div_sel        (0-fold substitutions per 0-fold site)
  dS = n_div_neut / n_sites_div_neut       (4-fold substitutions per 4-fold site)
  dN/dS                                    (lineage-specific, vs reconstructed ancestor)
  target ratio = n_sites_div_sel / n_sites_div_neut   (0-fold : 4-fold target sites)

Two horizontal-bar panels, one row per subspecies, coloured by taxon group.
"""
import glob
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from populations import Populations

try:
    summaries = list(snakemake.input.summaries)
    out = snakemake.output[0]
except NameError:
    summaries = sorted(glob.glob(
        "results/dfe/divergence/*/*/summary.gamma.div_on.8.json"))
    out = "results/graphs/divergence/catarrhini/dnds_target_ratio.8.png"

rows = []
for f in summaries:
    pop = f.split("/")[-2]
    d = json.load(open(f))
    dN = d["n_div_sel"] / d["n_sites_div_sel"]
    dS = d["n_div_neut"] / d["n_sites_div_neut"]
    rows.append(dict(
        population=pop,
        group=Populations.get_group_from_pop(pop),
        dNdS=dN / dS,
        target_ratio=d["n_sites_div_sel"] / d["n_sites_div_neut"],
    ))

df = pd.DataFrame(rows)
# order by taxon group (phylogenetic rank), then population name, so subspecies
# of a group sit together
df["rank"] = df["group"].map(Populations.get_label_rank)
df = df.sort_values(["rank", "population"], ascending=[False, False]).reset_index(drop=True)
colors = [Populations.get_color(g) for g in df["group"]]

y = np.arange(len(df))
fig, axes = plt.subplots(1, 2, figsize=(13, 0.32 * len(df) + 1.5), sharey=True)

axes[0].barh(y, df["dNdS"], color=colors, edgecolor="0.3", lw=0.3)
axes[0].set_xlabel(r"$d_N/d_S$  (0-fold / 4-fold, lineage-specific)")
axes[0].axvline(1.0, color="grey", lw=0.6, ls="--")

axes[1].barh(y, df["target_ratio"], color=colors, edgecolor="0.3", lw=0.3)
axes[1].set_xlabel("target-site ratio  (0-fold : 4-fold sites)")

axes[0].set_yticks(y)
axes[0].set_yticklabels(df["population"], fontsize=7)
axes[0].set_ylim(-0.6, len(df) - 0.4)

# legend by taxon group
seen = {}
for g, c in zip(df["group"], colors):
    seen.setdefault(g, c)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in seen.values()]
axes[1].legend(handles, list(seen.keys()), fontsize=7, title="taxon group",
               loc="lower right", frameon=False)

fig.suptitle("Subspecies-level divergence: dN/dS and target-site ratio (catarrhini)")
fig.tight_layout()
fig.savefig(out, dpi=180, bbox_inches="tight")
print("wrote", out)
print(df[["population", "group", "dNdS", "target_ratio"]].to_string(index=False))
