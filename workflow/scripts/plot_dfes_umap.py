"""
UMAP on unscaled DFEs using PDF mass in many s-ranges via Parametrizer.get_S_range.
"""
import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap

from parametrizer import Parametrizer
from populations import Populations

# ---- inputs ----
dfe_file = "results/dfe/catarrhini/dfe.unfolded.8.discrete.full.noeps.csv"
ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
out = "scratch/umap_unscaled_mass_ranges.png"

# ---- load DFEs ----
df = pd.read_csv(dfe_file)
dfes = {row.population: fd.DFE.from_json(row.json) for _, row in df.iterrows()}
pops = list(dfes.keys())

ne_df = pd.read_csv(ne_file)
ne_dict = dict(zip(ne_df["label"], ne_df["x"]))

# ---- define MANY unscaled s-ranges (edit as you like) ----
# deleterious: strong -> weak, then near-neutral, then beneficial weak -> strong
ranges = [
    (-np.inf, -10),
    (-10, -1),
    (-1, 0),
    (0, np.inf)
]
ranges = [
    (-np.inf, -1e-3),
    (-1e-3, -1e-5),
    (-1e-5, 0),
    (0, np.inf)
]

range_names = [f"mass_s_{lo:g}_{hi:g}" for lo, hi in ranges]

# ---- feature matrix: median bootstrap mass in each range ----
X = np.zeros((len(pops), len(ranges)), dtype=float)
labels = []

for i, p in enumerate(pops):
    dfe = dfes[p]
    Ne = ne_dict[p]
    for j, (lo, hi) in enumerate(ranges):
        mass_bs = Parametrizer.get_S_range(dfe, 4 * Ne * lo, 4 * Ne * hi)
        # mass_bs = Parametrizer.get_S_range(dfe, lo, hi)  # unscaled
        X[i, j] = np.median(mass_bs)
    labels.append(Populations.get_group_from_pop(p))

# optional: stabilize (UMAP often benefits from this)
X = np.clip(X, 1e-12, 1.0)
X = np.log(X)

# ---- UMAP ----
emb = umap.UMAP(
    n_neighbors=12,
    min_dist=0.1,
    n_components=2,
    metric="euclidean",
    random_state=0,
).fit_transform(X)

# ---- plot ----
plt.figure(figsize=(4.2, 4.2), dpi=300)

for lab in sorted(set(labels)):
    m = np.array(labels) == lab
    plt.scatter(emb[m, 0], emb[m, 1], s=35, alpha=0.85, label=lab)

plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.legend(fontsize=7, frameon=False, ncol=2)
plt.tight_layout()
plt.savefig(out)
plt.show()
