"""
Plot Ne against N/S with regression.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

from populations import Populations

try:
    testing = False
    counts_file = snakemake.input[0]
    x_file = snakemake.input.x
    y_file = snakemake.input.y
    labels = snakemake.params.labels
    include = snakemake.params.include
    title = snakemake.params.title
    out = snakemake.output[0]
except NameError:
    testing = True
    x_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    y_file = "results/stats/comp/ns/ingroup/population/NS.csv"
    include = None
    labels = None
    title = "$N_e$ vs $N/S$"
    out = "scratch/Ne_ns.png"

# -------------------- load --------------------

x_df = pd.read_csv(x_file)
y_df = pd.read_csv(y_file)

# merge on population label
df = x_df.merge(
    y_df,
    left_on="label",
    right_on="population",
    how="inner",
)

# optional filtering + ordering
if include is not None:
    df = df.set_index("label").loc[include].reset_index()

# group labels (e.g. great_apes, macaques, …)
if labels is not None:
    df["group"] = labels
else:
    df["group"] = df["label"].apply(Populations.get_group_from_pop)

# -------------------- colors & order --------------------

label_color = Populations.get_label_color_map(df["group"])
label_order = sorted(set(df["group"]), key=Populations.get_label_rank)

# -------------------- regression --------------------

x = df["x"]        # Ne
y = df["ns"]     # N/S

slope, intercept, r, pval, _ = linregress(x, y)

# -------------------- plot --------------------

plt.figure(figsize=(6, 4))

plt.scatter(
    x,
    y,
    color=[label_color[g] for g in df["group"]],
    s=60,
    alpha=0.8,
)

x_sorted = np.sort(x)
plt.plot(
    x_sorted,
    intercept + slope * x_sorted,
    "--",
    color="black",
    label=f"r={r:.2f}, p={pval:.3g}",
)

# legend
handles = [
    plt.Line2D(
        [0], [0],
        marker="o",
        linestyle="",
        color=label_color[l],
        label=Populations.label_to_text(l),
        markersize=6,
    )
    for l in label_order
]

handles.append(
    plt.Line2D([0], [0], linestyle="--", color="black",
               label=f"r={r:.2f}, p={pval:.3g}")
)

plt.legend(handles=handles, fontsize=9)

plt.xlabel(r"$N_e$")
plt.ylabel(r"$N/S$")
plt.title(title)
plt.tight_layout()

plt.savefig(out)

if testing:
    plt.show()
