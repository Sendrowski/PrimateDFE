"""
Plot Ne against pN/pS with regression.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress
from matplotlib.ticker import LogLocator, FuncFormatter, ScalarFormatter
from visualization import DFEvsNePlotter  # for log_label_pow
from populations import Populations
from regression import Regression

try:
    testing = False
    counts_file = snakemake.input[0]
    x_file = snakemake.input.x
    y_file = snakemake.input.y
    tree_file = snakemake.input.tree
    labels = snakemake.params.labels
    populations = snakemake.params.populations
    include = snakemake.params.include
    title = snakemake.params.title
    out = snakemake.output[0]
except NameError:
    testing = True
    x_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    y_file = "results/stats/pNpS/comp/original_ref/catarrhini/8.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    include = None
    labels = None
    populations = Populations.get_pops(8, "catarrhini")
    title = "$N_e$ vs $p_N/p_S$"
    out = "scratch/Ne_ns.png"

# -------------------- load --------------------

def to_species(name: str) -> str:
    """
    Convert a population name to a species name by removing the last underscore and suffix.
    """
    return name.rsplit("_", 1)[0] if name.count("_") >= 2 else name

x_df = pd.read_csv(x_file)
y_df = pd.read_csv(y_file)

# merge on population label
df = x_df.merge(
    y_df,
    left_on="label",
    right_on="label",
    how="inner"
)

df = df[df["label"].isin(populations)]

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

x = df["x_x"].to_numpy()
y = df["x_y"].to_numpy()

mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
x = x[mask]
y = y[mask]
pops = df.loc[mask, "label"].to_numpy()
groups = df.loc[mask, "group"].to_numpy()

intercept, slope, p, r = Regression.regress(
    "phylo",
    x=np.log10(x),
    y=np.log10(y),
    pops=pops,
    tree_file=tree_file,
)

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
y_fit = 10 ** (intercept + slope * np.log10(x_sorted))
plt.plot(
    x_sorted,
    y_fit,
    "--",
    color="black",
    label=f"r={r:.2f}, p={p:.3g}",
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
               label=f"r={r:.2f}, p={p:.3g}")
)

plt.legend(handles=handles, fontsize=9)

plt.xlabel(r"$N_e$")
plt.ylabel(r"$p_N/p_S$")
plt.title(title)
plt.tight_layout()

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")

ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
ax.xaxis.set_major_formatter(FuncFormatter(DFEvsNePlotter.log_label_pow))

# y-axis: plain decimals
ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.get_major_formatter().set_scientific(False)
ax.yaxis.get_major_formatter().set_useOffset(False)

plt.savefig(out)

if testing:
    plt.show()
