"""
Plot regression: Ne (Kuderna et al.) vs Ne (Watterson)
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import LogLocator, FuncFormatter

from populations import Populations
from regression import Regression
from visualization import DFEvsNePlotter

try:
    testing = False
    x_file = snakemake.input.x
    y_file = snakemake.input.y
    tree_file = snakemake.input.tree
    label_dict = snakemake.params.get("labels", None)
    out = snakemake.output[0]
except NameError:
    testing = True
    x_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    y_file = "resources/Kuderna/species.csv"
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    label_dict = {p: Populations.get_group_from_pop(p) for p in Populations.get_pops(8, "catarrhini")}
    out = "scratch/Ne_reg.png"

# -------------------- load --------------------

x_df = pd.read_csv(x_file)
x_df["display_name"] = x_df["label"].map(label_dict)  # label -> group
y_df = pd.read_csv(y_file)

merged = x_df.merge(
    y_df,
    left_on="label",
    right_on="SPECIES_BINOMIAL",
    how="inner",
)

merged = merged[merged.label.isin(label_dict.keys())]

labels = merged["display_name"].tolist()

label_to_color = Populations.get_label_color_map(labels)  # expects list of labels (groups)

x = merged["EFFECTIVE_POP_SIZE"]  # Kuderna Ne
y = merged["x"]  # Watterson Ne

intercept, slope, p, r = Regression.regress(
    'phylo',
    x=np.log10(x),
    y=np.log10(y),
    pops=merged["label"].to_numpy(),
    tree_file=tree_file
)

plt.figure(figsize=(5, 5), dpi=400)

plt.scatter(
    x,
    y,
    s=70,
    alpha=0.8,
    color=[label_to_color[l] for l in merged["display_name"]],
)

order = np.argsort(x.to_numpy())
x_sorted = x.iloc[order]
y_sorted = 10 ** (intercept + slope * np.log10(x_sorted))

reg_line, = plt.plot(x_sorted, y_sorted, "--", color="black", label=f"r={r:.2f}, p={p:.3g}")

handles = [
    plt.Line2D(
        [0], [0],
        marker="o",
        linestyle="",
        color=color,
        label=label,
    )
    for label, color in label_to_color.items()
]
handles.append(reg_line)

plt.legend(handles, [Populations.label_to_text(h.get_label()) for h in handles], fontsize=9)

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_aspect("equal", adjustable="box")

# y = x reference line
xmin = min(ax.get_xlim()[0], ax.get_ylim()[0])
xmax = max(ax.get_xlim()[1], ax.get_ylim()[1])

xx = np.logspace(np.log10(xmin), np.log10(xmax), 200)
plt.plot(xx, xx, ":", color="black", linewidth=1, label="$y=x$")

plt.xlabel("$N_e$ (Kuderna et al.)")
plt.ylabel("$N_e$ (Watterson)")
plt.title("$N_e$ (Kuderna et al.) vs $N_e$ (Watterson)")
plt.tight_layout()

ax.set_xlim(1e4, max(x) * 1.1)
ax.set_ylim(1e4, max(x) * 1.1)

ax.set_xscale("log")
for ax in [ax.xaxis, ax.yaxis]:
    ax.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax.set_major_formatter(FuncFormatter(DFEvsNePlotter.log_label_pow))

plt.savefig(out)

if testing:
    plt.show()
