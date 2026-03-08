"""
Plot DFEs as bars, with a tree on top.
"""

import fastdfe as fd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import toyplot.png
import toytree
from PIL import Image
from PIL import ImageDraw, ImageFont
from PIL import ImageOps
from fastdfe.visualization import Visualization
from matplotlib.collections import LineCollection
from matplotlib.container import BarContainer
from matplotlib.lines import Line2D

from populations import Populations

try:
    testing = False
    tree_file = snakemake.input.tree
    dfe_file = snakemake.input.dfe
    pops = snakemake.params.populations
    legend = snakemake.params.legend
    out_tree = snakemake.output.tree
    out_dfe = snakemake.output.dfe
    out_full = snakemake.output.full
except NameError:
    testing = True
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    dfe_file = "results/dfe/catarrhini/dfe.unfolded.8.gamma.full.noeps.csv"
    pops = Populations.get_pops(8, "catarrhini")
    legend = True
    out_tree = "scratch/tree.png"
    out_dfe = "scratch/dfe.png"
    out_full = "scratch/tree_plus_dfe.png"


def darken(color, factor=0.7):
    """
    Darken a color by multiplying its RGB values by a factor.
    """
    r, g, b = mcolors.to_rgb(color)
    return r * factor, g * factor, b * factor


def to_species(name: str) -> str:
    """
    Convert a population name to a species name by removing the last underscore and suffix.
    """
    return name.rsplit("_", 1)[0] if name.count("_") >= 2 else name


grey = "#666666"

# ----------------------------
# Load DFEs, collapse subspecies → species, discretize + average
# ----------------------------
df = pd.read_csv(dfe_file)

df["species"] = df["population"].map(to_species)

# ----------------------------
# Load tree (prune at species level)
# ----------------------------
tree = toytree.tree(tree_file)

# species present in DFEs
species_in_dfe = set(df["species"])

# keep tree tips whose species is present in DFEs
keep_species = [
    tip for tip in tree.get_tip_labels()
    if to_species(tip) in species_in_dfe
]

tree = toytree.mod.prune(tree, *keep_species)
tree = toytree.mod.ladderize(tree, direction=False)

# ----------------------------
# Load DFEs, collapse subspecies → species, discretize + average
# ----------------------------
intervals = np.array([-np.inf, -10, -1, 0, np.inf])

# species -> list of discretized results
species_vals = {}
species_errs = {}

for sp, g in df.groupby("species"):
    vals_list, lo_list, hi_list = [], [], []
    for j in g["json"]:
        d = fd.DFE.from_json(j)
        vals, errs = d.discretize(
            bins=intervals,
            confidence_intervals=True,
            ci_level=0.05,
            bootstrap_type="percentile",
            point_estimate="median",
        )

        vals_list.append(np.asarray(vals))
        lo_list.append(np.asarray(errs[0]))
        hi_list.append(np.asarray(errs[1]))

    V = np.vstack(vals_list)
    LO = np.vstack(lo_list)
    HI = np.vstack(hi_list)

    # mean across subspecies
    v_mean = V.mean(axis=0)
    lo_mean = LO.mean(axis=0)
    hi_mean = HI.mean(axis=0)

    species_vals[sp] = v_mean
    species_errs[sp] = np.vstack([lo_mean, hi_mean])  # (2, nbins) works with Visualization

# ----------------------------
# Canonical order from tree (species-level)
# ----------------------------
tree_order = [to_species(p) for p in tree.get_tip_labels()]

values = [species_vals[p] for p in tree_order]
errors = [species_errs[p] for p in tree_order]
labels_ordered = [Populations.get_group_from_pop(p) for p in tree_order]

# ----------------------------
# ---- DRAW TREE (toytree → PNG)
# ----------------------------
nnodes = tree.nnodes
node_colors = ["#999999"] * nnodes
node_mask = [False] * nnodes

label_to_color = {p: Populations.get_color(p) for p in tree.get_tip_labels()}

for node in tree.treenode.traverse("postorder"):
    if node.is_leaf():
        node_colors[node.idx] = label_to_color[node.name]
        node_mask[node.idx] = True

canvas, axes, _ = tree.draw(
    width=1600,
    height=600,
    tip_labels=[l.replace("_", " ") for l in tree.get_tip_labels()],
    tip_labels_align=True,
    node_sizes=15,
    node_colors=node_colors,
    node_mask=node_mask,
    edge_style={"stroke": 'black', "stroke-width": 2},
    scale_bar=True,
    tip_labels_style={
        "font-size": 16,
        "fill": "black"
    },
)

axes.x.spine.style["stroke-width"] = 3
axes.x.ticks.style["stroke-width"] = 3
axes.x.ticks.labels.style["font-size"] = 20

toyplot.png.render(canvas, out_tree, scale=2)

tree_img = Image.open(out_tree).convert("RGBA")

# crop using alpha channel (removes empty margins correctly)
alpha = tree_img.split()[-1]
bbox = alpha.getbbox()
tree_img = tree_img.crop(bbox)

# put on white background
bg = Image.new("RGBA", tree_img.size, (255, 255, 255, 255))
tree_img = Image.alpha_composite(bg, tree_img)

# add padding
tree_img = ImageOps.expand(tree_img, border=(180, 0, 0, 0), fill='white')

tree_img.save(out_tree)
tree_img = np.array(tree_img)

# ---- TREE FIGURE ----
fig_tree = plt.figure(figsize=(10, 5))
ax_tree = fig_tree.add_subplot(111)
ax_tree.imshow(tree_img)
ax_tree.axis("off")

fig_tree.savefig(out_tree, dpi=300, bbox_inches="tight")
plt.close(fig_tree)

# ---- DFE FIGURE ----
fig_dfe, ax_dfe = plt.subplots(figsize=(10, 2.5))

Visualization.plot_discretized(
    ax=ax_dfe,
    values=values,
    errors=errors,
    labels=labels_ordered,
    intervals=intervals,
    title="",
    show=False,
    kwargs_legend=dict(prop=dict(size=8)),
)

colors = [Populations.get_color(p) for p in tree_order]
bars = [c for c in ax_dfe.containers if isinstance(c, BarContainer)]

for container, color in zip(bars, colors):
    edge = darken(color, 0.7)
    for bar in container:
        bar.set_facecolor(color)
        bar.set_edgecolor(edge)
        bar.set_linewidth(0.6)

    eb = getattr(container, "errorbar", None)
    if eb:
        for cap in eb.lines[1]:
            cap.set_color(edge)
            cap.set_linewidth(2)

        lc = eb.lines[2][0]
        if isinstance(lc, LineCollection):
            lc.set_edgecolor(edge)
            lc.set_linewidth(2)

ax_dfe.get_legend().remove()

labels = []
for lab in ax_dfe.get_xticklabels():
    t = lab.get_text()
    t = t.replace("inf", "∞")
    labels.append(t)

labels = [t.get_text().replace("inf", "∞") for t in ax_dfe.get_xticklabels()]
ax_dfe.set_xticklabels(labels)
ax_dfe.set_xlabel("$S=4N_es$")
ax_dfe.set_ylabel("probability mass")

ax_dfe.spines["top"].set_visible(False)
ax_dfe.spines["right"].set_visible(False)

# ----------------------------
# Shared legend (groups)
# ----------------------------
labels_present = sorted(
    {Populations.get_group_from_pop(p) for p in tree_order},
    key=Populations.get_label_rank,
)

handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=6,
        markerfacecolor=Populations.get_label_color_map(labels_present)[lab],
        markeredgecolor="#444444",
        markeredgewidth=0.5,
        label=Populations.label_to_text(lab),
    )
    for lab in labels_present
]

leg = ax_dfe.legend(
    handles=handles,
    loc=(0.8, 0.67),
    frameon=True,
    fontsize=7,
)

fig_dfe.savefig(out_dfe, dpi=300, bbox_inches="tight")
plt.close(fig_dfe)

img1 = Image.open(out_tree)
img2 = Image.open(out_dfe)

# target width = minimum width
w = min(img1.width, img2.width)


def resize_to_width(img, w):
    h = int(img.height * (w / img.width))
    return img.resize((w, h), Image.LANCZOS)


img1 = resize_to_width(img1, w)
img2 = resize_to_width(img2, w)

h = img1.height + img2.height

canvas = Image.new("RGBA", (w, h), (255, 255, 255, 255))
canvas.paste(img1, (0, 0))
canvas.paste(img2, (0, img1.height))

draw = ImageDraw.Draw(canvas)

text = "Time (Ma)"
font = ImageFont.truetype("Arial.ttf", 40)

# center under the tree
tw, th = draw.textbbox((0, 0), text, font=font)[2:]
x = (canvas.width - tw) // 2
y = img1.height - 10  # slightly above boundary to DFE panel

draw.text((x, y), text, fill="black", font=font)

canvas.save(out_full)
canvas.show()
