"""
Plot a fossil-calibrated species tree using toytree.
"""
import matplotlib.pyplot as plt
import toyplot.png
import toytree
from matplotlib.lines import Line2D

from populations import Populations

try:
    testing = False
    tree_file = snakemake.input.tree
    pops = snakemake.params.populations
    out = snakemake.output[0]
except NameError:
    testing = True
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    #tree_file = "scratch/dfe_clustering.nwk"
    pops = Populations.get_pops(8, "catarrhini")
    # pops = ['Macaca_mulatta', 'Papio_anubis', 'Homo_sapiens', 'Semnopithecus_hypoleucos']
    out = "scratch/species_tree.png"

grey = '#666666'

# load tree (Newick)
tree = toytree.tree(tree_file)

# --- filter to populations ---
# keep only taxa present in both tree and pops
tips = set(tree.get_tip_labels())
keep = [p for p in pops if p in tips]

tree = toytree.mod.prune(tree, *keep)

tree = toytree.mod.ladderize(tree, direction=False)

# after pruning / ladderizing
nnodes = tree.nnodes
node_colors = ["#999999"] * nnodes
node_mask = [False] * nnodes

# map tip label -> color
label_to_color = {
    p: Populations.get_color(p)
    for p in tree.get_tip_labels()
}

# postorder: tips → internals
for node in tree.treenode.traverse("postorder"):
    idx = node.idx

    if node.is_leaf():
        node_colors[idx] = label_to_color[node.name]
        node_mask[idx] = True
    else:
        node_mask[idx] = False

# draw
canvas, axes, _ = tree.draw(
    width=1200,
    height=600,
    tip_labels=[l.replace("_", " ") for l in tree.get_tip_labels()],
    tip_labels_align=True,
    node_sizes=18,
    node_colors=node_colors,
    node_mask=node_mask,
    edge_style={"stroke": 'black', "stroke-width": 2},
    scale_bar=True,
    tip_labels_style={
        "font-size": 20,
        "fill": "black"
    },
)

axes.x.spine.style["stroke-width"] = 3
axes.x.ticks.style["stroke-width"] = 3
axes.x.ticks.labels.style["font-size"] = 20

scale = 2
toyplot.png.render(
    canvas,
    out,
    scale=scale
)

img = plt.imread(out)
h, w = img.shape[:2]

dpi = 500
fig = plt.figure(figsize=(w / dpi, h / dpi), dpi=dpi)
plt.imshow(img)
plt.axis("off")

# build legend entries (unique labels present)
labels_present = sorted(
    {Populations.get_group_from_pop(p) for p in keep},
    key=Populations.get_label_rank,
)

handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=1.8 * scale,
        markerfacecolor=Populations.get_label_color_map(labels_present)[label],
        markeredgecolor='#444444',  # ← edge color
        markeredgewidth=0.2 * scale,  # ← thin edge
        label=Populations.label_to_text(label),
    )
    for label in labels_present
]

# add x axis label as text
plt.text(
    0.5, -0.02, "Time (Ma)",
    ha='center', va='center',
    transform=plt.gca().transAxes,
    fontsize=6,
)

plt.legend(
    handles=handles,
    loc=(0.05, 0.65),
    frameon=True,
    framealpha=1.0,
    edgecolor="grey",
    fontsize=5,
    labelspacing=0.2,  # vertical spacing between entries
    handletextpad=0.3,  # space between marker and text
    borderpad=0.5,  # padding inside frame
)
plt.gca().get_legend().get_frame().set_linewidth(0.2 * scale)

plt.tight_layout()
plt.subplots_adjust(left=0, right=1.1, top=1, bottom=0.05)
plt.savefig(out, dpi=dpi)

if testing:
    plt.show()
