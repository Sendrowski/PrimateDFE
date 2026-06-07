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

# load tree (Newick)
tree = toytree.tree(tree_file)

# --- filter to populations ---
# keep only taxa present in both tree and pops
tips = set(tree.get_tip_labels())
keep = [p for p in pops if p in tips]
not_kept = set(pops) - set(keep)

tree = toytree.mod.prune(tree, *keep)

pass