"""
Write the species order used by the phylogeny figure, so that other figures can be laid out
in the same order without depending on the tree environment.
"""

import pandas as pd
import toytree

from populations import Populations

try:
    testing = False
    tree_file = snakemake.input.tree
    pops = snakemake.params.pops
    out = snakemake.output[0]
except NameError:
    testing = True
    tree_file = "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
    pops = Populations.get_pops(8, "catarrhini")
    out = "results/tables/phylogeny/catarrhini/species_order.8.csv"


def to_species(pop: str) -> str:
    return "_".join(pop.split("_")[:2])


species = {to_species(p) for p in pops}

tree = toytree.tree(tree_file)
keep = [tip for tip in tree.get_tip_labels() if to_species(tip) in species]
tree = toytree.mod.prune(tree, *keep)
tree = toytree.mod.ladderize(tree, direction=False)

order = [to_species(tip) for tip in tree.get_tip_labels()]
missing = sorted(species - set(order))

pd.DataFrame(dict(rank=range(len(order) + len(missing)), species=order + missing)).to_csv(out, index=False)

if testing:
    print(f"{len(order)} species from tree, {len(missing)} appended: {missing}")
