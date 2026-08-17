"""
Quantify how far the gene sets of the reference annotations overlap, and write the set of genes
named in every reference.

Genes are compared through the ortholog groups assigned by :mod:`map_orthologs`, since gene symbols
and GeneIDs are not comparable across annotations.
"""

import itertools
import os

import pandas as pd

try:
    testing = False
    gene_files = snakemake.input.genes
    refs = snakemake.params.refs
    out_overlap = snakemake.output.overlap
    out_shared = snakemake.output.shared
except NameError:
    testing = True
    refs = ["Homo_sapiens", "Pan_troglodytes"]
    gene_files = [f"results/tables/orthologs/{r}.genes.tsv" for r in refs]
    out_overlap = "results/tables/orthologs/catarrhini/overlap.8.csv"
    out_shared = "results/tables/orthologs/catarrhini/shared_groups.8.txt"

sets = {}
for ref, path in zip(refs, gene_files):
    sets[ref] = set(pd.read_csv(path, sep="\t", dtype=dict(group=str)).group)

rows = []
for a, b in itertools.combinations(sets, 2):
    shared = len(sets[a] & sets[b])
    rows.append(dict(
        reference_a=a,
        reference_b=b,
        genes_a=len(sets[a]),
        genes_b=len(sets[b]),
        shared=shared,
        fraction_of_smaller=shared / min(len(sets[a]), len(sets[b])),
        jaccard=shared / len(sets[a] | sets[b]),
    ))

overlap = pd.DataFrame(rows)
overlap.to_csv(out_overlap, index=False)

shared = set.intersection(*sets.values())
with open(out_shared, "w") as h:
    h.write("\n".join(sorted(shared)) + "\n")

if testing:
    print(overlap.to_string(index=False))
    print(f"genes named in all {len(sets)} references: {len(shared)}")
