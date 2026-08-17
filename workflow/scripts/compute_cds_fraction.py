"""
Fraction of coding sequence lying inside genes shared across all reference annotations.

The shared gene set covers only part of each reference's coding sequence, so spectra computed on it
must be given a correspondingly smaller mutational target. Coding intervals are merged before
measurement, since overlapping isoforms otherwise count the same bases repeatedly.
"""

import gzip
from bisect import bisect_right
from typing import Dict, List, Tuple

import pandas as pd

try:
    testing = False
    gffs = snakemake.input.gffs
    gene_files = snakemake.input.genes
    shared_file = snakemake.input.shared
    refs = snakemake.params.refs
    out = snakemake.output[0]
except NameError:
    testing = True
    refs = ["Homo_sapiens"]
    gffs = [f"resources/gff/{r}.gff.gz" for r in refs]
    gene_files = [f"results/tables/orthologs/{r}.genes.tsv" for r in refs]
    shared_file = "results/tables/orthologs/catarrhini/shared_groups.8.txt"
    out = "results/tables/orthologs/cds_fraction.csv"


def merge(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping or adjacent intervals.
    """
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    return merged


def index(genes: pd.DataFrame) -> Dict[str, Tuple[List[int], List[int]]]:
    """
    Merged gene intervals per contig, as parallel start and end lists for binary search.
    """
    out = {}
    for contig, group in genes.groupby("contig"):
        merged = merge(list(zip(group.start, group.end)))
        out[contig] = ([s for s, _ in merged], [e for _, e in merged])

    return out


shared = {line.strip() for line in open(shared_file) if line.strip()}

fractions = {}
for ref, gff, gene_file in zip(refs, gffs, gene_files):
    genes = pd.read_csv(gene_file, sep="\t", dtype=dict(group=str))
    genes = genes[genes.group.isin(shared)]
    idx = index(genes)

    cds: Dict[str, List[Tuple[int, int]]] = {}
    with gzip.open(gff, "rt") as h:
        for line in h:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue

            cds.setdefault(fields[0], []).append((int(fields[3]), int(fields[4])))

    total = kept = 0
    for contig, intervals in cds.items():
        starts, ends = idx.get(contig, ([], []))
        for start, end in merge(intervals):
            total += end - start + 1
            i = bisect_right(starts, start) - 1
            if i >= 0 and start <= ends[i]:
                kept += end - start + 1

    fractions[ref] = kept / total if total else 0.0
    if testing:
        print(f"{ref:26s} {fractions[ref]:.4f}")

pd.Series(fractions, name="cds_fraction_in_shared_genes").to_csv(out)
