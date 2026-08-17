"""
Extract the set of named genes from a GFF annotation.

Annotation sources differ in which attribute carries the gene symbol: GENCODE-style files use
``gene_name``, RefSeq-style files use ``gene`` (or ``Name``). Symbols are taken from the first of
these that is present, so gene sets are comparable across references. Unnamed model genes
(``LOC*``) are retained here and dropped downstream, since they are species-specific identifiers.
"""

import re

try:
    testing = False
    gff = snakemake.input.gff
    out = snakemake.output[0]
except NameError:
    testing = True
    gff = "resources/gff/Homo_sapiens.renamed.gff"
    out = "scratch/Homo_sapiens.genes.txt"

PATTERNS = [re.compile(p) for p in (r"gene_name=([^;]+)", r"gene=([^;]+)", r"Name=([^;]+)")]


def gene_name(attributes: str) -> str:
    """
    Return the gene symbol from a GFF attribute string, or an empty string if none is present.
    """
    for pattern in PATTERNS:
        match = pattern.search(attributes)
        if match:
            return match.group(1)

    return ""


genes = set()
with open(gff) as h:
    for line in h:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue

        name = gene_name(fields[8])
        if name:
            genes.add(name)

with open(out, "w") as h:
    h.write("\n".join(sorted(genes)) + "\n")

if testing:
    print(f"{len(genes)} named genes, {sum(1 for g in genes if not g.startswith('LOC'))} with a symbol")
