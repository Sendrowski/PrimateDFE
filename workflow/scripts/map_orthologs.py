"""
Assign each gene of a reference annotation to an ortholog group, so that gene sets are comparable
across references whose annotations use different identifiers.

NCBI GeneIDs are assigned per species and therefore do not match between annotations, and unnamed
RefSeq genes carry their ``LOC`` identifier as the symbol. Genes are instead anchored to their human
one-to-one ortholog, taken from NCBI ``gene_orthologs``. RefSeq annotations supply GeneIDs directly
through ``Dbxref``; GENCODE annotations supply Ensembl gene identifiers, which are joined to GeneIDs
through NCBI ``gene2ensembl``.
"""

import gzip
import re

import pandas as pd

try:
    testing = False
    gff = snakemake.input.gff
    orthologs = snakemake.input.orthologs
    gene2ensembl = snakemake.input.gene2ensembl
    out = snakemake.output[0]
except NameError:
    testing = True
    gff = "resources/gff/Macaca_nemestrina.gff.gz"
    orthologs = "resources/ncbi/gene_orthologs.gz"
    gene2ensembl = "resources/ncbi/gene2ensembl.gz"
    out = "scratch/Macaca_nemestrina.orthologs.txt"

HUMAN_TAXID = "9606"


def open_maybe_gzip(path: str):
    """
    Open a file that may or may not be gzipped.
    """
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_gene_anchors(path: str) -> dict:
    """
    Map each gene in the annotation to its own identifier: an NCBI GeneID where the annotation
    provides one, otherwise an Ensembl gene identifier.

    :param path: Path to the GFF annotation.
    :return: Gene interval key to identifier.
    """
    gene_id = re.compile(r"GeneID:(\d+)")
    ensembl = re.compile(r"gene_id=(ENS[A-Z]*G\d+)")

    anchors = {}
    with open_maybe_gzip(path) as h:
        for line in h:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            match = gene_id.search(fields[8]) or ensembl.search(fields[8])
            if match:
                anchors[(fields[0], int(fields[3]), int(fields[4]))] = match.group(1)

    return anchors


# human GeneID for every GeneID that has a one-to-one human ortholog, in both directions
pairs = pd.read_csv(orthologs, sep="\t", compression="gzip",
                    usecols=["#tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID"])
pairs = pairs[pairs.relationship == "Ortholog"]

to_human = {}
human_side = pairs[pairs["#tax_id"] == int(HUMAN_TAXID)]
to_human.update(dict(zip(human_side.Other_GeneID.astype(str), human_side.GeneID.astype(str))))
other_side = pairs[pairs.Other_tax_id == int(HUMAN_TAXID)]
to_human.update(dict(zip(other_side.GeneID.astype(str), other_side.Other_GeneID.astype(str))))

# the human reference anchors to itself, so its own GeneIDs are valid groups
human_ids = set(to_human.values())

# Ensembl gene identifier to GeneID, needed for GENCODE-style annotations
ens = pd.read_csv(gene2ensembl, sep="\t", compression="gzip",
                  usecols=["GeneID", "Ensembl_gene_identifier"])
ens_to_gene = dict(zip(ens.Ensembl_gene_identifier.astype(str), ens.GeneID.astype(str)))

rows = []
for (contig, start, end), anchor in read_gene_anchors(gff).items():
    if anchor.startswith("ENS"):
        anchor = ens_to_gene.get(anchor)

    group = to_human.get(anchor) or (anchor if anchor in human_ids else None)
    if group:
        rows.append(dict(contig=contig, start=start, end=end, group=group))

genes = pd.DataFrame(rows).drop_duplicates()
genes.to_csv(out, sep="\t", index=False)

if testing:
    print(f"{len(genes)} genes assigned to {genes.group.nunique()} ortholog groups")
