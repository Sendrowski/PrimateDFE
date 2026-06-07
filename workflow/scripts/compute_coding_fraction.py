"""
Compute the fraction of the reference genome that is annotated as coding (CDS).

Merges overlapping/duplicate CDS intervals per contig (so the numerator is the
union of CDS bases, not their multiset sum), restricts to contigs present in
the FASTA index, and divides by the total reference length.
"""

import gzip

import pandas as pd

try:
    gff = snakemake.input.gff
    fai = snakemake.input.fai
    ref = snakemake.wildcards.ref
    out = snakemake.output[0]
except NameError:
    # testing
    ref = "Homo_sapiens"
    gff = f"resources/gff/{ref}.gff.gz"
    fai = f"resources/ref/{ref}.fasta.fai"
    out = f"scratch/coding_fraction.{ref}.csv"


fai_df = pd.read_csv(
    fai,
    sep="\t",
    header=None,
    names=["contig", "length", "offset", "linebases", "linewidth"],
    usecols=["contig", "length"],
)
genome_bp = int(fai_df["length"].sum())
valid_contigs = set(fai_df["contig"])

intervals: dict[str, list[tuple[int, int]]] = {}
opener = gzip.open if gff.endswith(".gz") else open
with opener(gff, "rt") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 5 or cols[2] != "CDS":
            continue
        contig = cols[0]
        if contig not in valid_contigs:
            continue
        start = int(cols[3]) - 1
        end = int(cols[4])
        if end <= start:
            continue
        intervals.setdefault(contig, []).append((start, end))

cds_bp = 0
for ivs in intervals.values():
    ivs.sort()
    cur_end = -1
    for s, e in ivs:
        if s >= cur_end:
            cds_bp += e - s
            cur_end = e
        elif e > cur_end:
            cds_bp += e - cur_end
            cur_end = e

pd.DataFrame(
    [{
        "ref": ref,
        "cds_bp": cds_bp,
        "genome_bp": genome_bp,
        "fraction_coding": cds_bp / genome_bp if genome_bp else float("nan"),
    }]
).to_csv(out, index=False)
