"""
Annotate site degeneracy.
"""

import fastdfe as fd
import jsonpickle
import pandas as pd

try:
    testing = False
    vcf_in = snakemake.input.vcf
    gff = snakemake.input.gff
    fasta = snakemake.input.fasta
    vcf_out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf_in = "results/vcf/Pan_troglodytes/Pan_troglodytes.vcf.gz"
    gff = "resources/gff/Pan_troglodytes.genbank.gff.gz"
    fasta = "resources/ref/Pan_troglodytes.fasta"
    vcf_out = "scratch/chr22.degeneracy.vcf.gz"

ann = fd.Annotator(
    annotations=[fd.DegeneracyAnnotation()],
    vcf=vcf_in,
    gff=gff,
    fasta=fasta,
    output=vcf_out
)

ann.annotate()

pass
