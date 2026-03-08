"""
Subset a fasta file to specific number of contigs
"""
from Bio import SeqIO

try:
    testing = False
    file_in = snakemake.input[0]
    n = snakemake.params.n
    file_out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    file_in = "results/variants/chr1.vcf"
    n = 1000
    file_out = "results/variants/chr1_subset.vcf"

records = list(SeqIO.parse(file_in, "fasta"))
SeqIO.write(records[:n], file_out, "fasta")
