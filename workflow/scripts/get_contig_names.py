"""
Get contig names longer than a threshold number of base pairs.
"""

try:
    testing = False
    fasta = snakemake.input[0]
    threshold = snakemake.params.threshold
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    fasta = "resources/fasta_files/Atele_fusciceps.fasta"
    threshold = 1000000
    out = "scratch/Atele_fusciceps.txt"

from Bio import SeqIO
from tqdm import tqdm

with open(fasta, "r") as fasta_handle, open(out, "w") as out_handle:
    for record in tqdm(SeqIO.parse(fasta_handle, "fasta")):
        if len(record.seq) > threshold:
            out_handle.write(f"{record.id}\n")

    # add 'other' contig which denotes the sum of all contigs shorter than the threshold
    out_handle.write("other\n")