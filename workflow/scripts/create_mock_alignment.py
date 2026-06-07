"""
Creates a mock pairwise alignment FASTA file when no PAF file is provided.
Each sequence in the input FASTA is duplicated as both reference and query.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

try:
    testing = False
    ref_fasta = snakemake.input[0]
    out_fasta = snakemake.output[0]
except NameError:
    # testing
    testing = True
    ref_fasta = "resources/ref/Homo_sapiens.fasta"
    out_fasta = "scratch/mock_alignment.fasta"

ref_seqs = {record.id: record for record in tqdm(SeqIO.parse(ref_fasta, "fasta"), desc="Loading reference")}

fasta = open(out_fasta, 'w')
pbar = tqdm(total=len(ref_seqs), desc="Writing mock alignments")

for record in ref_seqs.values():
    aligned_ref = str(record.seq).upper()
    aligned_query = aligned_ref  # identical
    region = f"{record.id}:0-{len(record)}"
    records = [
        SeqRecord(Seq(aligned_ref), id=region, description="Reference"),
        SeqRecord(Seq(aligned_query), id=region, description=f"Query: {record.id}:0-{len(record)}")
    ]
    SeqIO.write(records, fasta, "fasta")
    pbar.update()

pbar.close()
fasta.close()