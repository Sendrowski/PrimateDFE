from collections import defaultdict

import pysam
from Bio import SeqIO
from tqdm import tqdm

try:
    ingroup_file = snakemake.input.ingroup
    ref_file = snakemake.input.ref
    out_stats = snakemake.output.stats
except NameError:
    ingroup_file = "results/vcf/ingroup/Homo_sapiens/Homo_sapiens.exons.biallelic.norm.vcf.gz"
    ref_file = 'resources/ref/Homo_sapiens.fasta'
    out_stats = "scratch/stats.txt"

ref = {seq.id: seq for seq in tqdm(SeqIO.parse(ref_file, 'fasta'), desc='Reading reference')}

inv = pysam.VariantFile(ingroup_file)
matches = defaultdict(int)

for rec in tqdm(inv, desc="Processing VCF"):

    if ref[rec.chrom][rec.pos - 1].upper() == rec.ref.upper():
        matches["match"] += 1
    else:
        matches["mismatch"] += 1
        matches[(ref[rec.chrom][rec.pos - 1].upper(), rec.ref.upper())] += 1

pass
