"""
Add outgroup genotypes to ingroup VCF file
"""

from cyvcf2 import VCF, Writer
from Bio import SeqIO
from tqdm import tqdm

try:
    ingroup_file = snakemake.input.ingroup
    outgroup_files = snakemake.input.outgroups
    out = snakemake.output[0]
except NameError:
    ingroup_file = "input.vcf"
    outgroup_files = ["outgroup1.fasta", "outgroup2.fasta"]
    out = "output.vcf"

# Load outgroup sequences into dict: {sample_name: {chrom: seq}}
outgroup_seqs = {}
for fasta_path in outgroup_files:
    sample = fasta_path.split("/")[-1].replace(".fasta", "")
    outgroup_seqs[sample] = {rec.id: rec.seq.upper() for rec in SeqIO.parse(fasta_path, "fasta")}

# Open VCF
vcf = VCF(ingroup_file)
vcf.add_format_to_header({'ID': 'GT', 'Description': 'Genotype', 'Type': 'String', 'Number': '1'})

# Add outgroup sample names
vcf.set_samples(list(vcf.samples) + list(outgroup_seqs.keys()))
out_vcf = Writer(out, vcf)

# Write records with added outgroup genotypes
for rec in tqdm(vcf, desc="Processing VCF"):
    genotypes = rec.genotypes  # existing ingroup genotypes
    chrom = rec.CHROM
    pos = rec.POS - 1  # 0-based for Python

    ref_allele = rec.REF
    alt_alleles = rec.ALT
    alleles = [ref_allele] + alt_alleles

    # Assign genotype for each outgroup
    for sample in outgroup_seqs:
        base = outgroup_seqs[sample].get(chrom, "N" * (pos + 1))[pos]
        if base == "N" or base not in alleles:
            gt = (-1, -1)  # missing
        else:
            idx = alleles.index(base)
            gt = (idx, idx)
        genotypes.append(gt)

    rec.set_genotypes(genotypes, vcf.samples)
    out_vcf.write_record(rec)

out_vcf.close()