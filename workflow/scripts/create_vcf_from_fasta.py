"""
Create a VCF file from a FASTA file
"""
from typing import Literal

import pysam
from Bio import SeqIO
from tqdm import tqdm
import numpy as np

try:
    testing = False
    ref_name = snakemake.params.ref_name.lower()
    chr = snakemake.params.get("chr", None)  # can be 'other' for all contigs not in chrs
    chrs = snakemake.params.chrs  # list of chromosomes
    mode: Literal["CG", "AT"] = snakemake.params.mode
    fasta_in = snakemake.input.fasta
    max_sites = snakemake.params.get("max_sites", np.inf)
    vcf_out = snakemake.output.vcf
except NameError:
    # testing
    testing = True
    ref_name = "Gorilla_gorilla_gorilla".lower()
    chr = "other"
    chrs = ["chr1", "chr2"]
    mode = "CG"
    fasta_in = f"resources/ref/{ref_name.capitalize()}.fasta"
    max_sites = 100
    vcf_out = f"scratch/{ref_name}.CG.vcf.gz"

# choose alleles based on mode
alleles = tuple(mode)

# Create a VCF file
vcf = pysam.VariantFile(vcf_out, "w")

# GT:DP:GQ:PGT:PID:PL:PS
vcf.header.add_meta('FORMAT', items=[('ID', 'GT'), ('Number', 1), ('Type', 'String'), ('Description', 'Genotype')])
vcf.header.add_meta('FORMAT', items=[('ID', 'DP'), ('Number', 1), ('Type', 'Integer'),
                                     ('Description', 'Depth of reference-supporting reads')])
vcf.header.add_meta('FORMAT',
                    items=[('ID', 'GQ'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Genotype Quality')])
vcf.header.add_meta('FORMAT',
                    items=[('ID', 'PGT'), ('Number', 1), ('Type', 'String'), ('Description', 'Phased Genotype')])
vcf.header.add_meta('FORMAT', items=[('ID', 'PID'), ('Number', 1), ('Type', 'String'), ('Description', 'Phase Set ID')])
vcf.header.add_meta('FORMAT', items=[('ID', 'PL'), ('Number', 'G'), ('Type', 'Integer'),
                                     ('Description', 'Phred-scaled Likelihoods')])
vcf.header.add_meta('FORMAT', items=[('ID', 'PS'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Phase Set')])
vcf.header.info.add("base", number=1, type='String', description="Base information")

vcf.header.add_sample(f"ref_{ref_name}")

# add all contigs to make sure chunks are compatible
[vcf.header.contigs.add(r.id, length=len(r.seq)) for r in tqdm(SeqIO.parse(fasta_in, "fasta"), desc="Adding contigs")]

# lengths = {rec.id: len(rec) for rec in records}
# f_chrs = sum(list(lengths.values())[:24]) / sum(lengths.values())

pbar = tqdm(desc="Creating VCF", total=max_sites)

# iterate over contigs
for record in SeqIO.parse(fasta_in, "fasta"):

    # only consider specified chromosome if specified
    if chr is None or record.id == chr or (chr == "other" and record.id not in chrs):

        # iterate over the sequence
        for i, base in enumerate(record.seq):

            # terminate if more than max_sites
            if pbar.n >= max_sites:
                break

            # Create a VCF record
            vcf_record = vcf.new_record(
                contig=record.id,
                start=i,
                stop=i + 1,
                alleles=alleles,
                qual=100,
                filter="PASS",
                info={},
            )

            vcf_record.info["base"] = str(base).upper()

            vcf_record.samples[f"ref_{ref_name}"]["GT"] = (0,)

            # Add the record to the VCF file
            vcf.write(vcf_record)

            pbar.update()

pbar.close()
vcf.close()
