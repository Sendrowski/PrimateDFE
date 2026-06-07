"""
Recode the genotype of reference allele to match the human reference genome and
adjust the genotype of the original reference accordingly.
We also correct for reverse complementarity between the human reference and the outgroup base
by considering the number of reverse complementarity bases within a window.
Intermediate bases are masked.
"""
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    input = snakemake.input.vcf
    ref = snakemake.input.ref
    window_size = snakemake.params.get('window_size', 6)
    n_complement = snakemake.params.get('n_complement', 5)
    n_mask = snakemake.params.get('n_mask', 2)
    out_vcf = snakemake.output.vcf
    out_stats = snakemake.output.stats
except NameError:
    # Testing
    input = "results/vcf/hg38/outgroup/Pan_troglodytes/AT/chunks/old_ref/0.exons.vcf.gz"
    ref = "resources/ref/hg38.fasta"
    window_size = 6  # window size for reverse complementarity
    n_complement = 5  # number of reverse complementarity bases within the window to take complement
    n_mask = 2  # number of reverse complementarity bases within the window to mask outgroup base
    out_vcf = "scratch/recoded.vcf.gz"
    out_stats = "scratch/recoded.stats.csv"

from cyvcf2 import VCF, Writer
from Bio import SeqIO
from tqdm import tqdm

# Load the reference genome
ref_seq = SeqIO.parse(ref, "fasta")
contig = next(ref_seq)

# Open the input VCF file using cyvcf2
vcf = VCF(input)

# Create a Writer object to write the modified VCF to the output file
w = Writer(out_vcf, vcf)

window = np.ones(window_size, dtype=bool)
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
stats = {'n_divergence': 0, 'n_masked': 0, 'n_complement': 0, 'n_total': 0, 'div_spectrum': defaultdict(int)}

# Process each variant in the VCF
for variant in tqdm(vcf):

    while contig.id != variant.CHROM:
        contig = next(ref_seq)
        window = np.ones(window_size, dtype=bool)

    human_base = contig.seq[variant.POS - 1].upper()
    outgroup_base = variant.INFO.get('base')

    window = np.roll(window, -1)
    window[-1] = (human_base, outgroup_base) in [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]

    # Update REF to match the human reference base
    variant.REF = human_base
    n = window.sum()

    if n < n_mask:
        outgroup_base_comp = outgroup_base
    elif n < n_complement:
        # Mask outgroup base if the window contains below n_mask and n_complement reverse complementarity bases
        outgroup_base_comp = 'N'
        stats['n_masked'] += 1
    else:
        # Take complement of outgroup base if the window contains at least n_complement reverse complementarity bases
        outgroup_base_comp = complement[outgroup_base]
        stats['n_complement'] += 1

    # Update genotype and ALT
    if outgroup_base_comp == human_base:
        variant.genotypes = [[0, 0, False]]
        variant.ALT = ['.']
    else:
        variant.ALT = [outgroup_base_comp]
        variant.genotypes = [[1, 1, False]]

        # Count the number of divergent sites
        if outgroup_base_comp != 'N':
            stats['n_divergence'] += 1
            stats['div_spectrum'][(human_base, outgroup_base_comp)] += 1

    # Write the modified variant to the output VCF
    w.write_record(variant)
    stats['n_total'] += 1

# Close the writer and VCF
w.close()
vcf.close()

stats['div'] = stats['n_divergence'] / stats['n_total'] if stats['n_total'] > 0 else 0
stats['masked'] = stats['n_masked'] / stats['n_total'] if stats['n_total'] > 0 else 0
stats['complement'] = stats['n_complement'] / stats['n_total'] if stats['n_total'] > 0 else 0
stats['div_spectrum'] = str(dict(stats['div_spectrum']))

stats = pd.DataFrame(stats, index=[0])
print(stats.T)

stats.to_csv(out_stats, sep='\t', index=False)