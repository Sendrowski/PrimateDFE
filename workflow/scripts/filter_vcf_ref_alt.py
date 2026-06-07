"""
Filter a VCF file by reference and alternative alleles.
"""

try:
    input = snakemake.input[0]
    ref = snakemake.params.ref
    alt = snakemake.params.alt
    output = snakemake.output[0]
except NameError:
    # testing
    input = "results/vcf/hg38/ingroup/Gorilla_gorilla_gorilla/chr/chr22.vcf.gz"
    ref = "A"
    alt = "T"
    output = f"scratch/variants.mut_{ref}_to_{alt}.vcf"

from cyvcf2 import VCF, Writer
from tqdm import tqdm

# Open the input VCF
vcf_in = VCF(input)

# Create a writer for the output VCF
vcf_out = Writer(output, vcf_in)

# Filter variants based on ref and alt alleles
for variant in tqdm(vcf_in):
    if variant.REF == ref and alt in variant.ALT:
        vcf_out.write_record(variant)

# Close files
vcf_in.close()
vcf_out.close()