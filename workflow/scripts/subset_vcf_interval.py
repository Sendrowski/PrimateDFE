"""
Subset a VCF file to a specific site interval using line index (not genomic coordinates).
"""
from cyvcf2 import VCF, Writer
from tqdm import tqdm

try:
    testing = False
    file_in = snakemake.input[0]
    file_out = snakemake.output[0]
    start = snakemake.params.start
    end = snakemake.params.end
except NameError:
    # testing
    testing = True
    file_in = "results/variants/sample.vcf.gz"
    file_out = "results/variants/sample.subset.100-200.vcf.gz"
    start = 100
    end = 200

vcf = VCF(file_in)
writer = Writer(file_out, vcf)

for i, variant in enumerate(tqdm(vcf)):
    if i < start:
        continue
    if i >= end:
        break
    writer.write_record(variant)

writer.close()
vcf.close()
