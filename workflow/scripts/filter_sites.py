"""
Remove sites with duplicate alleles or MNPs from a VCF file.
"""
from cyvcf2 import VCF, Writer
from tqdm import tqdm

try:
    testing = False
    vcf_in = snakemake.input[0]
    vcf_out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf_in = "resources/vcf/primates.subset.1000000.vcf.gz"
    vcf_out = "scratch/primates.subset.1000000.no_dups.vcf.gz"

vcf = VCF(vcf_in)
w = Writer(vcf_out, vcf)

filtered = 0
i_sites = 0
for record in tqdm(vcf):
    i_sites += 1

    # recoding didn't work, GATK kept complaining
    if record.REF in record.ALT:
        filtered += 1
        continue

    # we only need SNPs
    if not record.is_snp:
        filtered += 1
        continue

    w.write_record(record)

w.close()
vcf.close()

print(f"Filtered {filtered} / {i_sites} sites.")
