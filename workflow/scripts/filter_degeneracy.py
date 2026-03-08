"""
Filter VCF file to retain only 0, and 4-fold degenerate sites.
"""
from cyvcf2 import Reader, Writer
from tqdm import tqdm

try:
    testing = False
    vcf_in = snakemake.input[0]
    vcf_out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf_in = "results/vcf/hgdp.raw.vcf.gz"
    vcf_out = "scratch/hgdp.raw.deg.vcf.gz"

# only retain 0 and 4-fold degenerate sites and remove degeneracy tag
reader = Reader(vcf_in)
writer = Writer(vcf_out, reader)

# Filter and write variants
for variant in tqdm(reader):
    if variant.INFO.get('Degeneracy') in [0, 4]:
        writer.write_record(variant)

# Close writer and reader
writer.close()
reader.close()

pass
