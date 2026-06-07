"""
Scatter a VCF file into N approximately equal chunks using cyvcf2, with tqdm progress bars.
"""

from cyvcf2 import VCF, Writer
from tqdm import tqdm

try:
    input_vcf = snakemake.input[0]
    out_files = snakemake.output
    n_chunks = snakemake.params.n_chunks
except NameError:
    input_vcf = "scratch/test.vcf.gz"
    out_files = [f"scratch/chunk{i}.vcf.gz" for i in range(4)]
    n_chunks = 4

# First pass: count records
print("Counting records...")
vcf = VCF(input_vcf)
total = sum(1 for _ in tqdm(vcf, desc="Counting"))
vcf.close()

# Compute chunk sizes
base = total // n_chunks
extras = total % n_chunks
chunk_limits = [base + (1 if i < extras else 0) for i in range(n_chunks)]

# Second pass: write records
print("Writing chunks...")
vcf = VCF(input_vcf)
writers = [Writer(path, vcf) for path in out_files]

chunk_idx = 0
written = 0
with tqdm(total=total, desc="Writing") as pbar:
    for record in vcf:
        if written == chunk_limits[chunk_idx]:
            writers[chunk_idx].close()
            chunk_idx += 1
            written = 0
        writers[chunk_idx].write_record(record)
        written += 1
        pbar.update(1)

# Close final writer
if chunk_idx < len(writers):
    writers[chunk_idx].close()
vcf.close()