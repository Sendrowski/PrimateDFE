from cyvcf2 import VCF, Writer
from tqdm import tqdm
import os

try:
    # Snakemake parameters
    input_vcf = snakemake.input[0]
    output_dir = snakemake.output[0]  # Directory where chunk files will be saved
    size = snakemake.params.size  # Number of sites per chunk
except NameError:
    # Testing parameters
    input_vcf = "input.vcf.gz"
    output_dir = "chunks/"  # Directory for output chunks
    size = 1000  # Number of sites per chunk

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Open the input VCF file using cyvcf2
vcf_reader = VCF(input_vcf)

# Initialize chunk index and local record counter
chunk_index = 0
local_count = 0

# Prepare a writer for the first chunk
output_vcf = os.path.join(output_dir, f"{chunk_index}.vcf.gz")
vcf_writer = Writer(output_vcf, vcf_reader)

# Progress bar for tracking
pbar = tqdm(desc="Processing VCF Chunks")

# Iterate over the VCF entries and write them to chunk files
for j, record in enumerate(vcf_reader):
    if local_count < size:
        vcf_writer.write_record(record)
        local_count += 1
    else:
        # Close the current chunk writer and prepare for the next chunk
        vcf_writer.close()

        # Increment chunk index and reset local count
        chunk_index += 1
        local_count = 1  # Start counting for the new chunk

        # Prepare a new output file for the next chunk
        output_vcf = os.path.join(output_dir, f"{chunk_index}.vcf.gz")
        vcf_writer = Writer(output_vcf, vcf_reader)

        # Write the current record to the new chunk
        vcf_writer.write_record(record)

    # Update the progress bar
    pbar.update(1)

# Close the last chunk writer if there's any remaining records
if local_count > 0:
    vcf_writer.close()

pbar.close()
vcf_reader.close()
