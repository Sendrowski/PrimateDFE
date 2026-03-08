from Bio import SeqIO
import gzip
import sys

def create_identity_liftover(fasta_file, output_chain):
    """
    Generate a liftover chain file from a human FASTA reference genome.
    This maps each chromosome to itself with identical coordinates and compresses the output.
    
    Args:
        fasta_file (str): Path to the input FASTA file.
        output_chain (str): Path to the output chain file (will be .gz compressed).
    """
    temp_file = output_chain.rstrip(".gz")  # Temporary uncompressed file

    try:
        # Write the chain file to a temporary location
        with open(fasta_file, 'r') as fasta, open(temp_file, 'w') as chain:
            for seq_record in SeqIO.parse(fasta, "fasta"):
                chrom = seq_record.id
                length = len(seq_record)
                # Write a chain line mapping the chromosome to itself
                chain.write(f"chain 1000 {chrom} {length} + 0 {length} {chrom} {length} + 0 {length}\n")
                chain.write(f"{length}\n\n")

        # Compress the chain file
        with open(temp_file, 'rb') as temp, gzip.open(output_chain, 'wb') as gzipped_chain:
            gzipped_chain.writelines(temp)

        print(f"Identity liftover chain file created and compressed successfully: {output_chain}")
    except Exception as e:
        print(f"Error while processing FASTA: {e}")
        sys.exit(1)
    finally:
        # Clean up temporary uncompressed file if necessary
        import os
        if os.path.exists(temp_file):
            os.remove(temp_file)

# Example usage
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create an identity liftover chain file from a FASTA reference genome.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to the input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output chain file (e.g., output.chain.gz)")
    args = parser.parse_args()

    create_identity_liftover(args.fasta, args.output)
