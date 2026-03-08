"""
Add new INFO fields to a VCF file containing information about each variant.
"""

try:
    input = snakemake.input[0]
    output = snakemake.output[0]
except NameError:
    # testing
    input = "data/variants.vcf"
    output = "scratch/variants_annotated.vcf"

from cyvcf2 import VCF, Writer


def add_info(vcf_input: str, vcf_output: str):
    """
    Add new INFO fields to a VCF file containing information about each variant.
    """
    # Open the input VCF file using cyvcf2
    vcf = VCF(vcf_input)

    # Add new INFO fields to the VCF header for each piece of information
    vcf.add_info_to_header({'ID': 'CHR',
                            'Description': 'Chromosome',
                            'Type': 'String', 'Number': '1'})

    vcf.add_info_to_header({'ID': 'POS',
                            'Description': 'Position',
                            'Type': 'Integer', 'Number': '1'})

    vcf.add_info_to_header({'ID': 'REF',
                            'Description': 'Reference Allele',
                            'Type': 'String', 'Number': '1'})

    vcf.add_info_to_header({'ID': 'ALT',
                            'Description': 'Alternative Allele',
                            'Type': 'String', 'Number': '1'})

    # Create a Writer object to write the modified VCF to the output file
    w = Writer(vcf_output, vcf)

    # Process each variant in the VCF
    for variant in vcf:
        # Extract chromosome, position, reference allele, and alternative alleles
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0]  # We assume only one alternative allele for simplicity

        # Add these fields as separate INFO fields in the VCF
        variant.INFO["CHR"] = chrom
        variant.INFO["POS"] = pos
        variant.INFO["REF"] = ref
        variant.INFO["ALT"] = alt

        # Write the modified variant to the output VCF
        w.write_record(variant)

    # Close the writer
    w.close()
    vcf.close()


add_info(input, output)
