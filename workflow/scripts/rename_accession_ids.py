"""
Replace RefSeq contig names with accession ids using in the VCF file.
"""

import pandas as pd

try:
    testing = False
    gff_file = snakemake.input.gff
    report_file = snakemake.input.report
    output_file = snakemake.output[0]
except NameError:
    testing = True
    gff_file = "scratch/annotations.gff"
    report_file = "scratch/assembly_report.tsv"
    output_file = "scratch/annotations.renamed.gff"

# Load mapping
report = pd.read_csv(report_file, sep='\t')
mapping = report[['Source name', 'Target name']].dropna()
mapping = mapping.set_index('Source name')['Target name'].to_dict()

# Load GFF (skip comment lines)
gff = pd.read_csv(gff_file, sep='\t', header=None, comment='#')

# Replace contig names
gff[0] = gff[0].map(mapping).fillna(gff[0])

# Write modified GFF (no headers preserved)
gff.to_csv(output_file, sep='\t', header=False, index=False)