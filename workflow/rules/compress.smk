"""
Rules for compressing files using bgzip.
"""

# bgzip fasta file
rule bgzip_fasta:
    input:
        '{path}.fasta'
    output:
        '{path}.fasta.gz'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        bgzip -c {input} > {output}
        """