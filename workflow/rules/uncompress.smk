"""
Rules for uncompressing files using bgzip.
"""

# unzip fasta file
rule unzip_fasta:
    input:
        '{path}.fasta.gz'
    output:
        '{path}.fasta'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        gunzip -c {input} > {output}
        """