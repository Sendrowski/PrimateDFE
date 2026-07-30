"""
Filter VCF to drop sites in CpG dinucleotide context.

Uses a fastDFE :class:`Filterer` with fastDFE's :class:`CpGFiltration`
(available since v1.3.0; it pulls the ``±1`` FASTA context from the Filterer).
"""
import fastdfe as fd

try:
    vcf_in = snakemake.input.vcf
    fasta_path = snakemake.input.fasta
    vcf_out = snakemake.output[0]
except NameError:
    # testing
    vcf_in = "results/vcf/Pan_troglodytes/Pan_troglodytes.degeneracy.existing_outgroups.vcf.gz"
    fasta_path = "resources/ref/Pan_troglodytes.fasta"
    vcf_out = "scratch/Pan_troglodytes.degeneracy.existing_outgroups.no_cpg.vcf.gz"

f = fd.Filterer(
    vcf=vcf_in,
    output=vcf_out,
    fasta=fasta_path,
    filtrations=[fd.CpGFiltration()],
)

f.filter()
