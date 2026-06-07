"""
Filter VCF to drop sites in CpG dinucleotide context.

Uses a fastDFE :class:`Filterer` together with a project-local
:class:`CpGFiltration` (see ``utils.py``).
"""
import fastdfe as fd

from utils import CpGFiltration

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
    filtrations=[CpGFiltration(fasta=fasta_path)],
)

f.filter()
