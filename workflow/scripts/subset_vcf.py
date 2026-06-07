"""
Subset a VCF file to specific number of sites.
"""
import fastdfe as fd

try:
    testing = False
    file_in = snakemake.input[0]
    n = snakemake.params.n
    file_out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    file_in = "results/variants/chr1.vcf"
    n = 1000
    file_out = "results/variants/chr1_subset.vcf"

f = fd.Filterer(
    vcf=file_in,
    output=file_out,
    filtrations=[
        fd.NoFiltration(),
    ],
    max_sites=n
)

f.filter()
