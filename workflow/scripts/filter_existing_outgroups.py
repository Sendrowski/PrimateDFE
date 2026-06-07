"""
Filter VCF for existing outgroup genotypes.
"""

import fastdfe as fd
import pandas as pd
from cyvcf2 import Variant
from fastdfe.filtration import _count_filtered
from fastdfe.io_handlers import DummyVariant, get_called_bases

try:
    testing = False
    vcf_in = snakemake.input[0]
    outgroups = snakemake.params.outgroups
    vcf_out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    vcf_in = "results/vcf/Plecturocebus_cupreus/Cheracebus_lucifer.degeneracy.existing_outgroups.vcf.gz"
    outgroups = ['Pithecia_pithecia', 'Sapajus_apella']
    vcf_out = 'scratch/Pan_troglodytes.polarized.vcf.gz'

# initialize filterer
f = fd.Filterer(
    vcf=vcf_in,
    output=vcf_out,
    filtrations=[
        fd.ExistingOutgroupFiltration(
            outgroups=outgroups,
            n_missing=2
        )
    ]
)

f.filter()

pass
