"""
Merge ancestral allele annotations site-wise for different outgroup configurations.
"""
import ast

import cyvcf2
import numpy as np

try:
    testing = False
    vcf_in = snakemake.input
    vcf_out = snakemake.output
except NameError:
    # testing
    testing = True
    vcf_out = 'scratch/Pan_troglodytes.polarized.vcf.gz'
    vcf_in = [
        'results/vcf/primates.subset.10000.ingroup.Pan_troglodytes.outgroups.PD_0005.PD_0584.polarized.vcf.gz',
        'results/vcf/primates.subset.10000.ingroup.Pan_troglodytes.outgroups.PD_0134.PD_0867.polarized.vcf.gz',
        'results/vcf/primates.subset.10000.ingroup.Pan_troglodytes.outgroups.PD_0431.PD_0071.polarized.vcf.gz',
    ]

readers = [cyvcf2.VCF(vcf) for vcf in vcf_in]

i_annotations = 0

# merge annotations
while True:
    try:
        sites = np.array([reader.__next__() for reader in readers])
    except StopIteration:
        break

    if sites[0].is_snp:
        infos = np.array([site.INFO['AA_info'] for site in sites])

        # whether we have ancestral allele information
        has_ancestral = np.array([info.startswith('{') for info in infos])

        if any(has_ancestral):
            if not has_ancestral.all():
                pass

            has_called_outgroups = np.array(['.' not in ast.literal_eval(info)['outgroup_bases'] for info in infos])

            if has_called_outgroups.any():
                i_annotations += 1

            ancestral_alleles = np.array([s.INFO['AA'] for s in sites[has_ancestral]])

            if not np.all(ancestral_alleles[0] == ancestral_alleles):
                pass



pass
