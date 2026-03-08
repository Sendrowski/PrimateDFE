"""
Only retain site bi-allelic across the specified samples.
"""
import fastdfe as fd
import numpy as np
import pandas as pd
from fastdfe import SNPFiltration

try:
    testing = False
    file_in = snakemake.input[0]
    samples = snakemake.params.samples
    file_out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    file_in = "results/vcf/hg38/chr/chr1.french_samples_only.norm.french_snps.vcf.gz"
    samples = pd.read_csv('resources/PDGP_french.txt').iloc[:, 0].tolist()
    file_out = "scratch/filtered.vcf.gz"

class PatchedSNPFiltration(SNPFiltration):
    def _prepare_samples_mask(self) -> np.ndarray | None:
        """
        Prepare the samples mask.

        :return: The samples mask.
        """
        self._samples_mask = np.isin(self._handler._reader.samples, self.include_samples)

    def filter_site(self, variant) -> bool:
        """
        Filter site.

        :param variant: The variant to filter.
        :return: ``True`` if the variant is an SNP, ``False`` otherwise.
        """
        return SNPFiltration.filter_site(self, variant)


f = fd.Filterer(
    vcf=file_in,
    output=file_out,
    filtrations=[
        PatchedSNPFiltration(include_samples=samples)
    ]
)

f.filter()
