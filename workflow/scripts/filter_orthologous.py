"""
Filter a VCF down to sites lying inside genes that are annotated in every reference.

This makes the site sets comparable across species, so that spectra and DFE estimates can be
recomputed on a common gene set rather than on each reference's own annotation.
"""

from bisect import bisect_right
from typing import Dict, List, Tuple

import fastdfe as fd
import pandas as pd

try:
    testing = False
    vcf_in = snakemake.input.vcf
    genes_file = snakemake.input.genes
    shared_file = snakemake.input.shared
    vcf_out = snakemake.output[0]
except NameError:
    testing = True
    vcf_in = "results/vcf/Pan_troglodytes/Pan_troglodytes.degeneracy.polarized.vcf.gz"
    genes_file = "results/tables/orthologs/Pan_troglodytes.genes.tsv"
    shared_file = "results/tables/orthologs/catarrhini/shared_groups.8.txt"
    vcf_out = "scratch/Pan_troglodytes.orthologous.vcf.gz"


class GeneSetFiltration(fd.Filtration):
    """
    Retain only sites that fall within one of a given set of gene intervals.

    Intervals are merged per contig and queried by binary search. Sites on contigs the table does not
    cover are dropped, since they carry no gene assignment.
    """

    def __init__(self, genes: 'pd.DataFrame'):
        """
        :param genes: Gene intervals to retain, with columns ``contig``, ``start`` and ``end``.
        """
        super().__init__()

        self.genes = genes
        self._intervals: Dict[str, Tuple[List[int], List[int]]] = {}

    def _load(self):
        """
        Build merged, sorted gene intervals per contig.
        """
        for contig, group in self.genes.groupby("contig"):
            starts, ends = [], []
            for start, end in sorted(zip(group.start, group.end)):
                if starts and start <= ends[-1] + 1:
                    ends[-1] = max(ends[-1], end)
                else:
                    starts.append(start)
                    ends.append(end)

            self._intervals[contig] = (starts, ends)

        self._logger.info("Retaining %d genes across %d contigs", len(self.genes), len(self._intervals))

    def _setup(self, handler):
        """
        Build the interval index before filtering starts.

        :param handler: The handler.
        """
        self._load()

        super()._setup(handler)

    def filter_site(self, variant) -> bool:
        """
        Filter a site by whether it falls inside one of the retained genes.

        :param variant: The variant to filter.
        :return: ``True`` if the site lies within a retained gene, ``False`` otherwise.
        """
        keep = False
        for alias in self._handler.get_aliases(variant.CHROM):
            intervals = self._intervals.get(alias)
            if intervals is None:
                continue

            starts, ends = intervals
            i = bisect_right(starts, variant.POS) - 1
            keep = i >= 0 and variant.POS <= ends[i]
            break

        if not keep:
            self.n_filtered += 1

        return keep


shared = {line.strip() for line in open(shared_file) if line.strip()}
genes = pd.read_csv(genes_file, sep="\t", dtype=dict(group=str))
genes = genes[genes.group.isin(shared)]

f = fd.Filterer(
    vcf=vcf_in,
    output=vcf_out,
    filtrations=[GeneSetFiltration(genes=genes)],
)
f.filter()
