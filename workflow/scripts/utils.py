"""
Project-local utilities, including custom fastDFE filtrations.
"""
from typing import Union

import pysam
from fastdfe.filtration import Filtration, _count_filtered
from fastdfe.io_handlers import DummyVariant


class CpGFiltration(Filtration):
    """
    Drop sites whose REF base is in a CpG dinucleotide context, using a FASTA
    for ±1 lookup. A site is in CpG context iff:

      - REF == 'C' and the next base on the same strand is 'G', OR
      - REF == 'G' and the previous base on the same strand is 'C'.

    The FASTA path is held by the filtration itself rather than by the parent
    handler, since :class:`fastdfe.Filterer` does not currently expose a
    FASTA. Compatible with both ``cyvcf2.Variant`` and
    :class:`fastdfe.io_handlers.DummyVariant` so it works with
    :class:`fastdfe.parser.TargetSiteCounter` as well.
    """

    def __init__(self, fasta: str):
        super().__init__()
        self.fasta_path: str = fasta
        self._fasta: "pysam.FastaFile | None" = None
        self._cur_contig: str | None = None
        self._cur_contig_len: int = 0

    @property
    def fasta(self) -> "pysam.FastaFile":
        # lazy-open: TargetSiteCounter calls filter_site after Parser._teardown,
        # so we keep the handle alive for the lifetime of the object.
        if self._fasta is None:
            self._fasta = pysam.FastaFile(self.fasta_path)
        return self._fasta

    def _contig_length(self, chrom: str) -> int:
        if chrom != self._cur_contig:
            try:
                self._cur_contig_len = self.fasta.get_reference_length(chrom)
            except KeyError:
                self._cur_contig_len = 0
            self._cur_contig = chrom
        return self._cur_contig_len

    def _is_cpg(self, chrom: str, pos: int, ref: str) -> bool:
        if ref not in ('C', 'G'):
            return False
        clen = self._contig_length(chrom)
        if clen == 0:
            return False
        if ref == 'C' and pos < clen:
            try:
                next_b = self.fasta.fetch(chrom, pos, pos + 1).upper()
            except (ValueError, KeyError):
                return False
            return next_b == 'G'
        if ref == 'G' and pos >= 2:
            try:
                prev_b = self.fasta.fetch(chrom, pos - 2, pos - 1).upper()
            except (ValueError, KeyError):
                return False
            return prev_b == 'C'
        return False

    @_count_filtered
    def filter_site(self, variant: Union["cyvcf2.Variant", DummyVariant]) -> bool:
        ref = (variant.REF or '').upper()
        return not self._is_cpg(variant.CHROM, variant.POS, ref)
