"""
Lifts over a VCF file using a PAF file with CIGAR strings.
Outputs a new VCF where positions are mapped to the ref genome.
"""
import functools
import re
from collections import defaultdict
from typing import List

import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

try:
    testing = False
    paf_file = snakemake.input.paf
    vcf_in = snakemake.input.vcf
    ref_fasta = snakemake.input.ref
    query_fasta = snakemake.input.query
    vcf_out = snakemake.output.vcf
    stats_out = snakemake.output.stats
except NameError:
    testing = True
    paf_file = "results/minimap2/Pan_troglodytes/Gorilla_gorilla_gorilla.paf"
    vcf_in = "results/vcf/ingroup/Gorilla_gorilla_gorilla/Gorilla_gorilla_gorilla.exons.vcf.gz"
    ref_fasta = "resources/ref/Pan_troglodytes.fasta"
    query_fasta = "resources/ref/Gorilla_gorilla_gorilla.fasta"
    vcf_out = "scratch/lifted.vcf.gz"
    stats_out = "scratch/stats.csv"


class LazyFasta:
    def __init__(self, fasta_path):
        self.it = SeqIO.parse(fasta_path, "fasta").__iter__()
        self._cache = {}

    def __getitem__(self, chrom: str) -> Seq:
        while chrom not in self._cache:
            seq = self.it.__next__()
            self._cache[seq.id] = seq

        return self._cache[chrom]

    def load_all(self):
        for seq in tqdm(self.it, desc='Loading contigs'):
            self._cache[seq.id] = seq

    @functools.cached_property
    def contigs(self) -> List[str]:
        self.load_all()
        return list(self._cache.keys())


def parse_cigar(cigar):
    return [(op, int(length)) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def liftover_pos(pos, strand, query_start, ref_start, ref_end, cigar, query_len):
    ref_pos = ref_start
    query_pos = query_start if strand == '+' else query_len - ref_end

    for op, length in parse_cigar(cigar):
        if op in ('M', '=', 'X'):
            if query_pos <= pos < query_pos + length:
                return ref_pos + (pos - query_pos)
            ref_pos += length
            query_pos += length
        elif op == 'I':
            query_pos += length
        elif op == 'D':
            ref_pos += length
        else:
            raise ValueError(f"Unsupported CIGAR op: {op}")

    return None


# Load FASTA sequences
ref_seqs = LazyFasta(ref_fasta)
query_seqs = LazyFasta(query_fasta)

# Load and group PAF by query_name
paf_records = []
with open(paf_file) as f:
    for line in tqdm(f, desc="Loading PAF"):
        fields = line.strip().split('\t')
        cigar = next((x[5:] for x in fields[12:] if x.startswith("cg:Z:")), None)
        if not cigar:
            continue
        paf_records.append({
            "query_name": fields[0],
            "query_length": int(fields[1]),
            "query_start": int(fields[2]),
            "query_end": int(fields[3]),
            "strand": fields[4],
            "target_name": fields[5],
            "target_length": int(fields[6]),
            "target_start": int(fields[7]),
            "target_end": int(fields[8]),
            "cigar": cigar
        })

paf_df = pd.DataFrame(paf_records)
paf_by_chrom = {k: g.reset_index(drop=True) for k, g in paf_df.groupby("query_name")}

# Open input/output VCF
reader = pysam.VariantFile(vcf_in, 'r')
for seq_id in ref_seqs.contigs:
    if seq_id not in reader.header.contigs:
        reader.header.contigs.add(seq_id, length=len(ref_seqs[seq_id]))
writer = pysam.VariantFile(vcf_out, 'w', header=reader.header)

stats = defaultdict(int)
prev_hit = None

# Process each variant
for record in tqdm(reader, desc="Processing VCF"):
    chrom = record.chrom
    pos = record.pos - 1  # Convert to 0-based index

    if len(record.ref) != 1:
        stats['skipped_indel'] += 1
        continue

    if (
            prev_hit is not None and
            prev_hit.query_name == chrom and
            prev_hit.query_start <= pos < prev_hit.query_end
    ):
        hit = prev_hit
    else:
        if chrom not in paf_by_chrom:
            stats['missing_chrom'] += 1
            continue

        chrom_hits = paf_by_chrom[chrom]
        hits = chrom_hits[(chrom_hits["query_start"] <= pos) & (chrom_hits["query_end"] > pos)]
        if hits.empty:
            stats['missing_hit'] += 1
            continue

        hit = hits.iloc[0]
        prev_hit = hit

    lifted_pos = liftover_pos(
        pos=pos,
        strand=hit.strand,
        query_start=hit.query_start,
        ref_start=hit.target_start,
        ref_end=hit.target_end,
        query_len=hit.query_length,
        cigar=hit.cigar
    )

    if lifted_pos is None:
        stats['invalid_lift'] += 1
        continue

    query_base = query_seqs[hit.query_name][pos]
    ref_base = ref_seqs[hit.target_name][lifted_pos]

    if record.ref.upper() != query_base.upper():
        print(
            f"Mismatch at {hit.query_name}:{pos + 1}, ",
            f"ref: {record.ref.upper()}, query: {query_base.upper()}, ",
            f"alleles: {record.alleles}, strand: {hit.strand}, ",
            f"lifted_pos: {lifted_pos}, target_name: {hit.target_name}, ",
            f"query_length: {hit.query_length}"
        )

    if record.ref.upper() != query_base.upper():
        stats['ref_mismatch'] += 1

    query_base_next = query_seqs[hit.query_name][pos + 1] if pos + 1 < hit.query_length else 'N'
    ref_base_next = ref_seqs[hit.target_name][lifted_pos + 1] if lifted_pos + 1 < hit.target_length else 'N'

    if hit.strand == '-':
        query_base = reverse_complement(query_base)
        query_base_next = reverse_complement(query_base_next)

    if query_base.upper() == ref_base.upper():
        stats['matches'] += 1
    else:
        stats['mismatches'] += 1

    if query_base_next.upper() != ref_base_next.upper():
        stats['mismatches_next_base'] += 1

    record.chrom = hit.target_name
    record.pos = lifted_pos + 1  # Convert back to 1-based index
    writer.write(record)

reader.close()
writer.close()

# Output stats
n_lifted = stats['matches'] + stats['mismatches']
n_unlifted = stats['missing_hit'] + stats['invalid_lift'] + stats['missing_chrom']
stats['fraction_mismatches'] = stats['mismatches'] / n_lifted if n_lifted > 0 else 0
stats['fraction_mismatches_next_base'] = stats['mismatches_next_base'] / n_lifted if n_lifted > 0 else 0
stats['n_sites'] = n_lifted + n_unlifted
stats['fraction_lifted'] = n_lifted / stats['n_sites']

df = pd.DataFrame(stats, index=[0])
print(df.to_string())
df.to_csv(stats_out, index=False)

pass
