"""
Converts a PAF file with CIGAR strings to a fasta file with pairwise alignments.
"""
import re

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

try:
    testing = False
    paf_file = snakemake.input.paf
    ref_fasta = snakemake.input.ref
    query_fasta = snakemake.input.query
    n_chunks = snakemake.params.n_chunks
    i_chunk = snakemake.params.i_chunk
    out_vcf = snakemake.output.get("vcf", None)
    out_fasta = snakemake.output.get("fasta", None)
except NameError:
    # testing
    testing = True
    paf_file = "results/minimap2/Pan_paniscus/Pan_troglodytes.paf"
    ref_fasta = "resources/ref/Pan_paniscus.fasta"
    query_fasta = "resources/ref/Pan_troglodytes.fasta"
    n_chunks = 100
    i_chunk = 99
    out_vcf = None
    out_fasta = "scratch/pairwise_alignment.fasta"


def read_paf(paf_file):
    """Reads a PAF file and returns each line as a dictionary."""
    alignments = []
    with open(paf_file) as f:
        for line in tqdm(f, desc="Loading PAF"):
            fields = line.strip().split('\t')
            alignment = {
                "query_name": fields[0],
                "query_length": int(fields[1]),
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "strand": fields[4],
                "target_name": fields[5],
                "target_length": int(fields[6]),
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "match_bases": int(fields[9]),
                "block_length": int(fields[10]),
                "mapping_quality": int(fields[11]),
            }
            for field in fields[12:]:
                if field.startswith("cg:Z:"):
                    alignment["cigar"] = field[5:]
            alignments.append(alignment)
    return pd.DataFrame(alignments)


def get_chunk_indices(df, n_chunks):
    """Splits DataFrame into chunks with approximately equal total block_length."""
    total = df.block_length.sum()
    target = total / n_chunks
    chunks = []
    current_sum = 0
    start = 0
    for i, length in enumerate(df.block_length):
        current_sum += length
        if current_sum >= target:
            chunks.append((start, i + 1))
            start = i + 1
            current_sum = 0

    if start < len(df):
        chunks.append((start, len(df)))

    # if we only got n_chunks-1, add an empty chunk at the end
    if len(chunks) == n_chunks - 1:
        chunks.append((len(df), len(df)))
    return chunks


def parse_cigar(cigar: str) -> list:
    return [(op, int(length)) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]


alignments = read_paf(paf_file)

# Select only the i-th chunk
chunks = get_chunk_indices(alignments, n_chunks)
start, end = chunks[i_chunk]
alignments = alignments.iloc[start:end].reset_index(drop=True)

ref_fasta = {record.id: record for record in tqdm(SeqIO.parse(ref_fasta, "fasta"), desc="Loading reference FASTA")}
query_fasta = {record.id: record for record in tqdm(SeqIO.parse(query_fasta, "fasta"), desc="Loading query FASTA")}
pbar = tqdm(total=alignments.block_length.sum())


def pairwise_alignment(
        ref_start: int,
        query_start: int,
        query_end: int,
        query_length: int,
        ref_seq: Seq,
        query_seq: Seq,
        cigar: str,
        strand: str,
) -> (Seq, Seq, list):
    """
    Perform pairwise alignment based on the CIGAR string and return aligned sequences and variants.
    """
    if strand == '-':
        query_seq = query_seq.reverse_complement()

    aligned_ref = ''
    aligned_query = ''
    variants = []

    ref_pos = ref_start
    query_pos = query_start
    if strand == '-':
        query_pos = query_length - query_end

    for op, length in parse_cigar(cigar):
        if op == 'M':
            for i in range(length):
                if ref_seq[ref_pos + i] != query_seq[query_pos + i]:
                    variants.append((ref_pos + i + 1, ref_seq[ref_pos + i], query_seq[query_pos + i]))
            aligned_ref += ref_seq[ref_pos:ref_pos + length]
            aligned_query += query_seq[query_pos:query_pos + length]
            ref_pos += length
            query_pos += length
        elif op == '=':
            aligned_ref += ref_seq[ref_pos:ref_pos + length]
            aligned_query += query_seq[query_pos:query_pos + length]
            ref_pos += length
            query_pos += length
        elif op == 'X':
            for i in range(length):
                variants.append((ref_pos + i + 1, ref_seq[ref_pos + i], query_seq[query_pos + i]))
            aligned_ref += ref_seq[ref_pos:ref_pos + length]
            aligned_query += query_seq[query_pos:query_pos + length]
            ref_pos += length
            query_pos += length
        elif op == 'I':
            variants.append((ref_pos + 1, '-', query_seq[query_pos:query_pos + length]))
            aligned_ref += '-' * length
            aligned_query += query_seq[query_pos:query_pos + length]
            query_pos += length
        elif op == 'D':
            variants.append((ref_pos + 1, ref_seq[ref_pos:ref_pos + length], '-'))
            aligned_ref += ref_seq[ref_pos:ref_pos + length]
            aligned_query += '-' * length
            ref_pos += length
        else:
            raise ValueError(f"Unknown CIGAR operation: {op}")

        pbar.update(length)

    return aligned_ref.upper(), aligned_query.upper(), variants


if out_vcf is not None:
    vcf = open(out_vcf, 'w')
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

if out_fasta is not None:
    fasta = open(out_fasta, 'w')

for i, alignment in enumerate(alignments.itertuples()):
    pbar.set_description(f"Processing alignment {i + 1} / {len(alignments)}")

    ref_start = alignment.target_start
    ref_end = alignment.target_end
    query_length = alignment.query_length
    query_start = alignment.query_start
    query_end = alignment.query_end
    ref_id = alignment.target_name
    query_id = alignment.query_name
    cigar_string = alignment.cigar
    strand = alignment.strand

    ref_seq = ref_fasta[ref_id].seq
    query_seq = query_fasta[query_id].seq

    aligned_ref, aligned_query, variants = pairwise_alignment(
        ref_start, query_start, query_end, query_length,
        ref_seq, query_seq, cigar_string, strand
    )

    # make sure the reference coincides with the aligned reference without gaps
    ref_chars = np.array(list(aligned_ref), dtype='U1')
    assert ''.join(ref_chars[ref_chars != '-']) == str(ref_seq[ref_start:ref_end]).upper()

    if out_vcf is not None:
        for pos, ref, alt in variants:
            vcf.write(f"{ref_id}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

    if out_fasta is not None:
        records = [
            SeqRecord(Seq(aligned_ref), id=f"{ref_id}:{ref_start}-{ref_end}", description="Reference"),
            SeqRecord(Seq(aligned_query), id=f"{ref_id}:{ref_start}-{ref_end}",
                      description=f"Query: {query_id}:{query_start}-{query_end}")
        ]
        SeqIO.write(records, fasta, "fasta")

    pbar.update()

pbar.close()

if out_vcf is not None:
    vcf.close()

if out_fasta is not None:
    fasta.close()
