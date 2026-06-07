"""
Concatenate contigs of pairwise alignments to obtain a reference genome sequence.
"""
import gzip
import numpy as np
import re
from Bio import SeqIO
from collections import defaultdict
from tqdm import tqdm

try:
    testing = False
    fasta_in = snakemake.input.fasta
    ref_in = snakemake.input.ref
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    fasta_in = "results/minimap2/hg38/Gorilla_gorilla_gorilla.fasta"
    ref_in = "resources/ref/hg38.fasta"
    out = "scratch/concatenated.fasta.gz"

target = {}
query = {}
reader = SeqIO.parse(fasta_in, "fasta")
pbar = tqdm(reader, desc="Reading sequences")
while True:
    try:
        r = next(reader)
        q = next(reader)
    except StopIteration:
        break

    target[r.id] = r.seq
    query[q.id] = q.seq

    pbar.update(2)

pbar.close()

ref = {record.id: record for record in tqdm(SeqIO.parse(ref_in, "fasta"), desc="Reading reference")}

ids = list(query.keys())

ids = sorted(ids, key=lambda r: int(re.match(r'.+:(\d+)-\d+', r).groups()[0]))

id_dict = defaultdict(list)
for i in ids:
    id_dict[i.split(':')[0]].append(i)

seqs_ref = defaultdict(str)
seqs_query = defaultdict(str)
pbar = tqdm(total=sum(len(r) for r in ref.values()))
gap_length = 0
overlap_length = 0
n_mismatches = 0
n_gaps = 0

# write to file
with gzip.open(out, "wt") as f:
    for chr in ref.keys():

        if chr not in id_dict:
            continue

        pbar.set_description(f"Processing {chr}")
        prev_end = 0

        for i, contig_id in enumerate(id_dict[chr]):
            start, end = map(int, re.match(r'.+:(\d+)-(\d+)', contig_id).groups())

            # skip contigs that end before the previous one
            if prev_end > end:
                end = prev_end  # restore the previous end position
                overlap_length += end - start
                continue

            # handle overlapping contigs
            offset = prev_end - start if prev_end > start else 0
            overlap_length += offset

            # fill in the gap between the previous contig and this one
            seqs_ref[chr] += str(ref[chr].seq[prev_end:start])
            seqs_query[chr] += 'N' * (start - prev_end)
            gap_length += start - prev_end

            target_arr = np.array(list(target[contig_id]), dtype='U1')
            query_arr = np.array(list(query[contig_id]), dtype='U1')
            gaps_target = target_arr == '-'
            gaps_query = query_arr == '-'
            ref_arr = np.array(list(ref[chr].seq[start + offset:end].upper()), dtype='U1')

            # make sure the reference coincides with the target
            assert np.all(target_arr[~gaps_target][offset:] == ref_arr)

            seqs_ref[chr] += "".join(target_arr[~gaps_target][offset:])
            seqs_query[chr] += "".join(query_arr[~gaps_target][offset:])

            pbar.update(end - prev_end)

            gaps = gaps_target | gaps_query
            is_N = (target_arr[~gaps][offset:] == 'N') | (query_arr[~gaps][offset:] == 'N')
            n_mismatches += np.sum(target_arr[~gaps][offset:][~is_N] != query_arr[~gaps][offset:][~is_N])
            n_gaps += np.sum(gaps)

            prev_end = end

            # make sure the reference is correct up to this point
            assert seqs_ref[chr].upper() == ref[chr].seq[:end].upper()

        # add the rest of the reference
        seqs_ref[chr] += str(ref[chr].seq[end:])
        seqs_query[chr] += 'N' * (len(ref[chr]) - end)
        gap_length += len(ref[chr]) - end

        pbar.update(len(ref[chr]) - end)

        # make sure the reference is correct
        assert seqs_ref[chr].upper() == ref[chr].seq.upper()
        assert len(seqs_ref[chr]) == len(seqs_query[chr])

        f.write(f">{chr}\n")
        f.write(seqs_query[chr] + "\n")

    pbar.close()

print(f"Gap length between contigs: {gap_length} ({gap_length / pbar.total * 100:.2f}%)")
print(f"Overlap length: {overlap_length} ({overlap_length / pbar.total * 100:.2f}%)")
print(f"Number of mismatches: {n_mismatches} ({n_mismatches / pbar.total * 100:.2f}%)")
print(f"Number of gaps: {n_gaps} ({n_gaps / pbar.total * 100:.2f}%)")

pass
