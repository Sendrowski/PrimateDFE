import re
from collections import defaultdict, namedtuple

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from tqdm import tqdm

try:
    ingroup_file = snakemake.input.ingroup
    outgroup_files = snakemake.input.outgroups
    ref_file = snakemake.input.ref
    sample_names = snakemake.params.sample_names
    out_vcf = snakemake.output.vcf
    out_stats = snakemake.output.stats
except NameError:
    ingroup_file = "results/vcf/ingroup/Macaca_nemestrina/Macaca_tonkeana.exons.biallelic.vcf.gz"
    outgroup_files = [
        "results/minimap2/Macaca_nemestrina/Macaca_nemestrina.fasta",
        "results/minimap2/Macaca_nemestrina/Macaca_mulatta.fasta",
        "results/minimap2/Macaca_nemestrina/Papio_anubis.fasta",
    ]
    ref_file = 'resources/ref/Macaca_nemestrina.fasta'
    sample_names = [
        "Macaca_nemestrina",
        "Macaca_mulatta",
        "Papio_anubis"
    ]
    out_vcf = "scratch/output.vcf.gz"
    out_stats = "scratch/stats.txt"

IntervalSeq = namedtuple("IntervalSeq", ["chrom", "start", "end", "seq"])

query_dict = defaultdict(lambda: defaultdict(dict))
target_dict = defaultdict(lambda: defaultdict(dict))
interval_dict = defaultdict(dict)

ref = {seq.id: seq for seq in tqdm(SeqIO.parse(ref_file, 'fasta'), desc='Reading reference')}

for fasta_file, sample in zip(outgroup_files, sample_names):
    reader = SeqIO.parse(fasta_file, "fasta")
    pbar = tqdm(desc=f"Reading {sample}")
    while True:
        try:
            r = next(reader)
            q = next(reader)
        except StopIteration:
            break
        chrom, start, end = re.match(r'(.+):(\d+)-(\d+)', r.id).groups()
        key = (int(start), int(end))
        target_dict[sample][chrom][key] = r.seq
        query_dict[sample][chrom][key] = q.seq
        interval_dict[sample].setdefault(chrom, []).append(key)
        pbar.update(2)
    pbar.close()

for sample in interval_dict:
    for chrom in interval_dict[sample]:
        interval_dict[sample][chrom] = pd.DataFrame(interval_dict[sample][chrom])

inv = pysam.VariantFile(ingroup_file)
header = inv.header.copy()

if 'GT' not in header.formats:
    header.formats.add("GT", 1, "String", "Genotype")

for name in sample_names:
    header.add_sample(name)

writer = pysam.VariantFile(out_vcf, 'w', header=header)
matches = defaultdict(lambda: defaultdict(int))
current = {}
#curr_ref = {}

for rec in tqdm(inv, desc="Processing VCF"):
    chrom = rec.chrom
    pos = rec.pos - 1  # convert to 0-based index
    new_rec = writer.new_record()

    # Copy basic fields
    new_rec.chrom = rec.chrom
    new_rec.pos = rec.pos
    new_rec.id = rec.id
    new_rec.ref = rec.ref
    new_rec.alts = rec.alts
    new_rec.qual = rec.qual
    if rec.filter.keys():
        new_rec.filter.add(*rec.filter.keys())
    [new_rec.info.__setitem__(
        k,
        v + (0,) if isinstance(v, tuple) and len(v) == len(rec.alts) - 1 else
        (v + (0.0,)) if isinstance(v, tuple) and isinstance(v[0], float) and len(v) == len(rec.alts) - 1 else
        v if isinstance(v, tuple) else (v,)
    ) for k, v in rec.info.items()]

    # Ingroup genotypes
    for s in inv.header.samples:
        new_rec.samples[s]["GT"] = rec.samples[s]["GT"]

    # Outgroup genotypes
    for sample in sample_names:
        if sample not in current or current[sample].chrom != chrom or not (current[sample].start <= pos < current[sample].end):
            if chrom in interval_dict[sample]:
                df = interval_dict[sample][chrom]
                match = df[(df[0] <= pos) & (df[1] > pos)]
                if not match.empty:
                    start, end = match.iloc[0]

                    query_seq = np.array(list(query_dict[sample][chrom][(start, end)]), dtype='U1')
                    target_seq = np.array(list(target_dict[sample][chrom][(start, end)]), dtype='U1')

                    target_gaps = target_seq == '-'

                    assert ''.join(target_seq[~target_gaps]) == ref[chrom].seq[start:end].upper()

                    current[sample] = IntervalSeq(chrom, start, end, ''.join(query_seq[~target_gaps]))
                    #curr_ref[sample] = IntervalSeq(chrom, start, end, ''.join(target_seq[~target_gaps]))
                else:
                    current[sample] = IntervalSeq(chrom, pos, pos, "")
                    #curr_ref[sample] = IntervalSeq(chrom, pos, pos, "")
            else:
                current[sample] = IntervalSeq(chrom, pos, pos, "")
                #curr_ref[sample] = IntervalSeq(chrom, pos, pos, "")

        if not current[sample].seq:
            new_rec.samples[sample]["GT"] = (None, None)
            matches[sample]["no_data"] += 1
            continue

        base = current[sample].seq[pos - current[sample].start].upper()
        #base_ref = curr_ref[sample].seq[pos - curr_ref[sample].start].upper()

        alleles = list(new_rec.alleles)
        if base in ['A', 'C', 'G', 'T'] and base not in alleles:
            new_rec.alts += (base,)
            alleles.append(base)

        gt = (alleles.index(base), alleles.index(base)) if base in alleles else (None, None)
        if base in rec.alleles:
            matches[sample]["match"] += 1
        else:
            matches[sample]["mismatch"] += 1

        new_rec.samples[sample]["GT"] = gt

    writer.write(new_rec)

inv.close()
writer.close()

for sample in matches:
    total = matches[sample]["match"] + matches[sample]["mismatch"]
    matches[sample]["mismatch_fraction"] = matches[sample]["mismatch"] / total if total > 0 else 0

stats = pd.DataFrame(matches).T
print(stats)
stats.to_csv(out_stats)