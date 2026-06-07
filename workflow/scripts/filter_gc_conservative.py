"""
Filter VCF to retain only GC-conservative mutations (A<->T and C<->G),
which are unaffected by gBGC (biased gene conversion).
"""
from cyvcf2 import Reader, Writer
from tqdm import tqdm

try:
    testing = False
    vcf_in = snakemake.input[0]
    vcf_out = snakemake.output[0]
except NameError:
    testing = True
    vcf_in = "results/vcf/Macaca_nemestrina/Macaca_nigra.degeneracy.polarized.vcf.gz"
    vcf_out = "scratch/Macaca_nigra.degeneracy.polarized.gc_conservative.vcf.gz"

GC_CONSERVATIVE = {frozenset({'A', 'T'}), frozenset({'C', 'G'})}

reader = Reader(vcf_in)
writer = Writer(vcf_out, reader)

n_kept = 0
n_total = 0

for variant in tqdm(reader):
    n_total += 1
    ref = variant.REF.upper()
    alts = variant.ALT

    if len(alts) == 1 and frozenset({ref, alts[0].upper()}) in GC_CONSERVATIVE:
        writer.write_record(variant)
        n_kept += 1

writer.close()
reader.close()

if n_total > 0:
    print(f"Kept {n_kept}/{n_total} GC-conservative variants ({100 * n_kept / n_total:.1f}%)")
