"""
Correct a reference's divergence-site list for a focal population by removing
misattributed polymorphism.

The divergence sites (focal reference allele != reconstructed ancestor) come
from a single haploid reference genome, so some are not truly fixed in the
focal population: the reference simply carries the derived allele at a site that
is still segregating. Such sites are detectable because the focal population's
individuals are not all homozygous reference there.

A divergence site is kept iff the focal population is monomorphic for the
reference (= derived) allele at that position, i.e. none of the focal
individuals carries a non-reference allele. Anything segregating (or fixed for a
different allele) in the focal population is dropped. The result is a true
per-population count of lineage-specific fixed substitutions.
"""
import csv
import gzip
import json

import cyvcf2

try:
    sites_file = snakemake.input.sites
    vcf_file = snakemake.input.vcf
    individuals = list(snakemake.params.individuals)
    out_json = snakemake.output.json
except NameError:
    sites_file = "results/divergence/Homo_sapiens/divergence.JC69.sites.tsv.gz"
    vcf_file = "results/vcf/Homo_sapiens/Homo_sapiens.degeneracy.polarized.vcf.gz"
    individuals = None
    out_json = "scratch/divergence/Homo_sapiens.corrected.json"

# GC-conservative substitutions (A<->T, C<->G) are unaffected by gBGC; the
# corrected JSON carries both the all-substitution and GC-conservative-only
# counts so the divergence experiment can be re-run as a gBGC robustness check
# (paired with the filter_bgc polymorphism SFS) for direct comparability.
GC_CONSERVATIVE = {frozenset({"A", "T"}), frozenset({"C", "G"})}

# divergence positions (focal reference allele != reconstructed ancestor)
div = {}     # (chrom, pos) -> degeneracy {0,4}
div_gc = {}  # subset where the ancestral->focal substitution is GC-conservative
with gzip.open(sites_file, "rt") as fh:
    r = csv.reader(fh, delimiter="\t")
    next(r, None)
    for chrom, pos, deg, anc, foc in r:
        key = (chrom, int(pos))
        div[key] = int(deg)
        if frozenset({anc, foc}) in GC_CONSERVATIVE:
            div_gc[key] = int(deg)

raw_neut = sum(1 for d in div.values() if d == 4)
raw_sel = sum(1 for d in div.values() if d == 0)
raw_neut_gc = sum(1 for d in div_gc.values() if d == 4)
raw_sel_gc = sum(1 for d in div_gc.values() if d == 0)

# focal individuals present in the population VCF
vcf = cyvcf2.VCF(vcf_file, samples=(individuals if individuals else None))
n_samples_used = len(vcf.samples)

# divergence positions where the focal population carries any non-reference
# allele (heterozygous or homozygous-alt) -> not fixed for the derived allele.
excluded = set()
for v in vcf:
    key = (v.CHROM, v.POS)
    if key not in div:
        continue
    if v.num_het > 0 or v.num_hom_alt > 0:
        excluded.add(key)

exc_neut = sum(1 for k in excluded if div[k] == 4)
exc_sel = sum(1 for k in excluded if div[k] == 0)
exc_neut_gc = sum(1 for k in excluded if div_gc.get(k) == 4)
exc_sel_gc = sum(1 for k in excluded if div_gc.get(k) == 0)

result = dict(
    sites_file=sites_file,
    vcf_file=vcf_file,
    n_focal_samples=n_samples_used,
    n_div_neut_raw=raw_neut,
    n_div_sel_raw=raw_sel,
    n_excluded_neut=exc_neut,
    n_excluded_sel=exc_sel,
    n_div_neut=raw_neut - exc_neut,
    n_div_sel=raw_sel - exc_sel,
    # GC-conservative-only counts (gBGC robustness check)
    n_div_neut_gc_raw=raw_neut_gc,
    n_div_sel_gc_raw=raw_sel_gc,
    n_excluded_neut_gc=exc_neut_gc,
    n_excluded_sel_gc=exc_sel_gc,
    n_div_neut_gc=raw_neut_gc - exc_neut_gc,
    n_div_sel_gc=raw_sel_gc - exc_sel_gc,
)
with open(out_json, "w") as fh:
    json.dump(result, fh, indent=2)
print(json.dumps(result, indent=2))
