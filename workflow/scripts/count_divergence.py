"""
Count lineage-specific fixed (divergent) substitutions for a focal lineage,
stratified by degeneracy (0-fold = selected, 4-fold = neutral).

The focal lineage is represented by its single haploid reference genome
(``resources/ref/{ref}.fasta``). Outgroup alleles come from the minimap2
alignments onto the focal reference coordinates
(``results/minimap2/{ref}/{outgroup}.fasta``, same paired target/query layout
parsed by ``add_outgroups_vcf.py``), which retain *monomorphic* sites — the
sites that ``filter_biallelic_snps`` strips from the VCFs and that make the SFS
last bin structurally zero.

We restrict to coding sites (CDS in the GFF) only; no ancestral inference is run
on non-coding sites. CpG-context sites (matching ``fastdfe.CpGFiltration``) are
excluded throughout -- from both the divergence counts and the callable-site
target -- so divergence is consistent with the always-CpG-filtered polymorphism
SFS it is grafted onto. For every (non-CpG) coding site where the focal base
differs from at least one outgroup we emit a record, annotate its degeneracy with
fastDFE's own
``DegeneracyAnnotation`` (so the classification matches the SFS pipeline), and
reconstruct the maximum-a-posteriori ancestral allele with Ancestree's
phylogenetic (Felsenstein-pruning) inference over the outgroup ladder. A site is
a lineage-specific fixed substitution when the focal base differs from the
reconstructed ancestor.

Outputs a small JSON with ``n_div_neut`` / ``n_div_sel`` (the divergence counts)
plus diagnostics. ``n_sites_div`` (the divergence target size) is taken from the
SFS target-site counts downstream, consistent with the ``TargetSiteCounter``
used for the SFS.
"""
import json
import os
import re
from collections import defaultdict

import numpy as np
import pysam
from Bio import SeqIO
from tqdm import tqdm

import fastdfe as fd
from ancestree import Inference, JC69, K2, KingmanPolarizationPrior

try:
    testing = False
    ref_file = snakemake.input.ref
    fai_file = snakemake.input.fai
    gff_file = snakemake.input.gff
    outgroup_files = list(snakemake.input.outgroups)
    outgroup_names = list(snakemake.params.outgroup_names)
    focal_name = snakemake.params.focal_name
    model_name = snakemake.params.get("model", "JC69")
    work_vcf = snakemake.output.vcf
    out_json = snakemake.output.json
    out_sites = snakemake.output.sites
    chrom_filter = snakemake.params.get("chrom_filter", None)
except NameError:
    # testing: human focal lineage, single chromosome
    testing = True
    ref_file = "resources/ref/Homo_sapiens.fasta"
    fai_file = "resources/ref/Homo_sapiens.fasta.fai"
    gff_file = "resources/gff/Homo_sapiens.gff"
    outgroup_files = [
        "results/minimap2/Homo_sapiens/Pan_troglodytes.fasta",
        "results/minimap2/Homo_sapiens/Gorilla_gorilla.fasta",
        "results/minimap2/Homo_sapiens/Pongo_abelii.fasta",
    ]
    outgroup_names = ["Pan_troglodytes", "Gorilla_gorilla", "Pongo_abelii"]
    focal_name = "Homo_sapiens"
    model_name = "JC69"
    work_vcf = "scratch/divergence/Homo_sapiens.div.vcf.gz"
    out_json = "scratch/divergence/Homo_sapiens.divergence.json"
    out_sites = "scratch/divergence/Homo_sapiens.divergence.sites.tsv.gz"
    chrom_filter = "chr22"

BASES = ("A", "C", "G", "T")
BASE_SET = set(BASES)

os.makedirs(os.path.dirname(out_json) or ".", exist_ok=True)
os.makedirs(os.path.dirname(work_vcf) or ".", exist_ok=True)


# ---------------------------------------------------------------------------
# 1. CDS intervals per chromosome (0-based, half-open), from the GFF.
# ---------------------------------------------------------------------------
def parse_cds(gff_path, keep_chroms=None):
    starts = defaultdict(list)
    ends = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "CDS":
                continue
            chrom = f[0]
            if keep_chroms is not None and chrom not in keep_chroms:
                continue
            # GFF is 1-based inclusive -> 0-based half-open
            starts[chrom].append(int(f[3]) - 1)
            ends[chrom].append(int(f[4]))
    # Build a sorted unique array of CDS positions per chromosome.
    cds_positions = {}
    for chrom in starts:
        s = np.asarray(starts[chrom])
        e = np.asarray(ends[chrom])
        mask = np.zeros(int(e.max()), dtype=bool)
        for a, b in zip(s, e):
            mask[a:b] = True
        cds_positions[chrom] = np.flatnonzero(mask)  # 0-based positions
    return cds_positions


keep = {chrom_filter} if chrom_filter else None
print(f"Parsing CDS from {gff_file} ...", flush=True)
cds_positions = parse_cds(gff_file, keep_chroms=keep)
print(f"  {sum(len(v) for v in cds_positions.values()):,} coding positions "
      f"over {len(cds_positions)} contigs", flush=True)

# ---------------------------------------------------------------------------
# 2. Focal (reference) base at every coding position.
# ---------------------------------------------------------------------------
print(f"Reading focal reference {ref_file} ...", flush=True)
focal_base = {}  # chrom -> char array aligned to cds_positions[chrom]
cpg_mask = {}    # chrom -> bool array aligned to cds_positions[chrom]; True = CpG context
for rec in SeqIO.parse(ref_file, "fasta"):
    if rec.id not in cds_positions:
        continue
    seq = np.frombuffer(str(rec.seq).upper().encode(), dtype="S1").astype("U1")
    pos = cds_positions[rec.id]
    fb = seq[pos]
    focal_base[rec.id] = fb
    # CpG dinucleotide context, matching fastdfe.CpGFiltration: the reference base
    # is C with the next base G, or G with the previous base C (CpG sites are
    # hypermutable, so excluded -- consistently with the CpG-filtered SFS).
    nlen = seq.shape[0]
    nxt = np.full(pos.shape, "", dtype="U1")
    has_next = pos + 1 < nlen
    nxt[has_next] = seq[pos[has_next] + 1]
    prv = np.full(pos.shape, "", dtype="U1")
    has_prev = pos - 1 >= 0
    prv[has_prev] = seq[pos[has_prev] - 1]
    cpg_mask[rec.id] = ((fb == "C") & (nxt == "G")) | ((fb == "G") & (prv == "C"))
    del seq  # full-chromosome array; not retained (only coding bases kept)
print(f"  focal bases gathered for {len(focal_base)} contigs", flush=True)

# ---------------------------------------------------------------------------
# 3. Outgroup bases at coding positions, from the paired target/query FASTAs.
#    target = reference segment (with gaps), query = outgroup segment.
# ---------------------------------------------------------------------------
n_out = len(outgroup_names)
# og_base[chrom] -> (n_out, n_cds) char array; '.' = no data
og_base = {chrom: np.full((n_out, len(pos)), ".", dtype="U1")
           for chrom, pos in cds_positions.items()}

hdr_re = re.compile(r"(.+):(\d+)-(\d+)")
for oi, (fasta_file, sample) in enumerate(zip(outgroup_files, outgroup_names)):
    print(f"Reading outgroup {sample} <- {fasta_file} ...", flush=True)
    reader = SeqIO.parse(fasta_file, "fasta")
    nblocks = 0
    while True:
        try:
            r = next(reader)  # target (reference) segment
            q = next(reader)  # query (outgroup) segment
        except StopIteration:
            break
        m = hdr_re.match(r.id)
        if not m:
            continue
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        if chrom not in cds_positions:
            continue
        pos = cds_positions[chrom]
        # CDS positions falling inside this block
        lo = np.searchsorted(pos, start, "left")
        hi = np.searchsorted(pos, end, "left")
        if hi <= lo:
            continue
        target = np.frombuffer(str(r.seq).upper().encode(), dtype="S1").astype("U1")
        query = np.frombuffer(str(q.seq).upper().encode(), dtype="S1").astype("U1")
        nongap = target != "-"
        # ref-coordinate -> outgroup base (target gaps removed maps to ref bases)
        aligned = query[nongap]  # length == end - start
        idx = pos[lo:hi] - start  # 0-based offsets within block
        og_base[chrom][oi, lo:hi] = aligned[idx]
        nblocks += 1
    print(f"  {sample}: {nblocks:,} alignment blocks intersect CDS", flush=True)

# ---------------------------------------------------------------------------
# 4. Emit a VCF of variant coding columns (focal differs from >=1 outgroup).
#    Samples: focal (haploid) + each outgroup (haploid). Also tally callable
#    coding sites (focal base ACGT with >=1 outgroup with data) for the
#    Ancestree branch-rate calibration target.
# ---------------------------------------------------------------------------
# contig lengths from .fai for the VCF header
contig_len = {}
with open(fai_file) as fh:
    for line in fh:
        c, ln = line.split("\t")[:2]
        contig_len[c] = int(ln)

vcf_header = pysam.VariantHeader()
for c in cds_positions:
    length = contig_len.get(c) or int(cds_positions[c].max()) + 1
    vcf_header.contigs.add(c, length=length)
for s in [focal_name] + outgroup_names:
    vcf_header.add_sample(s)
vcf_header.formats.add("GT", 1, "String", "Genotype")

writer = pysam.VariantFile(work_vcf, "w", header=vcf_header)

n_considered = 0      # callable coding sites (focal ACGT, >=1 outgroup data, non-CpG)
n_variant = 0         # coding sites emitted (focal != >=1 outgroup)
n_cpg_dropped = 0     # callable coding sites excluded because they are in CpG context
for chrom in cds_positions:
    pos = cds_positions[chrom]
    fb = focal_base[chrom]
    ob = og_base[chrom]                       # (n_out, n_cds)
    focal_ok = np.isin(fb, BASES)
    og_ok = np.isin(ob, BASES)                # (n_out, n_cds)
    any_og = og_ok.any(axis=0)
    non_cpg = ~cpg_mask[chrom]
    callable_all = focal_ok & any_og
    n_cpg_dropped += int((callable_all & ~non_cpg).sum())
    # exclude CpG-context sites from both the callable target and the emitted
    # variants, so n_div is consistent with the CpG-filtered polymorphism SFS
    callable_mask = callable_all & non_cpg
    n_considered += int(callable_mask.sum())
    # variant where an outgroup has a valid base differing from focal
    diff = og_ok & (ob != fb[None, :])
    variant_mask = callable_mask & diff.any(axis=0)
    var_idx = np.flatnonzero(variant_mask)
    for i in var_idx:
        p0 = int(pos[i])
        ref_b = str(fb[i])
        og_bs = [str(ob[k, i]) for k in range(n_out)]
        alts = sorted({b for b in og_bs if b in BASE_SET and b != ref_b})
        alleles = [ref_b] + alts
        rec = writer.new_record()
        rec.chrom = chrom
        rec.start = p0                # 0-based; pysam sets POS = start+1
        rec.alleles = alleles
        rec.samples[focal_name]["GT"] = (0,)
        for k, name in enumerate(outgroup_names):
            b = og_bs[k]
            rec.samples[name]["GT"] = (alleles.index(b),) if b in alleles else (None,)
        writer.write(rec)
        n_variant += 1
writer.close()
print(f"Considered {n_considered:,} callable coding sites; "
      f"emitted {n_variant:,} variant coding columns -> {work_vcf}", flush=True)

pysam.tabix_index(work_vcf, preset="vcf", force=True)

# ---------------------------------------------------------------------------
# 5. Degeneracy annotation (fastDFE) -> (chrom,pos)->degeneracy {0,2,4}.
# ---------------------------------------------------------------------------
deg_vcf = work_vcf.replace(".vcf.gz", ".degeneracy.vcf.gz")
print("Annotating degeneracy with fastDFE ...", flush=True)
fd.Annotator(
    annotations=[fd.DegeneracyAnnotation()],
    vcf=work_vcf,
    gff=gff_file,
    fasta=ref_file,
    output=deg_vcf,
).annotate()

degeneracy = {}
import cyvcf2
for v in cyvcf2.VCF(deg_vcf):
    try:
        degeneracy[(v.CHROM, v.POS)] = int(v.INFO["Degeneracy"])
    except (KeyError, TypeError):
        pass
print(f"  degeneracy annotated for {len(degeneracy):,} sites", flush=True)

# ---------------------------------------------------------------------------
# 6. Ancestral reconstruction (Ancestree) + divergence counting.
# ---------------------------------------------------------------------------
model = JC69() if model_name == "JC69" else K2(fit_kappa=True)
print(f"Reconstructing ancestor with Ancestree ({model_name}) ...", flush=True)
inf = Inference.from_vcf(
    deg_vcf,
    ingroup_samples=[focal_name],
    outgroup_samples=outgroup_names,
    model=model,
    n_target_sites=int(n_considered),
    prior=KingmanPolarizationPrior(ingroup_samples=[focal_name]),
    progress=False,
)
inf.fit()

# GC-conservative substitutions (W<->W = A<->T, S<->S = C<->G) are unaffected by
# gBGC, matching fastdfe.BiasedGCConversionFiltration on the polymorphism side.
# We record the substitution (ancestral -> focal) per site so a downstream
# GC-conservative robustness check can sub-filter the identical fixed-tree JC69
# divergence site set, paired with the filter_bgc (CpG + GC-conservative) SFS.
GC_CONSERVATIVE = {frozenset({"A", "T"}), frozenset({"C", "G"})}


def is_gc_conservative(a, b):
    return frozenset({a, b}) in GC_CONSERVATIVE


n_div_neut = n_div_sel = 0          # 4-fold / 0-fold fixed substitutions
n_div_neut_gc = n_div_sel_gc = 0    # ... restricted to GC-conservative substitutions
n_recon_neut = n_recon_sel = 0     # 4-fold / 0-fold sites reconstructed
# Per-site divergence positions (focal reference allele != reconstructed
# ancestor), so a downstream per-population step can drop sites that are
# actually segregating in the focal population (misattributed polymorphism).
div_sites = []  # (chrom, pos, deg, ancestral, focal)
for site, post in inf.infer():
    deg = degeneracy.get((site.chrom, int(site.pos)))
    if deg not in (0, 4):
        continue
    focal_allele = site.tip_alleles.get(focal_name)
    if focal_allele is None:
        continue
    if deg == 4:
        n_recon_neut += 1
    else:
        n_recon_sel += 1
    if focal_allele != post.map_allele:
        anc = post.map_allele
        div_sites.append((site.chrom, int(site.pos), deg, anc, focal_allele))
        gc = is_gc_conservative(anc, focal_allele)
        if deg == 4:
            n_div_neut += 1
            n_div_neut_gc += int(gc)
        else:
            n_div_sel += 1
            n_div_sel_gc += int(gc)

# write per-site divergence positions (used by correct_divergence.py). Columns:
# chrom, pos, degeneracy, ancestral (reconstructed), focal (reference) allele.
import csv
import gzip
with gzip.open(out_sites, "wt", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["chrom", "pos", "degeneracy", "ancestral", "focal"])
    for chrom, pos, deg, anc, foc in div_sites:
        w.writerow([chrom, pos, deg, anc, foc])

div_mle = (inf.outgroup_divergence_mle.tolist()
           if getattr(inf, "outgroup_divergence_mle", None) is not None else None)

result = dict(
    focal=focal_name,
    model=model_name,
    outgroups=outgroup_names,
    chrom_filter=chrom_filter,
    filter_cpg=True,
    n_considered=int(n_considered),
    n_cpg_dropped=int(n_cpg_dropped),
    cpg_fraction=(n_cpg_dropped / (n_considered + n_cpg_dropped)
                  if (n_considered + n_cpg_dropped) else 0.0),
    n_variant=int(n_variant),
    n_recon_neut=int(n_recon_neut),
    n_recon_sel=int(n_recon_sel),
    n_div_neut=int(n_div_neut),
    n_div_sel=int(n_div_sel),
    n_div_neut_gc=int(n_div_neut_gc),
    n_div_sel_gc=int(n_div_sel_gc),
    outgroup_divergence_mle=div_mle,
)
with open(out_json, "w") as fh:
    json.dump(result, fh, indent=2)
print(json.dumps(result, indent=2), flush=True)
