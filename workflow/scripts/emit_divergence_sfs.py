"""
Graft lineage-specific divergence counts onto a population's polymorphism SFS,
producing a divergence-carrying (polyDFE-format) Spectra.

The polymorphism spectrum comes from the usual ``compute_sfs`` output; the last
(divergence) bin is replaced by ``n_div`` from ``count_divergence.py``. Because
divergence is reconstructed from the single haploid reference genome, it is a
property of the species/reference and is shared by all populations of that
species; only the polymorphism part differs between populations.

``Spectrum.from_polydfe`` recomputes the ancestral-monomorphic bin as
``n_sites - polymorphic - n_div``, so the total ``n_sites`` (the divergence
target size used as ``n_sites_div`` at inference time) is preserved.
"""
import json

import fastdfe as fd

try:
    sfs_file = snakemake.input.sfs
    div_file = snakemake.input.divergence
    gc_conservative = snakemake.params.get("gc_conservative", False)
    out = snakemake.output[0]
except NameError:
    sfs_file = "results/sfs/Homo_sapiens/Homo_sapiens/sfs.unfolded.8.csv"
    div_file = "results/divergence/Homo_sapiens/Homo_sapiens/divergence.corrected.json"
    gc_conservative = False
    out = "scratch/divergence/Homo_sapiens.sfs.unfolded.div.8.csv"

spectra = fd.Spectra.from_file(sfs_file)
with open(div_file) as fh:
    d = json.load(fh)

# GC-conservative robustness check: use the GC-conservative-only divergence
# counts grafted onto the GC-conservative (filter_bgc) polymorphism SFS, which
# is selected as ``sfs_file`` upstream. Both restrict to A<->T / C<->G, so the
# resulting alpha is directly comparable to the all-substitution version.
suffix = "_gc" if gc_conservative else ""
n_div_neut = float(d["n_div_neut" + suffix])
n_div_sel = float(d["n_div_sel" + suffix])

neut = spectra["neutral"]
sel = spectra["selected"]

neut_div = fd.Spectrum.from_polydfe(
    list(neut.polymorphic), n_sites=neut.n_sites, n_div=n_div_neut
)
sel_div = fd.Spectrum.from_polydfe(
    list(sel.polymorphic), n_sites=sel.n_sites, n_div=n_div_sel
)

out_spectra = fd.Spectra({"neutral": neut_div, "selected": sel_div})
out_spectra.to_file(out)
print("neutral:", neut_div.to_list())
print("selected:", sel_div.to_list())
