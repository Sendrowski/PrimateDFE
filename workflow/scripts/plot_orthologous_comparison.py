"""
Compare the spectra used throughout the study with spectra recomputed on the gene set shared across
all reference annotations.

The comparison addresses whether differences in gene content between annotations could drive the
cross-species patterns: if they could not, restricting every species to a common set of orthologous
genes should leave the spectra, and hence the inference, essentially unchanged.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from populations import Populations

try:
    testing = False
    sfs_files = snakemake.input.sfs
    ortho_files = snakemake.input.orthologous
    order_file = snakemake.input.order
    pops = snakemake.params.pops
    out = snakemake.output[0]
except NameError:
    testing = True


    def get_ref_from_pop(pop: str) -> str:
        s = "_".join(pop.split("_")[:2])
        df = pd.read_csv(f"resources/metadata/{s.split('_')[0]}_individuals.csv", sep="\t")
        return df[df["BAM_FOLDER"].str.contains(s)].iloc[0].REFERENCE_FOLDER.replace("_ssp", "")


    pops = Populations.get_pops(8, "catarrhini")
    sfs_files = [f"results/sfs/{get_ref_from_pop(p)}/{p}/sfs.unfolded.8.csv" for p in pops]
    ortho_files = [f"results/sfs/{get_ref_from_pop(p)}/{p}/sfs.unfolded.orthologous.8.csv" for p in pops]
    order_file = "results/tables/phylogeny/catarrhini/species_order.8.csv"
    out = "scratch/orthologous_comparison.png"


def proportions(counts: np.ndarray) -> np.ndarray:
    """
    Return the spectrum over the polymorphic bins, scaled to sum to one.
    """
    poly = counts[1:-1]
    return poly / poly.sum()


rows = []
for pop, full_path, ortho_path in zip(pops, sfs_files, ortho_files):
    try:
        full = pd.read_csv(full_path)
        ortho = pd.read_csv(ortho_path)
    except FileNotFoundError:
        continue

    for key in ("neutral", "selected"):
        a, b = full[key].to_numpy(), ortho[key].to_numpy()
        rows.append(dict(
            population=pop,
            type=key,
            sites_retained=b[0] / a[0],
            snps_retained=b[1:-1].sum() / a[1:-1].sum(),
            # total variation distance between the normalized spectra
            distance=0.5 * np.abs(proportions(a) - proportions(b)).sum(),
        ))

d = pd.DataFrame(rows)

rank = pd.read_csv(order_file).set_index("species")["rank"]
order = sorted(d.population.unique(),
               key=lambda p: (rank.get("_".join(p.split("_")[:2]), len(rank)), p))
position = {p: i for i, p in enumerate(order)}

fig, ax = plt.subplots(figsize=(6.2, 9.5), dpi=300)

height = 0.38
for offset, (key, color) in zip((height / 2, -height / 2),
                                dict(neutral="#1f77b4", selected="#ff7f0e").items()):
    sub = d[d.type == key]
    ax.barh([position[p] + offset for p in sub.population], sub.distance, height=height,
            color=color, label=key)

ax.set_yticks(range(len(order)))
ax.set_yticklabels([p.replace("_", " ") for p in order], fontsize=8, style="italic")
ax.set_ylim(len(order), -1)
ax.set_xlabel("change in the SFS (total variation distance)", fontsize=10)
ax.tick_params(axis="x", labelsize=9)
ax.legend(frameon=True, fontsize=9, loc="lower right")

fig.tight_layout()
fig.savefig(out, bbox_inches="tight")

d.to_csv(out.replace(".png", ".csv"), index=False)

if testing:
    for key in ("neutral", "selected"):
        sub = d[d.type == key]
        print(f"{key}: {len(sub)} populations, sites retained median {sub.sites_retained.median():.3f}, "
              f"SNPs retained median {sub.snps_retained.median():.3f}, "
              f"distance median {sub.distance.median():.4f} max {sub.distance.max():.4f}")
