"""
Plot Ne estimated from Watterson's theta against Ne estimated from nucleotide diversity pi.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import LogLocator, FuncFormatter

from populations import Populations
from visualization import DFEvsNePlotter

try:
    testing = False
    sfs_files = snakemake.input.sfs
    species_file = snakemake.input.species
    pops = snakemake.params.pops
    mus = snakemake.params.mus
    out = snakemake.output[0]
except NameError:
    testing = True


    def get_ref_from_pop(pop: str) -> str:
        s = "_".join(pop.split("_")[:2])
        df = pd.read_csv(f"resources/metadata/{s.split('_')[0]}_individuals.csv", sep="\t")
        return df[df["BAM_FOLDER"].str.contains(s)].iloc[0].REFERENCE_FOLDER.replace("_ssp", "")


    pops = Populations.get_pops(8, "catarrhini")
    refs = [get_ref_from_pop(p) for p in pops]
    sfs_files = [f"results/sfs/{r}/{p}/sfs.unpolarized.folded.8.csv" for r, p in zip(refs, pops)]
    species_file = "resources/Kuderna/species.csv"
    _k = pd.read_csv(species_file).set_index("SPECIES_BINOMIAL").MU_PER_GENERATION
    mus = [_k.get("_".join(p.split("_")[:2]), _k[r]) for p, r in zip(pops, refs)]
    out = "scratch/Ne_theta_vs_pi.png"


def theta_watterson(counts: np.ndarray, n: int) -> float:
    """
    Watterson's estimator per site, matching sfsutils.Spectrum.theta.
    """
    return counts[1:-1].sum() / np.sum(1 / np.arange(1, n)) / counts.sum()


def pi(counts: np.ndarray, n: int) -> float:
    """
    Nucleotide diversity per site. The pairwise weights are symmetric under i <-> n - i,
    so the same expression applies to folded and unfolded spectra.
    """
    i = np.arange(len(counts))
    return float((2 * i * (n - i) / (n * (n - 1)) * counts).sum() / counts.sum())


# -------------------- load --------------------

rows = []
for pop, mu, sfs_file in zip(pops, mus, sfs_files):
    counts = pd.read_csv(sfs_file)["neutral"].to_numpy()
    n = len(counts) - 1

    rows.append(dict(
        population=pop,
        group=Populations.get_group_from_pop(pop),
        Ne_watterson=theta_watterson(counts, n) / (4 * mu),
        Ne_pi=pi(counts, n) / (4 * mu),
    ))

df = pd.DataFrame(rows)

label_to_color = Populations.get_label_color_map(df.group.tolist())

# -------------------- plot --------------------

plt.figure(figsize=(5, 5), dpi=400)

plt.scatter(
    df.Ne_watterson,
    df.Ne_pi,
    s=70,
    alpha=0.8,
    color=[label_to_color[g] for g in df.group],
)

ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")

lim = [0.9 * min(df.Ne_watterson.min(), df.Ne_pi.min()), 1.1 * max(df.Ne_watterson.max(), df.Ne_pi.max())]
xx = np.logspace(np.log10(lim[0]), np.log10(lim[1]), 200)
plt.plot(xx, xx, ":", color="black", linewidth=1, label="$y=x$")

ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_aspect("equal", adjustable="box")

handles = [
    plt.Line2D([0], [0], marker="o", linestyle="", color=color, label=label)
    for label, color in label_to_color.items()
]
r = np.corrcoef(np.log10(df.Ne_watterson), np.log10(df.Ne_pi))[0, 1]

handles.append(plt.Line2D([0], [0], linestyle=":", color="black", label=f"$y=x$ ($r={r:.3f}$)"))

plt.legend(handles, [Populations.label_to_text(h.get_label()) for h in handles], fontsize=9)

plt.xlabel("$N_e$ (Watterson)")
plt.ylabel(r"$N_e$ ($\pi$)")
plt.title(r"$N_e$ (Watterson) vs $N_e$ ($\pi$)")

for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    axis.set_major_formatter(FuncFormatter(DFEvsNePlotter.log_label_pow))

plt.tight_layout()
plt.savefig(out, bbox_inches="tight")

if testing:
    plt.show()
