"""
Plot DFE statistic against Ne for multiple populations.
"""

import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

from utils import Parametrizer

try:
    testing = False
    json_files = snakemake.input.json
    ne_files = snakemake.input.ne
    populations = snakemake.params.populations
    stat = snakemake.params.stat
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    json_files = [
        "results/dfe/Homo_sapiens/Homo_sapiens/result.unfolded.10.full.eps.json",
        "results/dfe/Homo_sapiens/Pan_paniscus/result.unfolded.10.full.eps.json",
        "results/dfe/Homo_sapiens/Pongo_pygmaeus/result.unfolded.10.full.eps.json",
        "results/dfe/Homo_sapiens/Pongo_abelii/result.unfolded.10.full.eps.json",
    ]
    ne_files = [
        "results/stats/Ne/Homo_sapiens/Homo_sapiens/Ne.20.csv",
        "results/stats/Ne/Homo_sapiens/Pan_paniscus/Ne.20.csv",
        "results/stats/Ne/Homo_sapiens/Pongo_pygmaeus/Ne.20.csv",
        "results/stats/Ne/Homo_sapiens/Pongo_abelii/Ne.20.csv",
    ]
    populations = [
        "Homo_sapiens",
        "Pan_paniscus",
        "Pongo_pygmaeus",
        "Pongo_abelii",
    ]
    stat = 's_d'
    out = "scratch/dfe_Ne.png"


infs = {}
nes = {}
stats = {}
for pop, json_file, ne_file in zip(populations, json_files, ne_files):
    infs[pop] = fd.BaseInference.from_file(json_file)
    nes[pop] = float(open(ne_file).read())

    if stat == 's_d':
        stats[pop] = abs(Parametrizer.get_S_d(infs[pop].get_dfe()) / nes[pop])

    elif stat == 'S_d':
        stats[pop] = abs(Parametrizer.get_S_d(infs[pop].get_dfe()))

    elif stat == 'range_S_-1_0':
        stats[pop] = infs[pop].get_discretized(np.array([-np.inf, -1, 0, np.inf]))[0][1]

    elif stat == 'range_S_-10_-1':
        stats[pop] = infs[pop].get_discretized(np.array([-np.inf, -10, -1, np.inf]))[0][1]

    elif stat == 'range_S_inf_-10':
        stats[pop] = infs[pop].get_discretized(np.array([-np.inf, -10, np.inf]))[0][0]

    else:
        raise ValueError(f"Unsupported stat: {stat}")

fig, ax = plt.subplots()

# collect the data arrays
data = [stats[pop] for pop in populations]
positions = [nes[pop] for pop in populations]  # x-axis positions = Ne values

# scatter plot of means for each population
means = [np.mean(stats[pop]) for pop in populations]
cmap = plt.get_cmap("tab10")
for i, (pop, x, y) in enumerate(zip(populations, positions, means)):
    ax.scatter(x, y, color=cmap(i % cmap.N), label=pop, s=60)

slope, intercept, r_value, p_value, std_err = linregress(positions, means)
x_vals = np.linspace(min(positions), max(positions), 100)
y_vals = intercept + slope * x_vals
ax.plot(x_vals, y_vals, color="black", linestyle="--", label=f"Regression (r={r_value:.2f}, p={p_value:.3g})")

ax.legend(loc="best", fontsize=7)

ax.set_xlabel("Ne")
ax.set_ylabel(stat)
plt.tight_layout()
plt.savefig(out)

if testing:
    plt.show()
