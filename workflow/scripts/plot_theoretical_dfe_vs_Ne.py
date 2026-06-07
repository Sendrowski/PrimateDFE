"""
Plot theoretical DFE bin fractions vs. Ne for a fixed unscaled GammaExp DFE.

Unscaled DFE parameters are derived from inferred population-level estimates:
  - b and p_b are averaged directly across populations
  - S_d and S_b are first unscaled (divided by 4*Ne) then averaged
"""
import ast

import numpy as np
import pandas as pd
import fastdfe as fd
from fastdfe.parametrization import GammaExpParametrization

from visualization import DFEvsNePlotter

# Discretization used to compute alpha (proportion of beneficial substitutions)
_model = GammaExpParametrization()
_disc = fd.Discretization(n=20)

# --------------------------------------------------------------------------- #
# Derive unscaled DFE parameters from inferred per-population estimates        #
# --------------------------------------------------------------------------- #
dfe_file = "results/tables/dfe/catarrhini/dfe.unfolded.8.gamma.full.noeps.csv"
ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"

dfe_df = pd.read_csv(dfe_file)
dfe_df['params'] = dfe_df['params'].apply(ast.literal_eval)

ne_df = pd.read_csv(ne_file)
ne_map = dict(zip(ne_df['label'], ne_df['x']))

rows = []
for _, row in dfe_df.iterrows():
    Ne = ne_map.get(row['population'])
    if Ne is None:
        continue
    p = row['params']
    rows.append(dict(
        b=p['b'],
        p_b=p['p_b'],
        s_d=p['S_d'] / (4 * Ne),
        s_b=p['S_b'] / (4 * Ne),
    ))

scaled = pd.DataFrame(rows)
unscaled = scaled.median()

s_d = -0.1
b = 0.18
p_b = 0.2
s_b = 1e-5

print(f"Unscaled DFE parameters (mean across {len(rows)} populations):")
print(f"  s_d = {s_d:.4g}")
print(f"  b   = {b:.4g}")
print(f"  p_b = {p_b:.4g}")
print(f"  s_b = {s_b:.4g}")

# --------------------------------------------------------------------------- #
# Grid of Ne values                                                            #
# --------------------------------------------------------------------------- #
Ne_values = np.logspace(np.log10(1.3e4), np.log10(3e5), 40)
populations = [f"Ne_{int(ne)}" for ne in Ne_values]

ne_dict = dict(zip(populations, Ne_values))

# --------------------------------------------------------------------------- #
# Build a DFE object per Ne (scale s → S = 4 Ne s)                           #
# --------------------------------------------------------------------------- #
dfes = {}
for pop, Ne in zip(populations, Ne_values):
    params = dict(
        S_d=s_d * 4 * Ne,
        b=b,
        p_b=p_b,
        S_b=s_b * 4 * Ne,
    )
    # Compute alpha (proportion of beneficial substitutions) from the scaled params
    alpha = _disc.get_alpha(_model, params)

    # A single-row bootstrap DataFrame acts as the point estimate;
    # all percentile intervals collapse to zero (no sampling uncertainty).
    bootstraps = pd.DataFrame([{**params, 'alpha': alpha}])
    dfes[pop] = fd.DFE(params=params, bootstraps=bootstraps)

# --------------------------------------------------------------------------- #
# Statistics (population-scaled S bins)                                        #
# --------------------------------------------------------------------------- #
stat_list = [
    "range_S_-inf_-10",
    "range_S_-10_-1",
    "range_S_-1_0",
    "range_S_0_inf",
    "alpha",
]

labels = ["theoretical"] * len(populations)
label_color = {"theoretical": "steelblue"}

# --------------------------------------------------------------------------- #
# Plot                                                                         #
# --------------------------------------------------------------------------- #
plotter = DFEvsNePlotter(
    dfes=dfes,
    ne_dict=ne_dict,
    populations=populations,
    stat_list=stat_list,
    labels=labels,
    label_color=label_color,
    reg_type="linear",
)

fig = plotter.plot(
    file=None,
    show=False,
    show_legend=False,
    style="stacked",
)

# The stacked panels are drawn with zero vertical spacing, so matplotlib's
# default y-ticks are too dense and collide across panel boundaries. Cap the
# number of ticks per panel and prune the extreme ones (figure-local; does not
# affect the shared DFEvsNePlotter used by the data figures).
from matplotlib.ticker import MaxNLocator

for ax in fig.axes[:len(stat_list)]:
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))

fig.savefig("scratch/theoretical_dfe_vs_Ne.png", dpi=400)
