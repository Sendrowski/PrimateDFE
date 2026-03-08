import re

import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.container import BarContainer
from matplotlib.patches import Patch

param = dict(S_d=-300, b=0.3, p_b=0, S_b=1)
intervals_del = (-1e8, -1e-5, 1000)
intervals_ben = (1e-5, 1e4, 1000)

n = 8
theta = 1e-4
n_sites = 1e7
scales = [1, 2, 4]
parallelize = True

sim = fd.Simulation(
    sfs_neut=fd.Simulation.get_neutral_sfs(n=n, n_sites=n_sites, theta=theta),
    intervals_ben=intervals_ben,
    intervals_del=intervals_del,
    params=param,
    model=fd.GammaExpParametrization(),
    parallelize=parallelize
)
sim.run()
base_spectra = sim.get_spectra()

fig, (ax_sfs, ax_dfe) = plt.subplots(1, 2, figsize=(8, 3))

dfes = []

spectra = {}
for scale in scales:
    s = base_spectra.copy()
    s.data *= scale

    inf = fd.BaseInference(
        sfs_neut=s['neutral'],
        sfs_sel=s['selected'],
        intervals_ben=intervals_ben,
        intervals_del=intervals_del,
        fixed_params=dict(all=dict(h=0.5, eps=0, p_b=0, S_b=1)),
        n_runs=10,
        n_bootstraps=100,
        parallelize=parallelize,
        discretization=sim.discretization,
    )
    inf.run()

    dfes.append(inf.get_dfe())
    spectra[f"{scale}*neutral"] = s['neutral']
    spectra[f"{scale}*selected"] = s['selected']

fd.Spectra.from_spectra(spectra).plot(
    ax=ax_sfs,
    show=False
)

fd.DFE.plot_many(
    dfes,
    labels=[f"×{s}" for s in scales],
    point_estimate='median',
    ax=ax_dfe,
    show=False,
    intervals=[-np.inf, -100, -10, -1, 0]
)

ax_sfs.set_title("SFS")
ax_sfs.set_xlabel("derived allele count")
ax_sfs.set_ylabel("counts")
ax_sfs.legend(fontsize=8)

ax_dfe.set_title("Inferred DFE")

colors = ["#4C72B0",  # muted blue
          "#55A868",  # muted green
          "#C44E52"]  # muted red

scale2color = {str(s): c for s, c in zip(scales, colors)}
pat = re.compile(r'^(?P<scale>[\d.]+)\*(?P<kind>neutral|selected)$')

for container in ax_sfs.containers:
    lbl = container.get_label()
    m = pat.match(lbl)
    if not m:
        continue  # skip unlabeled containers, if any

    scale = m.group("scale")
    kind = m.group("kind")
    color = scale2color[scale]
    alpha = 0.6 if kind == "neutral" else 1.0

    for bar in container:
        bar.set_facecolor(color)
        bar.set_edgecolor("none")
        bar.set_hatch(None)
        bar.set_alpha(alpha)

bars = [c for c in ax_dfe.containers if isinstance(c, BarContainer)]

for container, color in zip(bars, colors):
    for bar in container:
        bar.set_facecolor(color)
        bar.set_edgecolor("none")
        bar.set_hatch(None)

handles = []
labels = []

for s, c in zip(scales, colors):
    handles.append(Patch(facecolor=c, alpha=0.6, edgecolor="none"))
    labels.append(f"{s}× neutral")
    handles.append(Patch(facecolor=c, alpha=1.0, edgecolor="none"))
    labels.append(f"{s}× selected")

ax_sfs.legend(handles, labels, fontsize=8)

handles = [
    Patch(facecolor=c, edgecolor="none", label=f"×{s}")
    for s, c in zip(scales, colors)
]

ax_dfe.legend(handles=handles, fontsize=8)

plt.tight_layout()
plt.savefig("scratch/dfe_invariance.png", dpi=300)
plt.show()

pass
