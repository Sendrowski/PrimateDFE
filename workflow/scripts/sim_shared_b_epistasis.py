"""
Does sharing the gamma shape b fabricate epistasis?

We simulate a *fixed unscaled* DFE that is deliberately NOT a single gamma
(a discrete step density over fastDFE's default DiscreteFractional bins), so
that no epistasis is present by construction: the unscaled effect distribution
is identical in every population, and the population-scaled DFE is just that
same distribution stretched by 4*Ne. The deterministic consequence is that the
true mean S_d is exactly proportional to Ne (log|S_d|-vs-log Ne slope = 1).

A gamma is closed under scaling (S = 4 Ne s keeps the shape b fixed), so if the
truth were gamma, sharing b across populations would be correctly specified and
recover slope 1. It is not gamma here, so a single shared b cannot match every
population's SFS; the question is whether the resulting misfit is absorbed into
an Ne-correlated scale (= spurious epistasis).

We therefore fit the simulated SFS jointly under two models -- shape b shared
across populations vs. b free per population -- and compare the inferred
log|S_d| vs. log Ne slope to the no-epistasis null of 1.
"""
import logging

import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np
from fastdfe.parametrization import DiscreteFractionalParametrization as DF

from parametrizer import Parametrizer

logging.getLogger("fastdfe").setLevel(logging.ERROR)

# --------------------------------------------------------------------------- #
# Configuration                                                               #
# --------------------------------------------------------------------------- #
n = 8                       # haploid sample size (matches the real data)
n_sites = 1e7               # number of selected sites simulated
mu = 1e-8                   # per-site mutation rate (sets theta = 4 Ne mu)
n_pops = 10
n_runs = 20                 # optimiser restarts per population
n_bootstraps = 100          # bootstraps for the per-population S_d CI
Ne = np.logspace(np.log10(1.3e4), np.log10(3e5), n_pops)   # real-data range
out = "scratch/sim_shared_b_epistasis.png"

# Fixed underlying unscaled DFE: deleterious-only step density over the default
# DiscreteFractional deleterious bins (-1e5,-100), (-100,-10), (-10,-1), (-1,0).
# Nominal masses (sum to 1); converted to the fractional S1..S5 encoding below.
nominal_del = [0.05, 0.20, 0.40, 0.35]
base_edges = np.array([-1e5, -100, -10, -1, 0, 1, 1000])   # default finite edges


def frac_params(a, b, c, d):
    """Fractional S1..S5 encoding a deleterious-only DiscreteFractional DFE
    with nominal bin masses [a,b,c,d] (sum 1) and zero beneficial mass."""
    assert abs(a + b + c + d - 1) < 1e-9
    return dict(S1=a, S2=b / (1 - a), S3=c / (1 - a - b), S4=1.0, S5=0.0)


free = frac_params(*nominal_del)

# Reference Ne at which the unscaled DFE coincides with the default bins; for any
# other Ne the *same* unscaled DFE corresponds to scaling the bin edges, because
# the bins live in S = 4 Ne s space. Geometric centre keeps the scale factors
# symmetric across the range.
Ne_ref = np.exp(np.mean(np.log(Ne)))
pops = [f"Ne_{int(x)}" for x in Ne]

# --------------------------------------------------------------------------- #
# Simulate SFS per population under the fixed unscaled DFE                     #
# --------------------------------------------------------------------------- #
spectra_dict = {}
for pop, x in zip(pops, Ne):
    neut = fd.Simulation.get_neutral_sfs(theta=4 * x * mu, n_sites=n_sites, n=n)
    sim = fd.Simulation(
        sfs_neut=neut,
        model=DF(intervals=base_edges * (x / Ne_ref)),   # same DFE, stretched by Ne
        params=free,
        parallelize=True,
    )
    sim.run()
    sp = sim.get_spectra()
    spectra_dict[f"neutral.{pop}"] = sp["neutral"]
    spectra_dict[f"selected.{pop}"] = sp["selected"]

spectra = fd.Spectra.from_spectra(spectra_dict)


# --------------------------------------------------------------------------- #
# Joint inference: shape shared vs. free                                       #
# --------------------------------------------------------------------------- #
def run_joint(share_b: bool) -> fd.JointInference:
    inf = fd.JointInference(
        sfs_neut=spectra["neutral.*"].merge_groups(level=1),
        sfs_sel=spectra["selected.*"].merge_groups(level=1),
        shared_params=[fd.SharedParams(params=["b"], types="all")] if share_b else [],
        model=fd.GammaExpParametrization(),
        fixed_params={"all": {"h": 0.5, "eps": 0, "p_b": 0, "S_b": 1}},
        n_bootstraps=n_bootstraps,
        n_runs=n_runs,
        parallelize=True,
    )
    inf.run()
    return inf


quantiles = [0.3, 0.5, 0.7]   # deleterious-mass probability levels for S_q


def collect(inf: fd.JointInference):
    """Median |S_d| (+ bootstrap CI), shape b, and median |S_q| per population."""
    Sd, lo, hi, b = [], [], [], []
    Sq = {f: [] for f in quantiles}
    for pop in pops:
        dfe = inf.joint_inferences[pop].get_dfe()
        boot = np.abs(dfe.bootstraps["S_d"].values)
        m = np.median(boot)
        Sd.append(m)
        lo.append(m - np.percentile(boot, 2.5))
        hi.append(np.percentile(boot, 97.5) - m)
        b.append(dfe.params["b"])
        qd = Parametrizer.get_S_quantiles(dfe, quantiles)
        for f in quantiles:
            vals = np.abs(qd[f])
            Sq[f].append(np.median(vals[np.isfinite(vals)]))
    return (np.array(Sd), np.array(lo), np.array(hi), np.array(b),
            {f: np.array(v) for f, v in Sq.items()})


def slope_of(y):
    return np.polyfit(np.log(Ne), np.log(y), 1)


results = {}
for share_b in (True, False):
    inf = run_joint(share_b)
    Sd, lo, hi, b, Sq = collect(inf)
    slope, intercept = slope_of(Sd)
    q_slopes = {f: slope_of(Sq[f])[0] for f in quantiles}
    results[share_b] = dict(Sd=Sd, lo=lo, hi=hi, b=b, Sq=Sq,
                            slope=slope, intercept=intercept, q_slopes=q_slopes)
    tag = "shared b" if share_b else "free b"
    print(f"{tag:10s}: mean |S_d| slope = {slope:.3f}   "
          f"quantile |S_q| slopes "
          + ", ".join(f"f={f:g}:{q_slopes[f]:.3f}" for f in quantiles)
          + f"   b = {np.round(b, 3)}")

# --------------------------------------------------------------------------- #
# Plot                                                                         #
# --------------------------------------------------------------------------- #
fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
styles = {True: ("#C44E52", "shared $b$"), False: ("#4C72B0", "free $b$")}
xx = np.array([Ne.min(), Ne.max()])

for share_b, (color, label) in styles.items():
    r = results[share_b]
    ax.errorbar(Ne, r["Sd"], yerr=[np.maximum(r["lo"], 0), np.maximum(r["hi"], 0)],
                fmt="o", ms=4, color=color, alpha=0.6, capsize=0, lw=0)
    ax.plot(xx, np.exp(r["intercept"]) * xx ** r["slope"], color=color, lw=2,
            label=f"{label} (slope $={r['slope']:.2f}$)")

# slope-1 (no-epistasis) reference through the geometric centre of the data
ymid = np.exp(np.mean([np.log(results[True]["Sd"]).mean(),
                       np.log(results[False]["Sd"]).mean()]))
ax.plot(xx, ymid * (xx / Ne_ref), color="0.55", lw=1.3, ls=(0, (6, 4)), zorder=0,
        label="no epistasis (slope $=1$)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$N_e$")
ax.set_ylabel(r"inferred $|S_d|$")
ax.set_title("Fixed unscaled (non-gamma) DFE: does sharing $b$ fabricate epistasis?")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.grid(True, which="major", ls=":", lw=0.5, color="0.8", alpha=0.7)
ax.set_axisbelow(True)
ax.legend(fontsize=9)
fig.savefig(out, bbox_inches="tight")
print(f"saved {out}")
