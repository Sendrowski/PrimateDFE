"""
Infer the marginal DFE for one population with divergence ON or OFF, using the
GammaExp parametrization and the ``full`` sub-model (free beneficial component,
required for a meaningful alpha). The same divergence-carrying SFS is used for
both modes:

* OFF (``use_divergence=False``): ``include_divergence=False``. fastDFE lumps the
  last bin into the monomorphic count, so the polymorphism likelihood is
  identical to a divergence-free fit (``n_monomorphic = data[0] + data[-1]``).
* ON (``use_divergence=True``): ``include_divergence=True`` with the divergence
  target sizes (``n_sites_div_*`` = each spectrum's ``n_sites``); the last bin
  enters the likelihood as the fixed-substitution count.

Writes the serialized inference (jsonpickle) plus a compact per-population
summary with alpha, its bootstrap CI, and the MK-style ``get_alpha_divergence``.
"""
import json

import fastdfe as fd
import numpy as np

try:
    sfs_file = snakemake.input.sfs
    use_divergence = bool(snakemake.params.use_divergence)
    parametrization = snakemake.params.get("parametrization", "gamma")
    sub_model = snakemake.params.get("sub_model", "full")
    n_bootstraps = snakemake.params.get("n_bootstraps", 100)
    n_runs = snakemake.params.get("n_runs", 100)
    n_bootstrap_retries = snakemake.params.get("n_bootstrap_retries", 10)
    parallelize = snakemake.params.get("parallelize", True)
    out_json = snakemake.output.get("json", None)
    out_summary = snakemake.output.get("summary", None)
    out_dfe = snakemake.output.get("dfe", None)
except NameError:
    sfs_file = "scratch/divergence/Homo_sapiens.sfs.unfolded.div.8.csv"
    use_divergence = True
    parametrization = "gamma"
    sub_model = "full"
    n_bootstraps = 100
    n_runs = 10
    n_bootstrap_retries = 3
    parallelize = True
    out_json = "scratch/divergence/Homo_sapiens.dfe.div_on.json"
    out_summary = "scratch/divergence/Homo_sapiens.dfe.div_on.summary.json"
    out_dfe = "scratch/divergence/Homo_sapiens.dfe.div_on.dfe.json"

spectra = fd.Spectra.from_file(sfs_file)
neut = spectra["neutral"]
sel = spectra["selected"]

if parametrization == "gamma":
    model = fd.GammaExpParametrization(bounds=dict(S_d=(-1e5, -1e-2)))
elif parametrization == "discrete":
    model = fd.DiscreteFractionalParametrization()
else:
    raise ValueError(f"Unknown parametrization {parametrization!r}")

# h fixed to 0.5 (semidominance), eps fixed to 0. The full sub-model leaves the
# beneficial component free (required for a meaningful alpha); the del sub-model
# restructures the DFE to deleterious-only via the parametrization's 'dele'
# submodel (gamma fixes p_b=0,S_b=1; discrete fixes S4=S5=1).
fixed_params = {"all": {"h": 0.5, "eps": 0}}
if sub_model == "del":
    fixed_params["all"].update(model.submodels["dele"])
elif sub_model != "full":
    raise ValueError(f"Unknown sub_model {sub_model!r}")

common = dict(
    model=model,
    fixed_params=fixed_params,
    n_bootstraps=n_bootstraps,
    n_runs=n_runs,
    n_bootstrap_retries=n_bootstrap_retries,
    parallelize=parallelize,
    do_bootstrap=True,
)

if use_divergence:
    inf = fd.BaseInference(
        sfs_neut=neut, sfs_sel=sel,
        include_divergence=True,
        n_sites_div_neut=neut.n_sites,
        n_sites_div_sel=sel.n_sites,
        **common,
    )
else:
    inf = fd.BaseInference(
        sfs_neut=neut, sfs_sel=sel,
        include_divergence=False,
        **common,
    )

inf.run()


def alpha_bootstrap_ci(inference, use_div):
    """alpha at the MLE plus a 95% bootstrap CI and its width."""
    alpha = inference.get_alpha(use_divergence=use_div)
    keys = list(inference.params_mle.keys())
    vals = []
    bs = inference.bootstraps
    if bs is not None:
        for _, row in bs.iterrows():
            p = {k: row[k] for k in keys if k in row}
            try:
                vals.append(inference.get_alpha(params=p, use_divergence=use_div))
            except Exception:
                pass
    vals = np.array([v for v in vals if v == v])
    if len(vals):
        lo, hi = float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
    else:
        lo = hi = float("nan")
    return float(alpha), lo, hi, (hi - lo)


# serialized marginal DFE (carries bootstraps) for the DFE-vs-Ne figure
if out_dfe is not None:
    inf.get_dfe().to_file(out_dfe)

if out_json is not None:
    inf.to_file(out_json)

# the alpha summary is only meaningful for the full sub-model (del fixes p_b=0)
if out_summary is not None:
    alpha, ci_lo, ci_hi, ci_width = alpha_bootstrap_ci(inf, use_divergence)

    # Consistency check (only meaningful when divergence is on): the MK-style
    # alpha estimated from observed divergence -- get_alpha(use_divergence=True)
    # -- vs the DFE-based alpha get_alpha(use_divergence=False), on the same fit.
    alpha_dfe_based = None
    alpha_divergence_based = None
    if use_divergence:
        try:
            alpha_dfe_based = float(inf.get_alpha(use_divergence=False))
            alpha_divergence_based = float(inf.get_alpha(use_divergence=True))
        except Exception:
            pass

    summary = dict(
        sfs_file=sfs_file,
        parametrization=parametrization,
        sub_model=sub_model,
        use_divergence=use_divergence,
        alpha=alpha,
        alpha_ci_lower=ci_lo,
        alpha_ci_upper=ci_hi,
        alpha_ci_width=ci_width,
        alpha_dfe_based=alpha_dfe_based,
        alpha_divergence_based=alpha_divergence_based,
        n_div_neut=float(neut.n_div),
        n_div_sel=float(sel.n_div),
        n_sites_div_neut=float(neut.n_sites),
        n_sites_div_sel=float(sel.n_sites),
        params_mle={k: float(v) for k, v in inf.params_mle.items()},
        log_likelihood=float(inf.likelihood) if hasattr(inf, "likelihood") else None,
    )
    with open(out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
