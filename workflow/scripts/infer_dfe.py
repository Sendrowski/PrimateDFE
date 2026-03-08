"""
Infer the DFE.
"""
import re

import fastdfe as fd
import numpy as np
from matplotlib import pyplot as plt

try:
    testing = False
    sfs_file = snakemake.input[0]
    parametrization = snakemake.params.parametrization
    sub_model = snakemake.params.sub_model
    use_eps = snakemake.params.use_eps
    n_bootstraps = snakemake.params.get("n_bootstraps", 100)
    n_runs = snakemake.params.get("n_runs", 100)
    n_bootstrap_retries = snakemake.params.get("n_bootstrap_retries", 10)
    parallelize = snakemake.params.get("parallelize", False)
    recessive = snakemake.params.get("recessive", False)
    out_discretized = snakemake.output.discretized
    out_sfs_comparison = snakemake.output.sfs_comparison
    out_json = snakemake.output.json
    out_dfe = snakemake.output.dfe
except NameError:
    # testing
    testing = True
    sfs_file = "results/sfs/Homo_sapiens/Homo_sapiens/sfs.unfolded.filter_bgc.8.csv"
    parametrization = "gamma"
    sub_model = "fixed_h=0.5"
    use_eps = False
    n_bootstraps = 100
    n_runs = 10
    n_bootstrap_retries = 2
    parallelize = True
    recessive = False
    out_discretized = "scratch/discretized.png"
    out_sfs_comparison = "scratch/sfs_comparison.png"
    out_json = "scratch/results.json"
    out_dfe = "scratch/dfe.json"

spectra = fd.Spectra.from_file(sfs_file)
fixed_params = {'all': {'h': 0.5}}

# results not affected by number of target sites
# spectra.data.iloc[0] *= 10

# choose model
if parametrization == "gamma":
    model = fd.GammaExpParametrization(bounds=dict(S_d=(-1e5, -1e-2)))
elif parametrization == "discrete":
    model = fd.DiscreteFractionalParametrization()
else:
    raise ValueError(f"No parametrization {parametrization}")

# defaults for deleterious-only models
if not 'full' in sub_model:
    if parametrization == "gamma":
        fixed_params["all"] |= {"p_b": 0, "S_b": 1}
    else:  # discrete
        fixed_params["all"] |= {"S4": 1, "S5": 0}

# fixed parameters: fixed_param=value_param=value
elif sub_model.startswith("fixed"):
    for k, v in re.findall(r"_(\w+)=([A-Za-z0-9+-.]+)", sub_model.replace("fixed", "")):
        if v == "None":
            fixed_params["all"].pop(k, None)
        else:
            fixed_params["all"][k] = float(v)

if not use_eps:
    fixed_params['all']['eps'] = 0

inf = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    # intervals_del=(-1.0e+10, -1.0e-5, 1000),
    n_bootstraps=n_bootstraps,
    n_runs=n_runs,
    n_bootstrap_retries=n_bootstrap_retries,
    parallelize=parallelize,
    model=model,
    fixed_params=fixed_params,
    h_callback=(lambda _, S: 0.4 * np.exp(-0.02 * abs(S))) if recessive else None
)

inf.run()

params = list(list(inf.params_mle.keys()) - fixed_params['all'].keys())
title = ', '.join([f"{k}={v}" for k, v in inf.bootstraps[params].mean().round(2).items()])
inf.plot_discretized(show=testing, title=title)
plt.savefig(out_discretized)

if out_dfe is not None:
    inf.get_dfe().to_file(out_dfe)

inf.plot_sfs_comparison(show=testing)
plt.xlabel('derived allele frequency')
plt.ylabel('number of variants')
plt.tight_layout()
plt.savefig(out_sfs_comparison)

# save inference
inf.to_file(out_json)

pass
