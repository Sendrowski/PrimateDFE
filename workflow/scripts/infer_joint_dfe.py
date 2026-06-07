"""
Infer the joint DFE from the SFS using fastDFE.
"""

__author__ = "Janek Sendrowski"
__contact__ = "sendrowski.janek@gmail.com"
__date__ = "2023-02-26"

import re

import fastdfe as fd

try:
    testing = False
    spectra_file = snakemake.input[0]
    parametrization = snakemake.params.parametrization
    sub_model = snakemake.params.sub_model
    use_eps = snakemake.params.use_eps
    parallelize = snakemake.params.get('parallelize', False)
    shared_params = snakemake.params.shared_params
    out_json = snakemake.output.json
    out_discretized = snakemake.output.discretized
    out_sfs_comparison = snakemake.output.sfs_comparison
except NameError:
    # testing
    testing = True
    spectra_file = "results/sfs/comp/Homo_sapiens/Castellano_data/sfs.unfolded.8.csv"
    parametrization = "gamma"
    sub_model = 'del'
    use_eps = False
    parallelize = True
    shared_params = ['b']
    out_json = "scratch/serialized.joint.json"
    out_discretized = "scratch/dfe.png"
    out_sfs_comparison = "scratch/spectra.png"

spectra = fd.Spectra.from_file(spectra_file)

fixed_params = {'all': {'h': 0.5}}

# choose model
if parametrization == "gamma":
    model = fd.GammaExpParametrization()
elif parametrization == "discrete":
    model = fd.DiscreteFractionalParametrization()
else:
    raise ValueError(f"No parametrization {parametrization}")

# defaults for deleterious-only models
if sub_model == "del":
    if parametrization == "gamma":
        fixed_params["all"] |= {"p_b": 0, "S_b": 1}
    else:  # discrete
        fixed_params["all"] |= {"S4": 1, "S5": 0}

# fixed parameters: fixed_param=value_param=value
elif sub_model.startswith("fixed"):
    fixed_params["all"] |= {
        k: float(v)
        for k, v in re.findall(
            r"_([A-Za-z_]+)=([-+]?[0-9]*\.?[0-9]+)", sub_model.replace("fixed", "")
        )
    }

if not use_eps:
    fixed_params['all']['eps'] = 0

inf = fd.JointInference(
    sfs_neut=spectra['neutral.*'].merge_groups(level=1),
    sfs_sel=spectra['selected.*'].merge_groups(level=1),
    shared_params=[fd.SharedParams(params=shared_params, types='all')],
    model=model,
    fixed_params=fixed_params,
    n_bootstraps=100,
    n_runs=100,
    n_bootstrap_retries=10,
    parallelize=parallelize
)

# perform inference
inf.run()

# save object in serialized form
inf.to_file(out_json)

# plot results
inf.plot_sfs_comparison(file=out_sfs_comparison, show=testing)
inf.plot_discretized(file=out_discretized, show=testing)

pass
