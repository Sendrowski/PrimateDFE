"""
Test bootstrapping inference on SFS data.
"""

import fastdfe as fd
import numpy as np
from tqdm import tqdm

try:
    testing = False
    sfs_file = snakemake.input[0]
    parametrization = snakemake.params.parametrization
except NameError:
    # testing
    testing = True
    sfs_file = "results/sfs/Homo_sapiens/Homo_sapiens/sfs.unfolded.8.csv"

spectra = fd.Spectra.from_file(sfs_file)
fixed_params = dict(all=dict(eps=0, S4=1, S5=0))
# fd.logger.setLevel('DEBUG')

inf_bs = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    fixed_params=fixed_params,
    model=fd.DiscreteFractionalParametrization(),
    do_bootstrap=True,
    parallelize=True,
    n_runs=10,
    n_bootstrap_retries=1,
    n_bootstraps=100,
    bootstrap_global_mode=True,
    seed=2
)

inf_bs.run()
inf_bs.plot_discretized()

pass

var = np.var([-r.fun for r in inf_bs.bootstrap_results])

inf = fd.BaseInference(
    sfs_neut=spectra['neutral'],
    sfs_sel=spectra['selected'],
    fixed_params=fixed_params,
    do_bootstrap=True,
    n_runs=10,
    n_bootstraps=2
)

inf.run()

for i in tqdm(range(100)):
    bs = fd.BaseInference(
        sfs_neut=spectra['neutral'].resample(seed=2 * i),
        sfs_sel=spectra['selected'].resample(seed=2 * i + 1),
        fixed_params=fixed_params,
        discretization=inf_bs.discretization,
        optimization=inf_bs.optimization,
        do_bootstrap=False,
        n_runs=1,
        parallelize=False,
        x0=dict(all=inf_bs.params_mle),
    )
    bs.run()

    inf.bootstraps.loc[i + 1] = bs.params_mle | dict(alpha=0, likelihood=inf.likelihood)

inf.plot_discretized()

pass
