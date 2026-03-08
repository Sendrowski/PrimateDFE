"""
Plot discretized DFE for a given population from a joint inference JSON file.
"""

import fastdfe as fd

try:
    testing = False
    json_file = snakemake.input[0]
    population = snakemake.params.population
    mode = snakemake.params.mode
    kind = snakemake.params.kind
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    json_file = "results/dfe/Castellano/joint/Pongo_abelii/b/unfolded.8.del.noeps.json"
    population = "Pongo_abelii"
    mode = "marginal"
    kind = "discretized"
    out = "scratch/dfe_Ne.png"

joint = fd.JointInference.from_file(json_file)
inf_dict = joint.marginal_inferences if mode == "marginal" else joint.joint_inferences
inf = inf_dict[population]

if kind == "discretized":
    title = ', '.join([f"{k}={v}" for k, v in inf.bootstraps.mean().round(2).items()])
    inf.plot_discretized(file=out, show=testing, title=title)
else:
    inf.plot_sfs_comparison(file=out, show=testing)

pass
