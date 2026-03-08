"""
Compare two nested DFE models on the same SFS data.
"""
import fastdfe as fd
from matplotlib import pyplot as plt
import re
import pandas as pd

try:
    testing = False
    inf_sub_file = snakemake.input.sub
    inf_sup_file = snakemake.input.sup
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    inf_sup_file = "results/dfe/Homo_sapiens/Homo_sapiens/dfe.unfolded.8.discrete.full.noeps.json"
    inf_sub_file = "results/dfe/Homo_sapiens/Homo_sapiens/dfe.unfolded.8.discrete.del.noeps.json"
    out = "scratch/lrt.csv"

inf_sub = fd.BaseInference.from_file(inf_sub_file)
inf_sup = fd.BaseInference.from_file(inf_sup_file)

p = inf_sub.compare_nested(inf_sup)

pd.DataFrame({
    "pvalue": [p],
    "lnl_sub": [inf_sub.likelihood],
    "lnl_sup": [inf_sup.likelihood]
}).to_csv(out, index=False)

pass
