"""
Perform LRT (df=2) and plot LR vs median bootstrap S_d.
"""
import fastdfe as fd
import pandas as pd
from scipy.stats import chi2
import numpy as np
from matplotlib import pyplot as plt

try:
    testing = False
    sub_file = snakemake.input.sub
    sup_file = snakemake.input.sup
    out_csv = snakemake.output.csv
    out_plot = snakemake.output.plot
except NameError:
    testing = True
    sub_file = "results/dfe/catarrhini/dfe.unfolded.8.gamma.del.noeps.csv"
    sup_file = "results/dfe/catarrhini/dfe.unfolded.8.gamma.full.noeps.csv"
    out_csv = "scratch/lrt.csv"
    out_plot = "scratch/lrt_vs_sd.png"

df_sub = pd.read_csv(sub_file)
df_sup = pd.read_csv(sup_file)

rows = []

for _, rsub in df_sub.iterrows():
    pop = rsub.population
    rsup = df_sup[df_sup.population == pop].iloc[0]

    dfe_sub = fd.DFE.from_json(rsub.json)
    dfe_sup = fd.DFE.from_json(rsup.json)

    lnl_sub = np.median(dfe_sub.bootstraps["likelihood"])
    lnl_sup = np.median(dfe_sup.bootstraps["likelihood"])

    lr = 2 * (lnl_sup - lnl_sub)
    p = chi2.sf(lr, df=2)

    sd = np.median(dfe_sup.bootstraps["S_d"])
    b = np.median(dfe_sup.bootstraps["b"])

    rows.append({
        "population": pop,
        "lnl_sub": lnl_sub,
        "lnl_sup": lnl_sup,
        "lr": lr,
        "pvalue": p,
        "S_d": sd,
        "b": b
    })

df = pd.DataFrame(rows)
df.to_csv(out_csv, index=False)

# plot
plt.scatter(df["S_d"], df["lr"])
plt.xlabel("Median bootstrap $S_d$")
plt.ylabel("LRT statistic (2ΔlnL)")
plt.tight_layout()
plt.show()

# fraction of significant tests
print("Fraction of significant tests:", np.mean(df["pvalue"] < 0.05))

pass