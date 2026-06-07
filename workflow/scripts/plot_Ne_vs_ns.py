"""
Plot regression plot
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

try:
    testing = False
    x_file = snakemake.input.x
    y_file = snakemake.input.y
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    x_file = "results/stats/Ne/comp/original_ref/all/8.csv"
    y_file = "results/stats/comp/ns/ingroup/population/ns_ratios.csv"
    out = "scratch/Ne_reg.png"

x_col = "x"
y_col = "ns_ratio"
x_label = "$N_e$"
y_label = "$N/S$"
x_join_key = "label"
y_join_key = "species"
title = "$N_e$ vs $N/S$"

# load data
x_df = pd.read_csv(x_file)
y_df = pd.read_csv(y_file)

x_df['join_key'] = x_df[x_join_key]
y_df['join_key'] = y_df[y_join_key]

# ensure consistent labels
merged = pd.merge(x_df, y_df, on='join_key')

# regression
slope, intercept, r_value, p_value, std_err = stats.linregress(
    merged[x_col], merged[y_col]
)

plt.figure(figsize=(6, 6))

# regression line
sns.regplot(data=merged, x=x_col, y=y_col, scatter=False, ci=None, line_kws={"color": "black", "linestyle": "--"})
# Color points by species label
sns.scatterplot(data=merged, x=x_col, y=y_col)

# Annotate each point with its species label
for _, row in merged.iterrows():
    plt.text(row[x_col], row[y_col], row["species"], fontsize=6, ha="left", va="bottom")

plt.text(0.05, 0.95, f"r={r_value:.2f}\np={p_value:.3g}", transform=plt.gca().transAxes, va="top")

plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(title)

if testing:
    plt.show()

plt.savefig(out)
