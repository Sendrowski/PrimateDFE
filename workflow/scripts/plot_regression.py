"""
Plot regression plot
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

try:
    testing = False
    x_file = snakemake.input.x
    y_file = snakemake.input.y
    x_col = snakemake.params.x_col
    y_col = snakemake.params.y_col
    x_label = snakemake.params.x_label
    y_label = snakemake.params.y_label
    x_join_key = snakemake.params.x_join_key
    y_join_key = snakemake.params.y_join_key
    include = snakemake.params.get("include", None)
    labels = snakemake.params.get("labels", None)
    title = snakemake.params.title
    figsize = snakemake.params.get("figsize", (6, 4))
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    x_file = "results/stats/Ne/comp/original_ref/all/8.csv"
    y_file = "resources/Kuderna/species.csv"
    x_col = "x"
    y_col = "EFFECTIVE_POP_SIZE"
    x_label = "Castellano"
    y_label = "New"
    x_join_key = "label"
    y_join_key = "SPECIES_BINOMIAL"
    labels = None
    include = None
    labels = [x.split('_')[0] for x in
              ['Homo_sapiens', 'Pan_troglodytes', 'Pan_paniscus', 'Pan_troglodytes', 'Gorilla_gorilla',
               'Gorilla_beringei', 'Pongo_abelii', 'Pongo_pygmaeus', 'Macaca_cyclopis', 'Macaca_fuscata',
               'Macaca_mulatta', 'Macaca_leonina', 'Macaca_maura', 'Macaca_nemestrina', 'Macaca_nigra',
               'Macaca_radiata',
               'Macaca_silenus', 'Macaca_tonkeana', 'Macaca_arctoides', 'Papio_anubis', 'Papio_cynocephalus',
               'Papio_hamadryas', 'Papio_kindae', 'Papio_papio', 'Papio_ursinus', 'Theropithecus_gelada',
               'Erythrocebus_patas', 'Rhinopithecus_bieti', 'Rhinopithecus_roxellana', 'Pygathrix_cinerea',
               'Pygathrix_nemaeus', 'Pygathrix_nigripes', 'Semnopithecus_entellus', 'Semnopithecus_hypoleucos',
               'Trachypithecus_francoisi', 'Trachypithecus_geei', 'Trachypithecus_phayrei']]
    title = "Population"
    figsize = (6, 4)
    out = "scratch/Ne_reg.png"

# load data
x_df = pd.read_csv(x_file)
y_df = pd.read_csv(y_file)

x_df = x_df.add_suffix("_x")
y_df = y_df.add_suffix("_y")

x_df['join_key'] = x_df[f"{x_join_key}_x"]
y_df['join_key'] = y_df[f"{y_join_key}_y"]

# ensure consistent labels
merged = x_df.merge(y_df, on="join_key")

# filter to keys
if include is not None:
    merged = merged[merged["join_key"].isin(include)]
    merged = merged.set_index("join_key").loc[include].reset_index()

if labels:
    cmap = plt.get_cmap("tab10")
    label_set = list(sorted(set(labels)))
    merged['hue'] = [cmap(label_set.index(labels[i]) / len(label_set)) for i in range(len(labels))]

# regression
slope, intercept, r_value, p_value, std_err = stats.linregress(
    merged[f"{x_col}_x"], merged[f"{y_col}_y"]
)

plt.figure(figsize=figsize)
# Color points by population label
plt.scatter(x=merged[f"{x_col}_x"], y=merged[f"{y_col}_y"], color=merged['hue'] if labels else None)

# regression line
x_vals = merged[f"{x_col}_x"]
order = np.argsort(x_vals)
x_sorted = x_vals.iloc[order]
y_sorted = (intercept + slope * x_vals).iloc[order]
reg_line, = plt.plot(
    x_sorted,
    y_sorted,
    color="black",
    linestyle="--",
    label=f"r={r_value:.2f}, p={p_value:.3g}"
)

legend_handles = []
legend_labels = []

if labels:
    unique = {lab: col for lab, col in zip(labels, merged['hue'])}

    for lab, col in unique.items():
        legend_handles.append(
            plt.Line2D([0], [0], marker='o', linestyle='', color=col, markersize=6)
        )
        legend_labels.append(lab)

legend_handles.append(reg_line)
legend_labels.append(reg_line.get_label())

plt.legend(legend_handles, legend_labels, fontsize=9)

plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(title)
plt.tight_layout()

if testing:
    plt.show()

plt.savefig(out)
