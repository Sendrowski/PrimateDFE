"""
Extract effective population size (Ne) estimates from Kuderna et al. (2022) for a set of species.
"""

import pandas as pd

try:
    testing = False
    x_file = snakemake.input.x
    out = snakemake.output[0]
except NameError:
    testing = True
    x_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    out = "scratch/kuderna_ne.csv"


def to_species(name: str) -> str:
    """
    Convert a population name to a species name by removing the last underscore and suffix.
    """
    return name.rsplit("_", 1)[0] if name.count("_") >= 2 else name


x_df = pd.read_csv(x_file)
k_df = pd.read_csv("resources/Kuderna/species.csv")[["SPECIES_BINOMIAL", "EFFECTIVE_POP_SIZE"]] \
    .rename(columns={"SPECIES_BINOMIAL": "species", "EFFECTIVE_POP_SIZE": "Ne"})

merged = (
    x_df.assign(species=x_df["label"].apply(to_species))
        .merge(k_df, on="species", how="left")
)

x_df.loc[merged["Ne"].notna(), "x"] = merged.loc[merged["Ne"].notna(), "Ne"]

x_df[["label", "x"]].to_csv(out, index=False)
