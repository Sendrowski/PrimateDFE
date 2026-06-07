"""
Combine multiple stats CSV files into a single summary CSV.

This script reads multiple CSV files containing variant statistics, aggregates 
their values, and writes a combined summary file.
"""

import pandas as pd

try:
    input_files = snakemake.input
    output_file = snakemake.output[0]
except NameError:
    # Testing
    input_files = ['scratch/recoded.stats.csv', 'scratch/recoded.stats.csv']
    output_file = "scratch/combined_stats.tsv"

# Load and concatenate all stats files
df_list = [pd.read_csv(f, sep="\t") for f in input_files]
df_combined = pd.concat(df_list, ignore_index=True)

# Aggregate numeric columns by summing them
agg_stats = df_combined.sum(numeric_only=True)

# Special handling for spectrum dictionary
div_spectra = [eval(d) for d in df_combined["div_spectrum"]]
combined_spectrum = {}
for spectrum in div_spectra:
    for key, value in spectrum.items():
        combined_spectrum[key] = combined_spectrum.get(key, 0) + value

# Recompute fractional stats
agg_stats["div"] = agg_stats["n_divergence"] / agg_stats["n_total"]
agg_stats["masked"] = agg_stats["n_masked"] / agg_stats["n_total"]
agg_stats["complement"] = agg_stats["n_complement"] / agg_stats["n_total"]

# Compute variance for fractional stats
agg_stats["div_var"] = df_combined["div"].var(ddof=1)
agg_stats["masked_var"] = df_combined["masked"].var(ddof=1)
agg_stats["complement_var"] = df_combined["complement"].var(ddof=1)

# Convert dictionary back to string format
agg_stats["div_spectrum"] = str(combined_spectrum)

# Convert aggregated stats to DataFrame and save
agg_stats = pd.DataFrame(agg_stats).T
agg_stats.to_csv(output_file, sep="\t", index=False)

print(agg_stats.T)