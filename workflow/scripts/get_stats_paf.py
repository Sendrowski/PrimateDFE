"""
Get alignment lengths from a PAF file.
"""

import pandas as pd

try:
    paf_file = snakemake.input[0]
    out = snakemake.output[0]
except NameError:
    # Testing
    paf_file = "results/minimap2/hg38/Gorilla_gorilla_gorilla.paf"
    out = "scratch/lengths.csv"

# Read the PAF file without column names
paf_df = pd.read_csv(paf_file, sep="\t", header=None, on_bad_lines='skip')

# Column 11: block length
block_lengths = paf_df[10]

# Compute sorted block lengths and N90/N99
sorted_lengths = block_lengths.sort_values(ascending=False).reset_index(drop=True)
cumulative = sorted_lengths.cumsum()
total = sorted_lengths.sum()
n90 = sorted_lengths[cumulative >= 0.9 * total].iloc[0]
n99 = sorted_lengths[cumulative >= 0.99 * total].iloc[0]

# Extract target reference intervals: ref name (5), start (7), end (8)
intervals = paf_df[[5, 7, 8]].copy()
intervals.columns = ["ref", "start", "end"]

# Sort intervals by reference and start
intervals.sort_values(by=["ref", "start"], inplace=True)

# Merge overlapping intervals and calculate total non-overlapping aligned length
merged_lengths = []

# Group intervals by reference name
for _, group in intervals.groupby("ref"):
    current_start, current_end = None, None

    # Iterate through each alignment (sorted by start)
    for row in group.itertuples(index=False):
        start, end = row.start, row.end

        if current_start is None:
            # First interval in group
            current_start, current_end = start, end
        elif start <= current_end:
            # Overlapping or adjacent interval → extend the current interval
            current_end = max(current_end, end)
        else:
            # Non-overlapping → save current interval length and start a new one
            merged_lengths.append(current_end - current_start)
            current_start, current_end = start, end

    # After loop, save the last interval
    if current_start is not None:
        merged_lengths.append(current_end - current_start)

# Save results
summary = pd.DataFrame({
    "total_length": [total],
    "overlap_length": [total - sum(merged_lengths)],
    "overlap_fraction": [(total - sum(merged_lengths)) / total],
    "N90": [n90],
    "N99": [n99],
})
summary.to_csv(out, index=False)