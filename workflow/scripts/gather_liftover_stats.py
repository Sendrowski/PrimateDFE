"""
Sum liftover stats across all VCF chunks.
"""

import pandas as pd

try:
    input_files = snakemake.input
    output_file = snakemake.output[0]
except NameError:
    input_files = [
        f"results/vcf/ingroup/Homo_sapiens/Pongo_abelii.liftover.stats.csv" for i in range(5)
    ]
    output_file = "scratch/liftover_combined.csv"

dfs = [pd.read_csv(f) for f in input_files]
total: pd.DataFrame = sum(dfs)

n_lifted = total['matches'] + total['mismatches']
total['fraction_mismatches'] = total['mismatches'] / n_lifted
total['fraction_mismatches_next_base'] = total['mismatches_next_base'] / n_lifted
total['n_sites'] = n_lifted + total['missing_hit'] + total['invalid_lift'] + total['missing_chrom']
total['fraction_lifted'] = n_lifted / total['n_sites']

total.to_csv(output_file, index=False)
print(total)