import pandas as pd

try:
    testing = False
    input_file = snakemake.input[0]
    output = snakemake.output[0]
except NameError:
    testing = True
    input_file = "results/tables/population_summary/catarrhini/population_summary.8.csv"
    output = "scratch/subspecies_outgroups.csv"

df = pd.read_csv(input_file, usecols=[0, 1], header=0)
df.columns = ["subspecies", "outgroups"]
df.to_csv(output, index=False)
