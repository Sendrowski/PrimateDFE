import pandas as pd

try:
    testing = False
    inputs = snakemake.input
    labels = snakemake.params.labels
    col_name = snakemake.params.get('colname', 'x')
    output = snakemake.output[0]
except NameError:
    testing = True
    inputs = ["scratch/Ne.txt", "scratch/Ne.txt"]
    labels = ["Homo_sapiens", "Pan_troglodytes"]
    col_name = 'Ne'
    output = "scratch/combined_Ne.txt"

records = [
    {"label": label, col_name: float(open(f).read().strip())} for label, f in zip(labels, inputs)
]

pd.DataFrame(records).to_csv(output, index=False)