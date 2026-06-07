"""
Combine DFE parameters from several runs into a single file.
"""
import fastdfe as fd
import pandas as pd

try:
    testing = False
    dfe_results = snakemake.input
    labels = snakemake.params.labels
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    dfe_results = [
        'results/dfe/Homo_sapiens/Homo_sapiens/sfs.7.json',
        'results/dfe/Homo_sapiens/Homo_sapiens/sfs.7.json',
    ]
    labels = ['Homo_sapiens1', 'Homo_sapiens2']
    out = 'scratch/combined_results.csv'

records = []
for json_file, label in zip(dfe_results, labels):
    inference = fd.BaseInference.from_file(json_file)
    params = inference.bootstraps.mean()
    params['label'] = label
    records.append(params)

df = pd.DataFrame(records).set_index('label')
df.to_csv(out)