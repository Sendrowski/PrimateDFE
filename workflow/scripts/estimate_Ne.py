"""
Estimate effective population size from SFS and mutation rate.
"""

import fastdfe as fd
import pandas as pd

try:
    testing = False
    sfs_file = snakemake.input.sfs
    liftover_stats = snakemake.input.get('liftover', None)
    no_lift = snakemake.params.no_lift
    mu = snakemake.params.mu
    key = snakemake.params.get('key', 'neutral')
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    sfs_file = "results/sfs/Homo_sapiens/Homo_sapiens/sfs.unfolded.8.csv"
    liftover_stats = "results/vcf/ingroup/Homo_sapiens/Homo_sapiens.liftover.stats.csv"
    no_lift = True
    metadata = pd.read_csv("resources/Kuderna/species.csv")
    mu = metadata[metadata.SPECIES_BINOMIAL == 'Homo_sapiens'].iloc[0].MU_PER_GENERATION
    key = 'neutral'
    out = "scratch/Ne.txt"

sfs = fd.Spectra.from_file(sfs_file)

if no_lift:
    fraction_lifted = 1
else:
    fraction_lifted = pd.read_csv(liftover_stats).fraction_lifted.iloc[0]

Ne = sfs[key].theta / fraction_lifted / (4 * mu)

with open(out, 'w') as h:
    h.write(str(Ne))
