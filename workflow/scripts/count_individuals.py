"""
Count the number of individuals in the dataset.
"""
import pandas as pd

from populations import Populations

print(Populations.get_pops(8, 'catarrhini'))

df = pd.Series(
    {p: len(Populations.get_individuals_from_pop(p))
     for p in Populations.get_pops(8, 'catarrhini')}
).to_frame(name="n")

pass