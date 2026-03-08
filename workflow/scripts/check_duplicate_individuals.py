"""
Check for duplicate individuals in populations
"""

from collections import defaultdict

import numpy as np

from populations import Populations

pops = np.array(Populations.get_pops(8, "catarrhini"))

ind_to_pops = defaultdict(list)

for pop in pops:
    for ind in Populations.get_individuals_from_pop(pop):
        ind_to_pops[ind].append(pop)

# duplicates = individuals appearing in >1 population
duplicates = {ind: ps for ind, ps in ind_to_pops.items() if len(ps) > 1}

print(duplicates)

pass