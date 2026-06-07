"""
Plot an example DFE statistic using GammaExpParametrization.
"""
import fastdfe as fd
import matplotlib.pyplot as plt
import numpy as np

p = fd.GammaExpParametrization
#x = p.x0
params = p.x0
dfe = p.get_pdf(**params)


x = np.linspace(-50, 10, 1000)  # adjust range to your data
y = dfe(x)

plt.figure(figsize=(6, 8))  # width=12, height=3 → flatter box
plt.plot(x, y, color="blue")
plt.fill_between(x, y, color="blue", alpha=0.3)
plt.xlabel("Selection coefficient (S)")
plt.ylabel("Density")
plt.title(f"GammaExpParametrization PDF\n{params}")
plt.tight_layout()
plt.show()