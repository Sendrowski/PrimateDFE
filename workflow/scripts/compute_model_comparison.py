"""Compute log-likelihood, AIC, and del-vs-full LRT p-values for the four
DFE-model variants (gamma x {full, del}, discrete x {full, del}) per
population, from per-inference BaseInference JSON files.
"""
import fastdfe as fd
import pandas as pd
from scipy.stats import chi2

try:
    inputs = list(snakemake.input)
    out = snakemake.output[0]
    pops = list(snakemake.params.populations)
    refs = list(snakemake.params.refs)
except NameError:
    from populations import Populations
    pops = Populations.get_pops(8, 'catarrhini')
    refs = [Populations.get_ref_from_pop(p) for p in pops]
    inputs = [
        f"results/dfe/{r}/{p}/result.unfolded.8.{para}.{sub}.noeps.json"
        for r, p in zip(refs, pops)
        for para in ('gamma', 'discrete')
        for sub in ('del', 'full')
    ]
    out = "scratch/model_comparison.csv"

K = {
    ('gamma', 'del'): 2,
    ('gamma', 'full'): 4,
    ('discrete', 'del'): 3,
    ('discrete', 'full'): 5,
}


def best_ll(path: str) -> float:
    bi = fd.BaseInference.from_file(path)
    return float(bi.runs['likelihood'].max())


def aic(ll: float, k: int) -> float:
    return 2 * k - 2 * ll


def lrt_p(ll_full: float, ll_del: float, df: int = 2) -> float:
    stat = 2 * (ll_full - ll_del)
    if stat <= 0:
        return 1.0
    return float(chi2.sf(stat, df))


by_pop: dict = {}
for path in inputs:
    fn = path.rsplit('/', 1)[-1]
    pop = path.rsplit('/', 2)[-2]
    para = fn.split('.')[3]
    sub = fn.split('.')[4]
    by_pop.setdefault(pop, {})[(para, sub)] = best_ll(path)

rows = []
for pop in pops:
    d = by_pop[pop]
    ll_gd, ll_gf = d[('gamma', 'del')], d[('gamma', 'full')]
    ll_dd, ll_df = d[('discrete', 'del')], d[('discrete', 'full')]
    rows.append({
        'population': pop,
        'LL_gamma_del': ll_gd,
        'LL_gamma_full': ll_gf,
        'AIC_gamma_del': aic(ll_gd, K[('gamma', 'del')]),
        'AIC_gamma_full': aic(ll_gf, K[('gamma', 'full')]),
        'LRT_p_gamma': lrt_p(ll_gf, ll_gd),
        'LL_discrete_del': ll_dd,
        'LL_discrete_full': ll_df,
        'AIC_discrete_del': aic(ll_dd, K[('discrete', 'del')]),
        'AIC_discrete_full': aic(ll_df, K[('discrete', 'full')]),
        'LRT_p_discrete': lrt_p(ll_df, ll_dd),
        'dAIC_full_gamma_minus_discrete': aic(ll_gf, 4) - aic(ll_df, 5),
    })

pd.DataFrame(rows).to_csv(out, index=False)
