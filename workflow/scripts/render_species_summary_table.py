"""Render a per-subspecies summary table (mu, Ne) as a LaTeX longtable
fragment for the appendix, addressing Reviewer 3's request for a tabulation
of mutation-rate estimates alongside the inferred effective population sizes.
"""
import pandas as pd

try:
    csv = snakemake.input[0]
    out = snakemake.output[0]
    pops = list(snakemake.params.populations)
except NameError:
    import sys
    sys.path.insert(0, 'workflow/scripts')
    from populations import Populations
    csv = "results/tables/population_summary/catarrhini/population_summary.8.csv"
    out = "scratch/species_summary.tex"
    pops = Populations.get_pops(8, 'catarrhini')


def species_label(pop: str) -> str:
    parts = pop.replace('_', ' ').split(' ')
    if len(parts) >= 2:
        return (
            r"\textit{" + ' '.join(parts[:2]) + "}"
            + ((' ' + ' '.join(parts[2:])) if len(parts) > 2 else '')
        )
    return r"\textit{" + pop.replace('_', ' ') + "}"


df = pd.read_csv(csv)
df = df[df.population.isin(pops)].copy()
# A handful of subspecies are missing from Kuderna; estimate_Ne falls back to
# the reference species's mu (see rule estimate_Ne in workflow/Snakefile),
# so we tabulate that mu here too.
kuderna = pd.read_csv('resources/Kuderna/species.csv')[['SPECIES_BINOMIAL', 'MU_PER_GENERATION']]
df = df.merge(
    kuderna.rename(columns={'SPECIES_BINOMIAL': 'reference', 'MU_PER_GENERATION': 'mu_ref'}),
    on='reference', how='left',
)
df['mu_e8'] = df['mu_ref'] * 1e8
df['ne_e4'] = df['ne_estimate'] / 1e4
df = df.sort_values('ne_estimate').reset_index(drop=True)

header = (
    r"\begin{center}" "\n"
    r"    \footnotesize" "\n"
    r"    \renewcommand{\arraystretch}{1.1}" "\n"
    r"    \rowcolors{2}{gray!8}{white}" "\n"
    "\n"
    r"    \begin{longtable}{l r r}" "\n"
    r"        \toprule" "\n"
    r"        \textbf{Subspecies} & $\boldsymbol{\mu\,(10^{-8}\text{/site/gen})}$ & $\boldsymbol{N_e\,(10^{4})}$ \\" "\n"
    r"        \midrule" "\n"
    r"        \endfirsthead" "\n"
    r"        \multicolumn{3}{c}{\tablename\ \thetable{} -- continued from previous page} \\" "\n"
    r"        \toprule" "\n"
    r"        \textbf{Subspecies} & $\boldsymbol{\mu\,(10^{-8}\text{/site/gen})}$ & $\boldsymbol{N_e\,(10^{4})}$ \\" "\n"
    r"        \midrule" "\n"
    r"        \endhead" "\n"
)

rows = []
for _, r in df.iterrows():
    rows.append(
        " " * 8 + " & ".join([
            species_label(r['population']),
            f"{r['mu_e8']:.2f}",
            f"{r['ne_e4']:.2f}",
        ]) + r" \\"
    )

footer = (
    "\n"
    r"        \bottomrule" "\n"
    r"    \end{longtable}" "\n"
    r"    \captionof{table}{Per-subspecies mutation rate and inferred effective population size. Mutation rates $\mu$ (per site per generation) are taken from \citet{kuderna2023} at the species level and applied to all subspecies sharing the same reference genome; for the five subspecies absent from \citet{kuderna2023} (\textit{M.\ hecki}, \textit{M.\ assamensis}, \textit{M.\ leucogenys}, \textit{R.\ brelichi}, \textit{T.\ poliocephalus}), $\mu$ is the value reported for the reference species onto which they were mapped. Effective population sizes $N_e$ are estimated as $\hat\theta_W / (4\mu)$ using Watterson's estimator on neutral (4-fold) sites (see Methods). Subspecies are ordered by $N_e$.}" "\n"
    r"    \label{tab:species_summary}" "\n"
    r"\end{center}" "\n"
)

with open(out, 'w') as f:
    f.write(header)
    f.write("\n".join(rows))
    f.write(footer)
