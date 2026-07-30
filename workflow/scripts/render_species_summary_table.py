"""Render a per-subspecies summary table (mu, Ne) as a LaTeX longtable
fragment for the appendix, addressing Reviewer 3's request for a tabulation
of mutation-rate estimates alongside the inferred effective population sizes.

Ne is read directly from the canonical Watterson estimates
(results/stats/Ne/comp/original_ref/{group}/8.csv) rather than the aggregated
population_summary table, so the tabulation cannot drift from the values used
elsewhere in the analysis.
"""
import pandas as pd

try:
    ne_file = snakemake.input.ne
    out = snakemake.output[0]
    pops = list(snakemake.params.populations)
    references = dict(snakemake.params.references)
except NameError:
    import sys
    sys.path.insert(0, 'workflow/scripts')
    from populations import Populations
    ne_file = "results/stats/Ne/comp/original_ref/catarrhini/8.csv"
    out = "scratch/species_summary.tex"
    pops = Populations.get_pops(8, 'catarrhini')

    def _ref_from_pop(pop: str) -> str:
        species = '_'.join(pop.split('_')[:2])
        meta = pd.read_csv(f"resources/metadata/{species.split('_')[0]}_individuals.csv", sep="\t")
        rows = meta[meta["BAM_FOLDER"].str.contains(species)]
        return rows.iloc[0].REFERENCE_FOLDER.replace('_ssp', '')

    references = {p: _ref_from_pop(p) for p in pops}


def species_label(pop: str) -> str:
    parts = pop.replace('_', ' ').split(' ')
    if len(parts) >= 2:
        return (
            r"\textit{" + ' '.join(parts[:2]) + "}"
            + ((' ' + ' '.join(parts[2:])) if len(parts) > 2 else '')
        )
    return r"\textit{" + pop.replace('_', ' ') + "}"


ne_df = pd.read_csv(ne_file)
ne_map = dict(zip(ne_df["label"], ne_df["x"]))
df = pd.DataFrame({'population': pops})
df['reference'] = df['population'].map(references)
df['ne_estimate'] = df['population'].map(ne_map)
df = df[df['ne_estimate'].notna()].copy()
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

# Subspecies whose own species is absent from Kuderna fall back to the reference
# species's mu; compute that list from the data so the caption never drifts.
kuderna_species = set(pd.read_csv('resources/Kuderna/species.csv')['SPECIES_BINOMIAL'])
absent = [p for p in df['population'] if '_'.join(p.split('_')[:2]) not in kuderna_species]
num_words = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five',
             6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten'}
n_absent = num_words.get(len(absent), str(len(absent)))
absent_labels = ', '.join(
    r"\textit{" + p.split('_')[0][0] + r".\ " + p.split('_')[1] + "}" for p in absent
)
caption = (
    r"        \rowcolor{white}\caption{Per-subspecies mutation rate and inferred effective "
    r"population size. Mutation rates $\mu$ (per site per generation) are taken from "
    r"\citet{kuderna2023} at the species level and applied to all subspecies sharing the same "
    r"reference genome; for the "
    + n_absent
    + r" subspecies absent from \citet{kuderna2023} ("
    + absent_labels
    + r"), $\mu$ is the value reported for the reference species onto which they were mapped. "
    + r"Effective population sizes $N_e$ are estimated as $\hat\theta / (4\mu)$ using "
    + r"Watterson's estimator on neutral (4-fold) sites (see Methods). Subspecies are ordered by "
    + r"$N_e$.}\label{tab:species_summary}"
)

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
    r"        \bottomrule" "\n"
    r"        \addlinespace[2ex]" "\n"
    + caption + "\n"
    r"        \endlastfoot" "\n"
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
    r"    \end{longtable}" "\n"
    r"\end{center}" "\n"
)

with open(out, 'w') as f:
    f.write(header)
    f.write("\n".join(rows))
    f.write(footer)
