"""Render the model-comparison CSV as a LaTeX table fragment for the appendix."""
import pandas as pd

try:
    csv = snakemake.input[0]
    out = snakemake.output[0]
except NameError:
    csv = "scratch/model_comparison.csv"
    out = "scratch/model_comparison.tex"


def fmt_ll(x: float) -> str:
    return f"{x:.1f}"


def fmt_aic(x: float) -> str:
    return f"{x:.1f}"


def fmt_p(p: float) -> str:
    if p >= 0.01:
        return f"{p:.2f}"
    if p < 1e-300:
        return r"$<10^{-300}$"
    exp = int(f"{p:e}".split('e')[1])
    mant = p / 10 ** exp
    return f"${mant:.1f}\\!\\times\\!10^{{{exp}}}$"


def species_label(pop: str) -> str:
    parts = pop.replace('_', ' ').split(' ')
    if len(parts) >= 2:
        return r"\textit{" + ' '.join(parts[:2]) + "}" + (' ' + ' '.join(parts[2:]) if len(parts) > 2 else '')
    return r"\textit{" + pop.replace('_', ' ') + "}"


df = pd.read_csv(csv)

header = (
    r"\begin{table}[H]" "\n"
    r"    \centering" "\n"
    r"    \rowcolors{2}{gray!8}{white}" "\n"
    r"    \renewcommand{\arraystretch}{1.15}" "\n"
    r"    \resizebox{\textwidth}{!}{%" "\n"
    r"    \begin{tabular}{l rr rr c rr rr c r}" "\n"
    r"        \toprule" "\n"
    r"         & \multicolumn{5}{c}{\textbf{GammaExp}} & \multicolumn{5}{c}{\textbf{Discrete}} &  \\" "\n"
    r"        \cmidrule(lr){2-6} \cmidrule(lr){7-11}" "\n"
    r"        \textbf{Subspecies} & LL$_{\text{del}}$ & LL$_{\text{full}}$ & AIC$_{\text{del}}$ & AIC$_{\text{full}}$ & LRT $p$ & LL$_{\text{del}}$ & LL$_{\text{full}}$ & AIC$_{\text{del}}$ & AIC$_{\text{full}}$ & LRT $p$ & $\Delta$AIC \\" "\n"
    r"        \midrule" "\n"
)

rows = []
for _, r in df.iterrows():
    rows.append(
        " " * 8 + " & ".join([
            species_label(r['population']),
            fmt_ll(r['LL_gamma_del']),
            fmt_ll(r['LL_gamma_full']),
            fmt_aic(r['AIC_gamma_del']),
            fmt_aic(r['AIC_gamma_full']),
            fmt_p(r['LRT_p_gamma']),
            fmt_ll(r['LL_discrete_del']),
            fmt_ll(r['LL_discrete_full']),
            fmt_aic(r['AIC_discrete_del']),
            fmt_aic(r['AIC_discrete_full']),
            fmt_p(r['LRT_p_discrete']),
            f"{r['dAIC_full_gamma_minus_discrete']:+.1f}",
        ]) + r" \\"
    )

# Mean row across subspecies; LRT p-values are not meaningfully averaged
# (they're bounded in [0,1] and combine nonlinearly), so report as em-dash.
mean_cells = [
    r"\textbf{Mean}",
    fmt_ll(df['LL_gamma_del'].mean()),
    fmt_ll(df['LL_gamma_full'].mean()),
    fmt_aic(df['AIC_gamma_del'].mean()),
    fmt_aic(df['AIC_gamma_full'].mean()),
    "---",
    fmt_ll(df['LL_discrete_del'].mean()),
    fmt_ll(df['LL_discrete_full'].mean()),
    fmt_aic(df['AIC_discrete_del'].mean()),
    fmt_aic(df['AIC_discrete_full'].mean()),
    "---",
    f"{df['dAIC_full_gamma_minus_discrete'].mean():+.1f}",
]
rows.append(" " * 8 + r"\midrule")
rows.append(" " * 8 + " & ".join(mean_cells) + r" \\")

footer = (
    "\n"
    r"        \bottomrule" "\n"
    r"    \end{tabular}}" "\n"
    r"    \caption{Per-subspecies model comparison across the four DFE-model variants of Table~\ref{tab:dfe_variants} (\texttt{GammaExpParametrization} and \texttt{DiscreteParametrization}, each deleterious-only (\emph{del}) or full (\emph{full})). $\mathrm{LL}$ is the log-likelihood at the MLE; $\mathrm{AIC} = 2k - 2\,\mathrm{LL}$ with free-parameter counts $k$ as in Table~\ref{tab:dfe_variants}. LRT $p$-values are from a $\chi^2_2$ test of the nested deleterious-only vs.\ full model within each parametrization: small $p$ ($p < 0.05$) indicates that adding the beneficial component significantly improves the fit, while values close to one indicate that the simpler deleterious-only model is preferred. The final column reports $\Delta\mathrm{AIC} = \mathrm{AIC}_{\text{$\gamma$ full}} - \mathrm{AIC}_{\text{discrete full}}$; positive values favor the discrete parametrization. The bottom row reports per-column means across subspecies (means are omitted for the LRT $p$-values, which are not meaningfully averaged). Overall, the deleterious-only model is preferred for 23/38 subspecies (LRT $p \geq 0.05$ in both parametrizations), and the gamma and discrete full-DFE fits are close in AIC (mean $\Delta\mathrm{AIC} = +4.0$ in favor of discrete; $|\Delta\mathrm{AIC}| < 5$ in 21/38 subspecies), consistent with the qualitative agreement between the two parametrizations reported in the main text. Inferences use unfolded SFS at $n=8$ haplotypes with $\varepsilon = 0$.}" "\n"
    r"    \label{tab:model_comparison}" "\n"
    r"\end{table}" "\n"
)

with open(out, 'w') as f:
    f.write(header)
    f.write("\n".join(rows))
    f.write(footer)
