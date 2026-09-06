"""Render the per-subspecies outgroup table as a LaTeX longtable fragment for
the appendix.

Outgroups are read from Populations.ingroup_outgroups, the same mapping the
alignment and ancestral-annotation rules consume, so the tabulation cannot
drift from the outgroups actually used. Subspecies inherit the entry of their
species.
"""
import sys

try:
    out = snakemake.output[0]
    pops = list(snakemake.params.populations)
except NameError:
    sys.path.insert(0, 'workflow/scripts')
    from populations import Populations
    out = "scratch/outgroups.tex"
    pops = Populations.get_pops(8, 'catarrhini')

sys.path.insert(0, 'workflow/scripts')
from populations import Populations  # noqa: E402


def species_of(pop: str) -> str:
    return '_'.join(pop.split('_')[:2])


def italic(name: str) -> str:
    return r"\textit{" + name.replace('_', ' ') + "}"


header = (
    r"\begin{center}" "\n"
    r"    \footnotesize" "\n"
    r"    \renewcommand{\arraystretch}{1.15}" "\n"
    r"    \rowcolors{2}{gray!8}{white}" "\n"
    "\n"
    r"    \begin{longtable}{p{0.34\textwidth} p{0.58\textwidth}}" "\n"
    r"        \toprule" "\n"
    r"        \textbf{Subspecies} & \textbf{Outgroups} \\" "\n"
    r"        \midrule" "\n"
    r"        \endfirsthead" "\n"
    r"        \multicolumn{2}{c}{\tablename~\thetable{} -- continued from previous page} \\" "\n"
    r"        \toprule" "\n"
    r"        \textbf{Subspecies} & \textbf{Outgroups} \\" "\n"
    r"        \midrule" "\n"
    r"        \endhead" "\n"
    r"        \bottomrule" "\n"
    r"        \addlinespace[2ex]" "\n"
    r"        \rowcolor{white}\caption{Outgroups used for each subspecies in the analysis.}"
    r"\label{tab:outgroups}" "\n"
    r"        \endlastfoot" "\n"
    "\n"
)

rows = []
for pop in pops:
    outgroups = Populations.ingroup_outgroups[species_of(pop)]
    rows.append(
        " " * 8 + italic(pop) + " & " + ", ".join(italic(o) for o in outgroups) + r" \\"
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
