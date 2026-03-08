import re
from collections import defaultdict
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import toytree
from matplotlib.axes import Axes
from matplotlib.ticker import LogLocator, FuncFormatter, ScalarFormatter

from parametrizer import Parametrizer
from populations import Populations
from regression import Regression


class DFEvsNePlotter:
    """
    End-to-end engine for plotting DFE-derived statistics against effective
    population size.

    This class owns:
      * taxonomic grouping
      * colors and legend ordering
      * bootstrap extraction from DFEs
      * summary statistics
      * regression and plotting
    """

    def __init__(
            self,
            dfes,
            ne_dict,
            populations,
            stat_list,
            labels,
            label_color=None,
            tree_file=None,
            reg_type: Literal["phylo", "linear"] = "phylo",
    ):
        """
        :param dfes: Mapping ``population -> fastdfe.DFE``
        :param ne_dict: Mapping ``population -> N_e``
        :param populations: Population identifiers
        :param stat_list: Statistics to compute and plot
        :param labels: Optional group labels (otherwise inferred via Populations)
        :param tree_file: Newick file for phylogenetic tree
        :param label_color: Optional mapping ``label -> color``
        :param reg_type: Regression type, either 'phylo' or 'linear'
        """
        self.dfes = dfes
        self.ne = ne_dict
        self.tree_file = tree_file
        self.populations = populations
        self.stat_list = stat_list
        self.labels = labels
        self.reg_type = reg_type

        if label_color is None:
            self.label_color = Populations.get_label_color_map(self.labels)
        else:
            self.label_color = label_color

        self.label_order = sorted(set(self.labels), key=Populations.get_label_rank)

        self._boot = defaultdict(dict)
        self._stats = defaultdict(dict)

    def plot(
            self,
            file=None,
            fig=None,
            show=True,
            show_legend=True,
            qlo=5,
            qhi=95,
            legend_n_cols=None,
            style: Literal["default", "stacked", "dfe"] = "default",
            title: str = ''
    ) -> plt.Figure:
        """
        Compute, plot, and optionally save.

        :param file: If given, save figure to this path
        :param fig: If given, plot on this figure
        :param show: Whether to display the figure
        :param show_legend: Whether to draw the taxonomy legend
        :param qlo: Lower percentile for uncertainty
        :param qhi: Upper percentile for uncertainty
        :param legend_n_cols: Number of columns for legend (if None, auto-compute)
        :param style: Plot style, either 'default' or 'stacked'
        :param title: Overall figure title (only for 'default' style)
        :return: (figure, axes)
        """
        self._boot.clear()
        self._stats.clear()

        self._compute_bootstrap(self.stat_list)
        self._summarize(qlo, qhi)

        stat_labels = [self.make_label(s) for s in self.stat_list]

        if style == "default":
            fig = self._plot_axes_default(fig, self.stat_list, stat_labels, show_legend, legend_n_cols)
        elif style == "stacked":
            fig = self._plot_axes_stacked(fig, self.stat_list, stat_labels, show_legend, legend_n_cols)
        elif style == "dfe":
            fig = self._plot_axes_dfe(fig, self.stat_list, stat_labels, show_legend, legend_n_cols)
        else:
            raise ValueError(style)

        #plt.title(title, y=-0.25)

        if file is not None:
            fig.savefig(file)

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    @staticmethod
    def log_label_pow(x, pos):
        """
        Format log-scaled x-axis labels as powers of 10, with 1, 2, and 5 multiples.
        """
        if x <= 0:
            return ""
        k = int(np.floor(np.log10(x)))
        n = x / 10 ** k
        if np.isclose(n, 1.0):
            return rf"$10^{{{k}}}$"
        return rf"${int(round(n))}\times10^{{{k}}}$"

    def plot_two_datasets_stacked(
            self,
            datasets: dict,
            file=None,
            fig=None,
            show=True,
            qlo=5,
            qhi=95,
    ) -> plt.Figure:
        """
        Plot two datasets in stacked style, with a shared legend.
        """
        if fig is None:
            fig = plt.subplots(
                len(self.stat_list) + 1,
                1,
                figsize=(10, 0.8 * len(self.stat_list) + 1),
                sharex=True,
                dpi=400
            )[0]
        axes = fig.axes

        legend_handles = {}

        for name, (dfes, ne_dict, labels, color) in datasets.items():
            # temporarily swap in dataset-specific state
            self.dfes = dfes
            self.ne = ne_dict
            self.labels = labels
            self.label_color = defaultdict(lambda: color)

            self._boot.clear()
            self._stats.clear()
            self._compute_bootstrap(self.stat_list)
            self._summarize(qlo, qhi)

            for ax, stat in zip(axes[:len(self.stat_list)], self.stat_list):
                for pop in self.populations:
                    m, lo, hi = self._stats[stat][pop]
                    x = self.ne[pop]

                    ax.scatter(x, m, color=color, s=60, alpha=0.6)
                    ax.errorbar(
                        x, m,
                        yerr=[[max(lo, 0)], [max(hi, 0)]],
                        fmt="none", ecolor=color, capsize=3, alpha=0.6
                    )

                self._add_regression(ax, stat, color=color)

            # one proxy handle per dataset
            legend_handles[name] = plt.Line2D(
                [0], [0], marker="o", linestyle="",
                color=color, label=name
            )

        # y-axis labels
        for ax, stat in zip(axes[:len(self.stat_list)], self.stat_list):
            ax.set_ylabel(
                self.make_label(stat),
                fontsize=15,
                rotation=0,
                labelpad=60,
                va="center",
                ha="center",
            )
            ax.ticklabel_format(style="plain", axis="x")

        # x-axis only on bottom data panel
        for ax in axes[:len(self.stat_list) - 1]:
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

        bottom = axes[len(self.stat_list) - 1]
        bottom.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
        bottom.set_xlabel("$N_e$")

        # legend panel
        axes[-1].axis("off")
        axes[-1].legend(
            legend_handles.values(),
            legend_handles.keys(),
            loc="lower center",
            fontsize=9,
            ncol=len(legend_handles),
        )

        fig.tight_layout()
        fig.subplots_adjust(hspace=0.0)

        ax = fig.axes[-1]

        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(FuncFormatter(self.log_label_pow))

        if file is not None:
            fig.savefig(file)

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    def plot_alpha_three_datasets_stacked(
            self,
            datasets: dict,
            file=None,
            fig=None,
            show=True,
            qlo=5,
            qhi=95,
    ) -> plt.Figure:
        """
        Plot alpha vs Ne as three vertically stacked subplots
        (one per dataset) using clade-based colouring.
        """

        if len(datasets) != 3:
            raise ValueError("Exactly three datasets required.")

        stat = "alpha"

        if fig is None:
            fig = plt.subplots(
                4, 1,  # 3 data panels + 1 legend
                figsize=(7, 4.5),
                sharex=True,
                dpi=400
            )[0]

        axes = fig.axes
        data_axes = axes[:3]
        ax_leg = axes[-1]

        legend_handles = {}

        # global Ne range
        all_ne = np.concatenate([
            np.array(list(ne_dict.values()))
            for (_, (_, ne_dict, _, _)) in datasets.items()
        ])
        xmin, xmax = all_ne.min(), all_ne.max()

        for ax, (name, (dfes, ne_dict, labels, _)) in zip(data_axes, datasets.items()):

            self.dfes = dfes
            self.ne = ne_dict
            self.labels = labels

            # restore clade color mapping
            self.label_color = Populations.get_label_color_map(self.labels)

            self._boot.clear()
            self._stats.clear()
            self._compute_bootstrap([stat])
            self._summarize(qlo, qhi)

            for pop, lab in zip(self.populations, self.labels):
                m, lo, hi = self._stats[stat][pop]
                x = self.ne[pop]

                h = ax.scatter(x, m,
                               color=self.label_color[lab],
                               s=60,
                               alpha=0.8)

                if lab not in legend_handles:
                    legend_handles[lab] = h

                ax.errorbar(
                    x, m,
                    yerr=[[max(lo, 0)], [max(hi, 0)]],
                    fmt="none",
                    ecolor=self.label_color[lab],
                    capsize=3,
                    alpha=0.8,
                )

            self._add_regression(ax, stat)
            ax.legend(fontsize=8, loc="lower right")

            ax.set_ylabel(r"$\alpha$", fontsize=14)
            ax.set_title(name, fontsize=11)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-0.1, 1.1)

        # shared x-axis formatting (same as other plots)
        bottom = data_axes[-1]
        bottom.set_xscale("log")
        bottom.set_xlabel(r"$N_e$")
        bottom.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        bottom.xaxis.set_major_formatter(FuncFormatter(self.log_label_pow))

        for ax in data_axes[:-1]:
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

        # taxonomy legend (same logic as stacked DFE plot)
        ax_leg.axis("off")
        ax_leg.legend(
            [legend_handles[l] for l in self.label_order],
            [l.replace("_", " ") for l in self.label_order],
            loc="center",
            fontsize=8,
            ncol=5,
        )

        fig.tight_layout()
        fig.subplots_adjust(hspace=0)

        if file is not None:
            fig.savefig(file)

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    def _compute_bootstrap(self, stat_list):
        """
        Compute bootstrap distributions for all requested statistics.

        :param stat_list: Statistics to extract from each DFE
        """
        for pop in self.populations:
            dfe = self.dfes[pop]
            Ne = self.ne[pop]

            S_d = np.abs(Parametrizer.get_S_d(dfe))
            S_b = np.abs(Parametrizer.get_S_b(dfe))
            p_b = Parametrizer.get_p_b(dfe)

            for stat in stat_list:
                if stat == "S_d":
                    self._boot[stat][pop] = S_d
                elif stat == "s_d":
                    self._boot[stat][pop] = S_d / (4 * Ne)
                elif stat == "S_b":
                    self._boot[stat][pop] = S_b
                elif stat == "p_b":
                    self._boot[stat][pop] = p_b
                elif stat == "S_b*p_b":
                    self._boot[stat][pop] = S_b * p_b
                elif stat in ['h', 'b', 'alpha']:
                    self._boot[stat][pop] = dfe.bootstraps[stat]
                elif stat.startswith("range_S_"):
                    lo, hi = self._parse_range(stat)
                    self._boot[stat][pop] = Parametrizer.get_S_range(dfe, lo, hi)
                elif stat.startswith("range_s_"):
                    lo, hi = self._parse_range(stat)
                    self._boot[stat][pop] = Parametrizer.get_S_range(dfe, lo * 4 * Ne, hi * 4 * Ne)
                else:
                    raise ValueError(f"Unknown statistic: {stat}")

    def _summarize(self, qlo, qhi):
        """
        Summarize bootstrap distributions into medians and percentile intervals.

        :param qlo: Lower percentile
        :param qhi: Upper percentile
        """
        for stat in self._boot:
            for pop, vals in self._boot[stat].items():
                m = np.median(vals)
                lo = m - np.percentile(vals, qlo)
                hi = np.percentile(vals, qhi) - m
                self._stats[stat][pop] = (m, lo, hi)

    def _plot_axes_stacked(self, fig, stat_list, stat_labels, show_legend, legend_n_cols=None) -> plt.Figure:
        """
        Draw all statistic panels and the consolidated legend.

        :param axes: Matplotlib axes
        :param stat_list: Statistics to plot
        :param stat_labels: Axis titles
        :param show_legend: Whether to show taxonomy legend
        """
        if fig is None:
            fig = plt.subplots(
                len(stat_list) + 1,
                1,
                figsize=(10, 0.8 * len(stat_list) + 1),
                sharex=True,
                dpi=400
            )[0]
        axes = fig.axes

        legend_handles = {}

        for ax, stat, title in zip(axes[:len(stat_list)], stat_list, stat_labels):
            for pop, lab in zip(self.populations, self.labels):
                m, lo, hi = self._stats[stat][pop]
                x = self.ne[pop]

                h = ax.scatter(x, m, color=self.label_color[lab], s=60, alpha=0.8)
                if stat == stat_list[0] and lab not in legend_handles:
                    legend_handles[lab] = h

                ax.errorbar(x, m, yerr=[[max(lo, 0)], [max(hi, 0)]], fmt="none",
                            ecolor=self.label_color[lab], capsize=3, alpha=0.8)

            self._add_regression(ax, stat)
            ax.legend(fontsize=8, loc="right")
            ax.set_ylabel(title, fontsize=15, rotation=0, labelpad=60, va="center", ha="center")

            ax.ticklabel_format(style="plain", axis="x")

        axes[-2].set_xlabel("$N_e$")
        axes[-1].axis("off")
        if show_legend:
            axes[-1].legend(
                [legend_handles[l] for l in self.label_order],
                [l.replace("_", " ") for l in self.label_order],
                loc="lower center",
                fontsize=8,
                ncol=legend_n_cols or len(self.label_order) // 10 + 1,
            )

        for ax in axes[:len(stat_list) - 1]:
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

        fig.tight_layout()
        fig.subplots_adjust(hspace=0.0)

        bottom = axes[len(stat_list) - 1]
        bottom.tick_params(axis="x", which="both", bottom=True, labelbottom=True)

        ax = fig.axes[-1]
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(FuncFormatter(self.log_label_pow))

        return fig

    def _plot_axes_dfe(self, fig, stat_list, stat_labels, show_legend, legend_n_cols=None):

        pops_sorted = sorted(self.populations, key=lambda p: self.ne[p])
        n = len(pops_sorted)

        if fig is None:
            fig, axes = plt.subplots(
                1, n,
                figsize=(0.2 * n, 4),
                #sharey=True,
                dpi=400
            )
        else:
            axes = fig.axes

        if n == 1:
            axes = [axes]

        centers = np.arange(len(stat_list))

        for ax, pop in zip(axes, pops_sorted):
            ax.invert_yaxis()

            vals = np.array([
                np.median(self._boot[stat][pop])
                for stat in stat_list
            ])

            ax.barh(
                centers,
                vals,
                color="steelblue",
                height=0.8
            )

            ax.set_title(pop.replace("_", " "), fontsize=7, rotation=90, fontstyle="italic")

            ax.set_xticks([])
            ax.set_xlim(0, 1)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        axes[0].set_yticks(centers)
        axes[0].set_yticklabels([l.replace('S \\in', '') for l in stat_labels], fontsize=9)
        axes[0].tick_params(axis="y", left=True, labelleft=True)

        for ax in axes[1:]:
            ax.set_yticks([])

        axes[0].set_ylabel(r"$S$", rotation=0, fontsize=12, labelpad=10)

        fig.tight_layout()
        fig.subplots_adjust(top=0.55, wspace=0)

        return fig

    def _plot_axes_default(self, fig, stat_list, stat_labels, show_legend, legend_n_cols=None) -> plt.Figure:
        """
        Draw all statistic panels and the consolidated legend.

        :param axes: Matplotlib axes
        :param stat_list: Statistics to plot
        :param stat_labels: Axis titles
        :param show_legend: Whether to show taxonomy legend
        """
        n = len(self.stat_list)

        if fig is None:
            fig, _ = plt.subplots(int(np.ceil((n + 1) / 2)), 2, figsize=(8, 7), dpi=400)

        axes = fig.axes

        legend_handles = {}

        for ax, stat, title in zip(axes[:len(stat_list)], stat_list, stat_labels):
            for pop, lab in zip(self.populations, self.labels):
                m, lo, hi = self._stats[stat][pop]
                x = self.ne[pop]

                h = ax.scatter(x, m, color=self.label_color[lab], s=60, alpha=0.8)
                if stat == stat_list[0] and lab not in legend_handles:
                    legend_handles[lab] = h

                ax.errorbar(x, m, yerr=[[max(lo, 0)], [max(hi, 0)]], fmt="none",
                            ecolor=self.label_color[lab], capsize=3, alpha=0.8)

            self._add_regression(ax, stat)
            ax.set_title(title, fontsize=15)

            # force scientific notation for S_d panels
            if stat in ("S_d",):
                ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
                ax.yaxis.get_major_formatter().set_useMathText(True)

            # show x-axis only on the bottom data row (just before legend panel)
            bottom_axes = axes[-3:-1]  # last two data panels
            if ax in bottom_axes:
                ax.set_xlabel("$N_e$")
            else:
                ax.set_xlabel("")
                ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

            ax.ticklabel_format(style="plain", axis="x")
            ax.legend(fontsize=8)

        axes[-1].axis("off")
        if show_legend:
            axes[-1].legend(
                [legend_handles[l] for l in self.label_order],
                [l.replace("_", " ") for l in self.label_order],
                loc="center",
                fontsize=8,
                ncol=legend_n_cols or len(self.label_order) // 10 + 1,
            )

        for ax in fig.axes:
            ax.set_xscale("log")

            ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.get_major_formatter().set_scientific(False)
            ax.xaxis.get_major_formatter().set_useMathText(False)

        fig.tight_layout()

        return fig

    def _add_regression(self, ax: Axes, stat: str, color="black"):
        """
        Add phylogenetic regression line to axis.
        """
        Nes = np.array([self.ne[p] for p in self.populations])
        meds = np.array([self._stats[stat][p][0] for p in self.populations])
        pops = np.array(self.populations)

        # fit on log(Ne)
        intercept, slope, p, r = Regression.regress(
            self.reg_type,
            x=np.log(Nes),
            y=meds,
            pops=pops,
            tree_file=self.tree_file,
        )

        xx = np.linspace(np.log(Nes.min()), np.log(Nes.max()), 200)
        yy = intercept + slope * xx

        label = f"r={r:.2f}, p={p:.3g}" if not np.isnan(p) else ''

        ax.plot(np.exp(xx), yy, "--", color=color, label=label)

    @staticmethod
    def _format_bound(x, is_s):
        if np.isneginf(x):
            return r"\text{-}\infty"
        if np.isposinf(x):
            return r"\infty"
        if is_s and x != 0:
            k = int(np.log10(abs(x)))
            return rf"\text{{{'-' if x < 0 else ''}}}10^{{\text{{-}}{abs(k)}}}"
        else:
            return rf"\text{{{'-' if x < 0 else ''}}}{abs(x):g}"

    @classmethod
    def make_label(cls, key):
        if key in ('alpha'):
            return f"$\{key}$"
        elif key in ("S_d", "s_d", "S_b", "p_b", "h", "b"):
            return f"${key}$"
        elif key == "S_b*p_b":
            return r"$S_b p_b$"

        m = re.match(r"range_(S|s)_(.+)_(.+)", key)
        if not m:
            raise ValueError(key)

        var, lo, hi = m.groups()
        lo = -np.inf if lo == "-inf" else float(lo)
        hi = np.inf if hi == "inf" else float(hi)

        lo_s = cls._format_bound(lo, var == "s")
        hi_s = cls._format_bound(hi, var == "s")

        left = "(" if not np.isneginf(lo) else "("
        right = "]" if not np.isposinf(hi) else ")"

        return rf"${var} \in {left}{lo_s}, {hi_s}{right}$"

    @staticmethod
    def _parse_range(name):
        """
        Parse a ``range_S_<lo>_<hi>`` statistic key.

        :param name: Statistic key
        :return: (lo, hi)
        """
        _, _, lo, hi = name.split("_")
        lo = -np.inf if lo == "-inf" else float(lo)
        hi = np.inf if hi == "inf" else float(hi)
        return lo, hi
