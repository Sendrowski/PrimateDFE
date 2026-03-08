import pathlib
import subprocess
import tempfile
from typing import Tuple, Literal
import toytree
import numpy as np
import pandas as pd
from scipy.stats import linregress


class Regression:
    """
    Regression methods.
    """

    @staticmethod
    def regress(
            kind: Literal["linear", "phylo"],
            *,
            x: np.ndarray,
            y: np.ndarray,
            tree_file: pathlib.Path | str | None = None,
            pops: np.ndarray | None = None,
    ) -> Tuple[float, float, float, float]:
        """
        Unified regression interface.
        """
        if kind == "linear":
            return Regression.linear(x, y)

        if kind == "phylo":
            if tree_file is None or pops is None:
                raise ValueError("tree_file and pops are required for phylo regression")
            return Regression.phylo(tree_file, pops, x, y)

        raise ValueError(f"Unknown regression type: {kind}")

    @staticmethod
    def phylo(
            tree_file: pathlib.Path | str,
            pops: np.ndarray,
            x: np.ndarray,
            y: np.ndarray
    ) -> Tuple[float, float, float, float]:
        """
        PGLS via nlme::gls with corBrownian or corGrafen (ML).
        Returns: intercept, slope, pval(slope), r(fitted vs observed)
        """
        from populations import Populations

        species = np.array([Populations.to_species(p) for p in pops])

        # average within species (mirror R)
        df = pd.DataFrame({"species": species, 'population': pops, "x": x, "y": y})
        df = (
            df.groupby("species", as_index=False)
            .agg(x=("x", "mean"), y=("y", "mean"))
            .dropna()
        )

        tree = toytree.tree(tree_file)

        tips = tree.get_tip_labels()

        keep_species = [tip for tip in tree.get_tip_labels() if tip in df["species"].values]

        tree = toytree.mod.prune(tree, *keep_species)

        df = df[df["species"].isin(tips)]

        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            data_file = tmp / "data.csv"
            out_file = tmp / "out.txt"
            tree_file_pruned = tmp / "tree.nwk"

            tree.write(tree_file_pruned)

            np.savetxt(
                data_file,
                np.column_stack([df["species"].values, df["x"].values, df["y"].values]),
                delimiter=",",
                header="species,Ne,y",
                comments="",
                fmt="%s",
            )

            r_code = f"""
            library(nlme)
            library(ape)

            d <- read.csv("{data_file}")
            d$species <- factor(d$species, levels=read.tree("{tree_file_pruned}")$tip.label)
            tree <- read.tree("{tree_file_pruned}")

            fit <- gls(y ~ Ne, data=d,
                       correlation=corGrafen(phy=tree, form=~species, value=0.8, fixed=FALSE),
                       method="ML")
            s <- summary(fit)

            r <- cor(d$Ne, d$y, use = "complete.obs")

            writeLines(
              paste(coef(fit)[1], coef(fit)[2],
                    s$tTable["Ne","p-value"],
                    r,
                    sep=","),
              "{out_file}"
            )
            """
            try:
                print("Running R code for phylogenetic regression...")
                subprocess.run(["Rscript", "-e", r_code], check=True)
            except subprocess.CalledProcessError:
                print("R code failed, falling back to linear regression")
                return Regression.linear(x, y)
            else:
                intercept, slope, p, r = map(float, out_file.read_text().split(","))

        return intercept, slope, p, r

    @staticmethod
    def linear(
            x: np.ndarray,
            y: np.ndarray,
    ) -> Tuple[float, float, float, float]:
        """
        Ordinary least squares linear regression.
        """
        res = linregress(x, y)
        return res.intercept, res.slope, res.pvalue, res.rvalue
