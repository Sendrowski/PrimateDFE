import fastdfe as fd
import numpy as np


class Parametrizer:
    """
    Parametrizer utils.
    """

    discretization = fd.Discretization(n=10)

    @classmethod
    def get_S_d(cls, dfe: fd.DFE) -> np.ndarray:
        """
        Get the mean deleterious selection coefficient from the DFE.
        """
        s = cls.discretization.s
        ds = cls.discretization.interval_sizes
        m = (s < 0)

        x = np.zeros(len(dfe.bootstraps), dtype=float)

        for i, bs in enumerate(dfe.get_bootstrap_dfes()):
            p_del = bs.cdf(0)
            mean_S_del = np.sum(bs.pdf(s)[m] * s[m] * ds[m])

            x[i] = mean_S_del / p_del if p_del > 0 else 0

        return x

    @classmethod
    def get_S_b(cls, dfe: fd.DFE) -> np.ndarray:
        """
        Get the mean beneficial selection coefficient from the DFE.
        """
        s = cls.discretization.s
        ds = cls.discretization.interval_sizes
        m = (s > 0)

        x = np.zeros(len(dfe.bootstraps), dtype=float)

        for i, bs in enumerate(dfe.get_bootstrap_dfes()):
            p_ben = 1 - bs.cdf(0)
            mean_S_ben = np.sum(bs.pdf(s)[m] * s[m] * ds[m])

            x[i] = mean_S_ben / p_ben if p_ben > 0 else 0

        return x

    @staticmethod
    def get_p_b(dfe: fd.DFE) -> float:
        """
        Get the proportion of beneficial mutations from the DFE.
        """
        x = np.zeros(len(dfe.bootstraps))

        for i, bs in dfe.bootstraps.iterrows():
            cdf = dfe.model.get_cdf(**bs)
            x[int(i)] = 1 - cdf(0)

        return x

    @classmethod
    def get_S_range(
            cls,
            dfe: fd.DFE,
            S_min: float,
            S_max: float,
    ) -> np.ndarray:
        """
        Get the proportion of mutations in the given selection coefficient range.
        """
        x = np.zeros(len(dfe.bootstraps), dtype=float)

        for i, bs in enumerate(dfe.get_bootstrap_dfes()):
            x[i] = bs.cdf(S_max) - bs.cdf(S_min)

        return x
