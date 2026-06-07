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

    @staticmethod
    def _cdf_on(bs, x) -> np.ndarray:
        """
        Evaluate a bootstrap DFE's scaled CDF on an array of S values,
        vectorising when supported and falling back to a per-point loop.
        """
        try:
            y = np.asarray(bs.cdf(x), dtype=float)
            if y.shape == np.shape(x):
                return y
        except (TypeError, ValueError):
            pass
        return np.array([bs.cdf(float(v)) for v in np.atleast_1d(x)], dtype=float)

    @classmethod
    def get_cdf_at(cls, dfe: fd.DFE, S: float) -> np.ndarray:
        """
        Get the cumulative scaled-DFE probability P[S' < S] per bootstrap.
        Used for the iso-threshold epistasis test, where S = 4 N_e s* for a
        fixed unscaled threshold s*.
        """
        x = np.zeros(len(dfe.bootstraps), dtype=float)

        for i, bs in enumerate(dfe.get_bootstrap_dfes()):
            x[i] = bs.cdf(S)

        return x

    @classmethod
    def get_S_quantiles(cls, dfe: fd.DFE, fs) -> dict:
        """
        Get deleterious-mass quantiles S_q for several fractions ``fs`` at once,
        reusing the per-bootstrap CDF evaluation. Returns {f: array}.
        """
        neg = cls.discretization.s[cls.discretization.s < 0]  # ascending, < 0
        n = len(dfe.bootstraps)
        out = {f: np.full(n, np.nan) for f in fs}

        for i, bs in enumerate(dfe.get_bootstrap_dfes()):
            p_del = bs.cdf(0)
            if p_del <= 0:
                continue
            cdf_neg = cls._cdf_on(bs, neg)
            for f in fs:
                target = f * p_del
                if cdf_neg[0] <= target <= cdf_neg[-1]:
                    out[f][i] = np.interp(target, cdf_neg, neg)

        return out

    @classmethod
    def get_S_quantile(cls, dfe: fd.DFE, f: float) -> np.ndarray:
        """
        Get the scaled coefficient S_q below which fraction ``f`` of the
        *deleterious* probability mass lies, per bootstrap (f in (0, 1)).

        This is a tail-robust alternative to the mean S_d: a quantile rescales
        linearly with N_e under no epistasis (slope 1 of log|S_q| vs log N_e),
        whereas the mean is dominated by the unobserved strong/lethal tail.
        """
        neg = cls.discretization.s[cls.discretization.s < 0]  # ascending, < 0

        x = np.zeros(len(dfe.bootstraps), dtype=float)

        for i, bs in enumerate(dfe.get_bootstrap_dfes()):
            p_del = bs.cdf(0)
            if p_del <= 0:
                x[i] = np.nan
                continue

            cdf_neg = cls._cdf_on(bs, neg)  # ascending in s
            target = f * p_del
            # guard against clamping at the grid edges: a target beyond the
            # near-zero resolution (or below the most-negative point) is not
            # resolvable and would otherwise pin S_q to the grid floor.
            if target > cdf_neg[-1] or target < cdf_neg[0]:
                x[i] = np.nan
                continue
            # invert: find S_q with cdf(S_q) = f * p_del
            x[i] = np.interp(target, cdf_neg, neg)

        return x

    @classmethod
    def _omega_components(cls, dfe: fd.DFE) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute omega (model dN/dS) and omega_a (its adaptive part) per bootstrap.

        Mirrors fastDFE's own alpha computation: y = phi(S) * q(S), where q(S) is the
        per-site fixation rate relative to neutral (S / (1 - exp(-S)) for h=0.5). Then
        omega = sum(y) and omega_a = sum(y[S > 0]), so alpha = omega_a / omega exactly.
        """
        disc = cls.discretization
        s = disc.s
        n = len(dfe.bootstraps)
        omega = np.zeros(n, dtype=float)
        omega_a = np.zeros(n, dtype=float)

        for i, bs in enumerate(dfe.get_bootstrap_dfes()):
            y = bs.model._discretize(bs.params, disc.bins) * disc.get_counts_fixed(bs.params.get('h', 0.5))
            omega[i] = y.sum()
            omega_a[i] = y[s > 0].sum()

        return omega, omega_a

    @classmethod
    def get_omega(cls, dfe: fd.DFE) -> np.ndarray:
        """
        Get omega (the model-predicted dN/dS) per bootstrap.
        """
        return cls._omega_components(dfe)[0]

    @classmethod
    def get_omega_a(cls, dfe: fd.DFE) -> np.ndarray:
        """
        Get omega_a (the adaptive component of dN/dS) per bootstrap.
        """
        return cls._omega_components(dfe)[1]
