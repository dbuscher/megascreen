"""
Calculate theoretical predictions to compare against simulations
"""
import numpy as np
from numpy import pi
from scipy.special import gamma, jn, j1
import scipy.integrate
from MegaScreen import VonKarmanSpectrum


def winker_variance_quad(n, R=1, r0=2, L0=1e8):
    """Evaluate the variance of a Zernike of radial order n using the method of Winker 1991"""
    return (
        scipy.integrate.quad(
            lambda f: jn(n + 1, (2 * np.pi) * f) ** 2
            * VonKarmanSpectrum(f / R, r0, L0)
            / f,
            0.0,
            np.inf,
        )[0]
        * (n + 1)
        * (2 / np.pi)
        / R ** 2
        * (2 * R / r0) ** (5 / 3)
    )


def jinc_filter(f, fmin):
    w = 2 * pi * f
    return (1 - (4 * (j1(w) / (w)) ** 2)) if f > fmin else (w ** 2 / 4)


def winker_piston_residual(r0=1.0, L0=1e8, lower_limit=0, fmin=1e-4):
    return (
        scipy.integrate.quad(
            lambda f: (
                (1 - (4 * (j1(2 * pi * f) / (2 * pi * f)) ** 2))
                if f > fmin
                else (2 * pi * f) ** 2 / 4
            )
            * f
            * VonKarmanSpectrum(f, r0, L0),
            lower_limit,
            np.inf,
        )[0]
        * (2 * pi)
        * 2 ** (-5 / 3)
    )


def interf_spectrum_quad(baseline, frequencies, r0, L0, eps1=1.0):
    return np.array(
        [
            8
            * (
                scipy.integrate.quad(
                    lambda u, f, b, r0, L0: VonKarmanSpectrum(
                        np.sqrt(u ** 2 + f ** 2), r0, L0
                    )
                    * (1.0 - np.cos(2 * np.pi * u * b)),
                    0,
                    eps1 / baseline,
                    args=(f, baseline, r0, L0),
                    epsrel=1e-3,
                    limit=400,
                )[0]
                + scipy.integrate.quad(
                    lambda u, f, b, r0, L0: VonKarmanSpectrum(
                        np.sqrt(u ** 2 + f ** 2), r0, L0
                    ),
                    eps1 / baseline,
                    np.inf,
                    args=(f, baseline, r0, L0),
                    epsrel=1e-4,
                    limit=100,
                )[0]
                - scipy.integrate.quad(
                    lambda u, f, b, r0, L0: VonKarmanSpectrum(
                        np.sqrt(u ** 2 + f ** 2), r0, L0
                    ),
                    eps1 / baseline,
                    np.inf,
                    args=(f, baseline, r0, L0),
                    epsrel=1e-4,
                    weight="cos",
                    wvar=2 * np.pi * baseline,
                    limit=100,
                )[0]
            )
            for f in frequencies
        ]
    )


def sf_integrated(x, r0=1.0, L0=2e9):
    return [
        2
        * scipy.integrate.quad(
            lambda f, r, r0, L0: VonKarmanSpectrum(f, r0, L0)
            * (1.0 - scipy.special.jn(0, 2 * np.pi * f * r))
            * 2
            * pi
            * f,
            0,
            np.inf,
            args=(r, r0, L0),
            epsrel=1e-3,
        )[0]
        for r in x
    ]
