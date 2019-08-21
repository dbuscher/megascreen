"""
Statistical tests of Zernike decompositions of phase screens
Includes tests against the Noll (1976) results for the residual variance
of a wavefront with various orders of Zernike removal. 
"""
import numpy as np
from numpy import pi
import MegaScreen
from Zernike import ZernikeGrid, jtonm
from astropy.table import Table
import time

from scipy.special import gamma, jn, j1
import scipy.integrate


def noll_covariance_analytic(ni, nj, m=None, const=0.0229):
    """Evaluate covariance of Zernikes using equation A2 from Noll 1976

    Returns covariance of radial order ni and nj  and azimuthal
    order m.  Note that if m_j != m_i then the correlation is zero. 
    Does not check if m is allowed."""
    if m == None:
        m = ni % 2  # Lowest valid m
    return (
        const
        * 2 ** (-5 / 3)
        * (-1) ** ((ni + nj - 2 * m) / 2)
        * np.sqrt((ni + 1) * (nj + 1))
        * pi ** (8.0 / 3)
        * gamma(14.0 / 3.0)
        * gamma((ni + nj - 5.0 / 3) / 2)
        / (
            gamma((ni - nj + 17 / 3.0) / 2)
            * gamma((nj - ni + 17.0 / 3) / 2)
            * gamma((ni + nj + 23 / 3.0) / 2)
        )
    )


# Table IV from Noll 1976
noll_table_iv = [
    1.0299,
    0.582,
    0.134,
    0.111,
    0.0880,
    0.0648,
    0.0587,
    0.0525,
    0.0463,
    0.0401,
    0.0377,
    0.0352,
    0.0328,
    0.0304,
    0.0279,
    0.0267,
    0.0255,
    0.0243,
    0.0232,
    0.0220,
    0.0208,
]


def winker_variance_quad(n, R=1, r0=2, L0=1e8):
    """Evaluate the variance of a Zernike of radial order n using the method of Winker 1991"""
    return (
        scipy.integrate.quad(
            lambda f: jn(n + 1, (2 * np.pi) * f) ** 2
            * MegaScreen.VonKarmanSpectrum(f / R, r0, L0)
            / f,
            0.0,
            np.inf,
        )[0]
        * (n + 1)
        * (2 / np.pi)
        / R ** 2
        * (2 * R / r0) ** (5 / 3)
    )


def noll_variance_quad(n):
    return (
        scipy.integrate.quad(
            lambda f: jn(n + 1, (2 * np.pi) * f) ** 2 * f ** (-8 / 3) / f ** 2,
            0.0,
            np.inf,
        )[0]
        * (n + 1)
        * (0.0229 * 2 / np.pi)
        * 2 ** (-5 / 3)
    )


def noll_piston_residual_quad():
    return (
        scipy.integrate.quad(
            lambda f: (1 - (4 * (j1(2 * pi * f) / (2 * pi * f)) ** 2)) * f ** (-8 / 3),
            1e-8,
            np.inf,
        )[0]
        * (0.0229 * 2 * pi)
        * 2 ** (-5 / 3)
    )


def noll_piston_residual_analytic(p=8 / 3, const=0.0229):
    """Evaluates the analytic solution for the piston-removed
    residual atmospheric aberrations using the last equation in 
    the Appendix of Noll 1976

    This should give :math:`\Delta_1` from Table IV of that paper
    but the values deviate by more than 1 significant figure.
    """
    return (
        pi
        * gamma(p + 2)
        / (
            2 ** p
            * (gamma((p + 3) / 2)) ** 2
            * gamma((p + 5) / 2)
            * gamma((1 + p) / 2)
            * np.sin((pi / 2) * (p - 1))
        )
        * (2 * pi) ** (5 / 3)
        * (const * 2 * pi)
        * 2 ** (-5 / 3)
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
            * MegaScreen.VonKarmanSpectrum(f, r0, L0),
            lower_limit,
            np.inf,
        )[0]
        * (2 * pi)
        * 2 ** (-5 / 3)
    )


def CookieCutterSubScreens(screenGenerator, nx, ny):
    """Generate a a sequence of rectangular screens cut
    out from a larger screen"""
    for motherScreen in screenGenerator:
        for iy in range(0, motherScreen.shape[0] - ny + 1, ny):
            for ix in range(0, motherScreen.shape[1] - nx + 1, nx):
                yield motherScreen[iy : iy + ny, ix : ix + nx]


def NollVariance(screen, zernike):
    """ Return the values for the variance of a phase screen as increasing
    numbers of Zernikes are removed"""
    screen = screen * zernike[0]  # Cut out a circle
    normalisation = np.sum(zernike[0])
    result = np.zeros(len(zernike))
    for i in range(len(zernike)):
        screen = screen - zernike[i] * np.sum(zernike[i] * screen) / normalisation
        result[i] = np.sum(screen ** 2) / normalisation
    return result


def test1(diameter=32, screenSize=256, numIter=100):
    NollPrint(
        CookieCutterSubScreens(
            McGlameryScreen(nfft=screenSize, r0=float(diameter)), diameter, diameter
        ),
        diameter,
        numIter,
    )


def test(diameter=32, L0=2000, numIter=100):
    NollPrint(
        MegaScreen.MegaScreen(
            r0=float(diameter), L0=L0, windowShape=[diameter, diameter], dx=diameter
        ),
        diameter,
        numIter,
    )


def Winker(
    diameter=32,
    L0Min=16,
    L0Max=8000,
    numL0=20,
    numIter=100,
    maxRadial=2,
    nfftOuter=256,
    nfftInner=256,
    randomSeed=12345,
):
    """Derive simulated data corresponding to Figure 2 of Winker 1991"""

    L0s = np.logspace(np.log10(L0Min), np.log10(L0Max), numL0)
    variance = np.array(
        [
            NollTest(
                MegaScreen.MegaScreen(
                    r0=float(diameter),
                    L0=L0,
                    nfftWoofer=nfftInner,
                    nfftTweeter=nfftOuter,
                    windowShape=[diameter, diameter],
                    dx=diameter,
                ),
                diameter,
                numIter,
                maxRadial=maxRadial,
            )
            for L0 in L0s
        ]
    ).transpose()
    t = Table(
        [L0s] + list(variance),
        names=["L0"] + ["Z" + str(i) for i in range(len(variance))],
    )
    return t


def NollTest(generator, diameter, numIter, maxRadial=5):
    """Test Noll theory on simulated turbulence"""
    zernike = ZernikeGrid(diameter, maxRadial=maxRadial)
    numScreen = 0
    variance = np.zeros(len(zernike))
    for i in range(numIter):
        screen = next(generator)
        variance = variance + NollVariance(screen, zernike)
    return variance / numIter


def NollPrint(generator, diameter, numIter):
    variance = NollTest(generator, diameter, numIter)
    print("Zernikes removed     Residual variance    Noll Result")
    for i in range(len(variance)):
        print("%16d %16f %16f" % (i, variance[i], Noll[i]))


def SaveTable(results, prefix, args):
    for key in sorted(args):
        results.meta[key] = args[key]
    results.meta["MegaScreenVersion"] = MegaScreen.__version__
    results.write(
        time.strftime("tmp/" + prefix + "%y%m%d-%H%M.dat"), format="ascii.ecsv"
    )


if __name__ == "__main__":
    Winker(numIter=10000)
