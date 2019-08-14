"""
Statistical tests of Zernike decompositions of phase screens
Includes tests against the Noll (1976) results for the residual variance
of a wavefront with various orders of Zernike removal. 
"""
import numpy as np
from MegaScreen import *
from MegaScreen.Zernike import ZernikeGrid
import MegaScreen as ms
from astropy.table import Table
import time

Noll = [
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
        MegaScreen(
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
                MegaScreen(
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
    results.meta["MegaScreenVersion"] = ms.__version__
    results.write(
        time.strftime("tmp/" + prefix + "%y%m%d-%H%M.dat"), format="ascii.ecsv"
    )


if __name__ == "__main__":
    Winker(numIter=10000)
