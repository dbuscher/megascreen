import numpy as np
import scipy.integrate
import MegaScreen as ms
import matplotlib.pyplot as plt
import sys
from astropy.table import Table
import time


def SaveTable(results,prefix,args):
    for key in sorted(args):
        results.meta[key]=args[key]
    results.write(time.strftime("data/"+prefix+"%y%m%d-%H%M.dat"),
                format='ascii.ecsv')

def CookieCutterSubScreens(screenGenerator, nx, ny):
    """Generate a a sequence of rectangular screens cut
    out from a larger screen"""
    for motherScreen in screenGenerator:
        for iy in range(0, motherScreen.shape[0] - ny + 1, ny):
            for ix in range(0, motherScreen.shape[1] - nx + 1, nx):
                yield motherScreen[iy:iy + ny, ix:ix + nx]


def structure_function_brute_force(sig):
    l = len(sig)
    a = [np.sum((sig[:-lag] - sig[lag:]) ** 2, axis=0) / (l - lag) for lag in range(1, l)]
    return np.array(a)


def average_sf(sig):
    a = structure_function_brute_force(sig)
    return np.mean(a, axis=-1)


def average_sf_xy(screen, decimate):
    return average_sf(screen[:, ::decimate]), average_sf(screen.transpose()[:, ::decimate])


def screen_sf_xy(generator, decimate, numScreen):
    sf = []
    for i in range(numScreen):
        screen = next(generator)
        sf.append(average_sf_xy(screen, decimate))
    return np.mean(sf, axis=0)


def mcglamery_sf(r0=5, L0=10000, diameter=200, nfft=1024, decimate=20, numScreen=4000):
    args = locals()
    args["MegaScreenVersion"] = ms.__version__
    generator = CookieCutterSubScreens(ms.McGlameryScreen(r0, L0, nfft), diameter, diameter)
    sf_x, sf_y = screen_sf_xy(generator, decimate, numScreen)
    r = np.arange(len(sf_x)) + 1
    model = sf_integrated(r, r0, L0)
    results = Table([r, model, sf_x, sf_y], names=["r", "model", "sf_x", "sf_y"])
    SaveTable(results,prefix,args)
    return results


def multiscreen_sf(r0=5, L0=10000, diameter=128,
                   decimate=20, numScreen=10000,
                   frequencyOverlap=4.0, fractionalSupport=1.0,
                   nfftWoofer=512, nfftTweeter=256,
                   prefix="multi"):
    args = locals()
    args["MegaScreenVersion"] = ms.__version__
    generator = ms.MegaScreen(r0, L0, [diameter, diameter], dx=diameter,
                              frequencyOverlap=frequencyOverlap, fractionalSupport=fractionalSupport,
                              nfftWoofer=nfftWoofer, nfftTweeter=nfftTweeter,
                              debug=True)
    sf = []
    for i in range(numScreen):
        screens = next(generator)
        sf.append([average_sf_xy(screens[i], decimate) for i in range(3)])
    results = np.mean(sf, axis=0)
    r = np.arange(results.shape[2]) + 1
    model = sf_integrated(r, r0=r0, L0=L0)
    t = Table({"r": r, "model": model})
    for i in range(3):
        sf_x, sf_y = results[i]
        t["sf_x" + str(i)] = sf_x
        t["sf_y" + str(i)] = sf_y
    SaveTable(t,prefix,args)
    return t


def sf_integrated(x, r0=1.0, L0=2e9):
    return [2 * scipy.integrate.quad(lambda f, r, r0, L0:
                                     VonKarmanSpectrum(f, r0, L0) *
                                     (1.0 - scipy.special.jn(0, 2 * np.pi * f * r)) * 2 * pi * f,
                                     0, np.inf, args=(r, r0, L0), epsrel=1e-3)[0]
            for r in x]



# Generic command-line interface to run test function or give full control
if __name__ == '__main__':
    if len(sys.argv) == 1:
        multiscreen_sf()
    else:
        exec(sys.argv[1])
