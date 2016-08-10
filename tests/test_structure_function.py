import numpy as np
from MegaScreen import *
import MegaScreen as ms
from scipy import signal, special
import matplotlib.pyplot as plt
import sys
from astropy.io import ascii
from astropy.table import Table
import time
import inspect

def CookieCutterSubScreens(screenGenerator, nx, ny):
    """Generate a a sequence of rectangular screens cut
    out from a larger screen"""
    for motherScreen in screenGenerator:
        for iy in range(0,motherScreen.shape[0]-ny+1,ny):
            for ix in range(0,motherScreen.shape[1]-nx+1,nx):
                yield motherScreen[iy:iy+ny,ix:ix+nx]

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


def multiscreen_sf(r0=5, L0=10000, diameter=128,
                   decimate=20, numScreen=10000,
                   frequencyOverlap=4.0, fractionalSupport=1.0,
                   nfftWoofer=512,nfftTweeter=256,
                   prefix="multi"):
    args=locals()
    args["MegaScreenVersion"] = ms.__version__
    generator=MegaScreen(r0, L0, [diameter, diameter], dx=diameter,
                         frequencyOverlap=frequencyOverlap,fractionalSupport=fractionalSupport,
                         nfftWoofer=nfftWoofer,nfftTweeter=nfftTweeter,
                         debug=True)
    sf = []
    for i in range(numScreen):
        screens = next(generator)
        sf.append([average_sf_xy(screens[i], decimate) for i in range(3)])
    results = np.mean(sf, axis=0)
    r=np.arange(results.shape[2]) + 1
    model = sf_integrated(r, r0=r0, L0=L0)
    t=Table({"r": r, "model": model})
    for i in range(3):
        sf_x, sf_y = results[i]
        t["sf_x"+str(i)]=sf_x
        t["sf_y"+str(i)]=sf_y
    for key in sorted(args):
        t.meta[key]=args[key]
    t.write(time.strftime("tmp/"+prefix+"%y%m%d-%H%M.dat"),
            format='ascii.ecsv')
    plt.figure(figsize=(12, 4))
    for i in range(3):
        plt.subplot(1, 3, i + 1)
        sf_x, sf_y = results[i]
        plt.loglog(r, sf_x)
        plt.loglog(r, sf_y)
        plt.loglog(r, model)
    plt.show()


def mcglamery_sf(r0=5, L0=10000, diameter=200, nfft=1024, decimate=20, numScreen=4000):
    args=locals()
    args["MegaScreenVersion"] = ms.__version__
    generator=CookieCutterSubScreens(McGlameryScreen(r0, L0, nfft),diameter,diameter)
    t=sf_test(generator, r0=r0, L0=L0,
              decimate=decimate, numScreen=numScreen,args=args,prefix="mcglamery")
    plot_table(t)

def test2(r0=10, L0=2000, diameter=128, decimate=10, numScreen=1000):
    sf_test(SingleGenerator(r0, L0, diameter, dx=diameter, which=1), decimate=decimate, numScreen=numScreen)


def test_mega(r0=10, L0=10000, diameter=128, decimate=10, numScreen=100):
    args=locals()
    sf_test(MegaScreen(r0, L0, [diameter, diameter], dx=diameter,frequencyOverlap=2.0),
            r0=r0,L0=L0,
            decimate=decimate, numScreen=numScreen,
            args=args,prefix="mega")





def sf_integrated_single(r, r0, L0):
    return 2 * scipy.integrate.quad(lambda f, r, r0, L0:
                                    VonKarmanSpectrum(f, r0, L0) *
                                    (1.0 - scipy.special.jn(0, 2 * np.pi * f * r)) * 2 * pi * f,
                                    0, np.inf, args=(r, r0, L0), epsrel=1e-3)[0]



def sf_integrated(x, r0=1.0, L0=2e9):
    return [sf_integrated_single(r, r0, L0) for r in x]

def plot_integral(r0=1.0, L0=2e9):
    x = np.logspace(-1, 2)
    y = sf_integrated(x, r0, L0)
    plt.loglog(x, y)
    plt.loglog(x, 6.88 * (x / r0) ** (5. / 3.))
    plt.show()

def SingleGenerator(r0, L0, diameter, dx, which=0):
    for screens in MegaScreen(r0, L0, windowShape=[diameter, diameter], dx=dx, debug=True):
        yield screens[which]


def sf_test(generator, r0, L0, decimate=10, numScreen=4000,args={},prefix="tmp"):
    sf_x, sf_y = screen_sf_xy(generator, decimate, numScreen)
    r=np.arange(len(sf_x))+1
    model=sf_integrated(r,r0,L0)
    results = Table([r, model, sf_x, sf_y], names=["r","model", "sf_x", "sf_y"])
    for key in sorted(args):
        results.meta[key]=args[key]
    results.write(time.strftime("tmp/"+prefix+"%y%m%d-%H%M.dat"),
                format='ascii.ecsv')
    return results

def plot_table(t):
    r=t["r"]
    sf_x=t["sf_x"]
    sf_y=t["sf_y"]
    model=t["model"]
    plt.loglog(r,sf_x)
    plt.loglog(r,sf_y)
    plt.loglog(r,model)
    plt.show()


# Generic command-line interface to run test function or give full control
if __name__ == '__main__':
    if len(sys.argv) == 1:
        multiscreen_sf()
    else:
        exec(sys.argv[1])
