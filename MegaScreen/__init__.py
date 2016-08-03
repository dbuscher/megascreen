# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
from numpy import sqrt, fft, random, pi
import functools
import scipy.interpolate

__version__ = "0.4.1"


def FrequencyGrid(shape, pixelSize=1.0):
    """Return a 2-d grid with absolute frequency relevant to an FFT of a grid of pixels of size pixelSize."""
    return sqrt(np.add.outer(np.fft.fftfreq(shape[0], pixelSize) ** 2,
                             np.fft.fftfreq(shape[1], pixelSize) ** 2,
                             ))


def VonKarmanSpectrum(f, r0, L0=1e6, alpha=11. / 3.):
    """Phase spectrum of atmospheric seeing with Von Karman turbulence"""
    return 0.0229 * r0 ** (2.0 - alpha) * (f ** 2 + 1 / L0 ** 2) ** (-alpha / 2)


def FftScreen(spectrum, shape, pixelSize=1.0):
    """Generate infinite sequence of screens based on filtered 2D white noise

    Parameters
    ----------
    spectrum : function(f)
       spectrum of fluctuations, assumed radially symmetric
    shape: Sequence[int]
        Size of output grid
    pixelSize: float
        Pixel size of output grid

    Yields
    ------
    out : ndarray
        2D array of phase disturbances
    """
    f = FrequencyGrid(shape, pixelSize)
    filter = sqrt(spectrum(f) * f[0, 1] * f[1, 0])
    while 1:
        sample = random.normal(size=(2,) + filter.shape)
        result = fft.fft2(filter * (sample[0] + 1j * sample[1]))
        yield result.real
        yield result.imag


def McGlameryScreen(r0, L0=1e5, nfft=256):
    return FftScreen(functools.partial(VonKarmanSpectrum, r0=r0, L0=L0), shape=(nfft, nfft))


def SplineTiles(tileGenerator):
    """Generate a sequence of splined tiles with shape (n0/2,n1) from a
    sequence of tiles with shape (n0,n1)
    """
    previous = next(tileGenerator)
    n0 = previous.shape[0] // 2
    assert (n0 * 2 == previous.shape[0])
    cspline = np.cos(np.linspace(0, pi / 2, n0, endpoint=False))
    sspline = np.sin(np.linspace(0, pi / 2, n0, endpoint=False))
    for current in tileGenerator:
        yield previous[n0:] * cspline[:, np.newaxis] + current[:n0] * sspline[:, np.newaxis]
        previous = current


def GridInterpolator(grid):
    # print(grid.shape)
    xgrid = np.arange(grid.shape[0])
    ygrid = np.arange(grid.shape[1])
    return scipy.interpolate.RectBivariateSpline(xgrid, ygrid, grid)


def SlidingPixels(tileGenerator, x, y, dx):
    """Return phase values from a set of pixel coordinates sliding along an infinite ribbon.

    Parameters
    ----------
    tileGenerator: iterator
        A sequence of 2D phase screens which are stiched together to form a ribbon
    x,y: 1D arrays
        starting pixel coordinates in units of the tile grid size
    dx: float
        increment of pixel `x` coordinate on each iteration
    Yields
    -------
    1D array
        phase values at each pixel
    """
    tiles = [next(tileGenerator)]
    xtile = tiles[0].shape[0]
    assert xtile >= dx
    xmin = np.amin(x)
    xmax = np.amax(x)
    numTile = int(np.ceil((xmax - xmin + dx) / xtile))
    for i in range(numTile):
        tiles.append(next(tileGenerator))
    interpolator = GridInterpolator(np.concatenate(tiles))
    xoffset = -xmin
    ynew = y - np.amin(y)
    while True:
        yield interpolator(x + xoffset, ynew, grid=False)
        xoffset += dx
        if xoffset + xmin > xtile:
            tiles.pop(0)
            tiles.append(next(tileGenerator))
            interpolator = GridInterpolator(np.concatenate(tiles))
            xoffset -= xtile


def PixelCoords(origin, shape, pixelSize=1, theta=0):
    """Return x and y coodinates of a grid of pixels in rectangular region
    given by *origin* and *shape*, in a frame scaled to *pixelSize* and
    rotated by angle *theta*
    """
    c = np.cos(theta)
    s = np.sin(theta)
    print(origin, shape, pixelSize)
    x = (origin[0] + np.arange(shape[0])) * pixelSize
    y = (origin[1] + np.arange(shape[1])) * pixelSize
    return np.add.outer(c * x, s * y).flatten(), \
           np.add.outer(-s * x, c * y).flatten()


def SlidingWindows(tileGenerator, shape, dx, origins=((0.0, 0.0),), pixelSize=1, theta=0.0):
    """Return phase values from a set of rectangular windows sliding along an infinite ribbon.

    Parameters
    ----------
    tileGenerator: iterator
        A sequence of 2D phase screens which are stiched together to form a ribbon
    origins: sequence of pairs of floats
        The origins of each of the rectangular windows
    shape: tuple
        Shape of rectangular window (same for all windows)
    dx: float
        Increment of pixel `x` coordinate on each iteration
     Yields
     -------
     2D array
         phase values at each pixel
     """

    coords = [PixelCoords(origin=origin, shape=shape, pixelSize=pixelSize, theta=theta)
              for origin in origins]
    # print(coords)
    coords = np.array(coords)
    x = coords[:, 0, :].flat
    y = coords[:, 1, :].flat
    numWindow = len(origins)
    if numWindow == 1:
        newshape = shape
    else:
        newshape = [numWindow] + list(shape)
    for screen in SlidingPixels(tileGenerator, x, y, dx):
        yield np.reshape(screen, newshape)


def NestedSpectra(spectrum, f0, eps=1e-6):
    grad = (spectrum(f0 * (1 + eps)) - spectrum(f0 * (1 - eps))) / (2 * f0 * eps)
    c1 = spectrum(f0)

    def OuterSpectrum(f):
        s = spectrum(f)
        s1 = np.where(f < f0, c1 - 2 * grad * f0 / np.pi * np.cos(np.pi * f / (2 * f0)), s)
        return np.where(s1 < s, s1, s)

    def InnerSpectrum(f):
        return spectrum(f) - OuterSpectrum(f)

    return InnerSpectrum, OuterSpectrum


def NestedScreen(spectrum, windowShape, dx,
                 windowOrigins=((0.0, 0.0),), pixelSize=1.0, theta=0.0,
                 nfftWoofer=256, nfftTweeter=256, frequencyOverlap=4.0, fractionalSupport=0.5,
                 debug=False):
    """Generate a sequence of phase screens using nested FFTs

    Parameters
    ----------
    spectrum: function
       Returns the spectral power of the phase perturbations at a given frequency
    windowShape: pair of ints
       Shape in pixels of rectangular windows
    windowOrigins: sequence of pairs of floats
        The origins of each of the rectangular windows
    pixelSize : float
       The size of the output pixels in terms of the pixel size of
       the tweeter screen
    theta: float
        Rotation of the window coordinate system wrt the x axis, to allow different
        'wind' directions
    dx : float
       Increment in the x coordinate of all windows between subsequent calls,
       in units of the outer spectrum grid size
    nfftWoofer, nfftTweeter: float
        Size of the square FFTs used to produce the woofer and tweeter
        screens
    frequencyOverlap : float
       The Nyquist frequency of the woofer spectrum in units of the fundamental
       frequency of the tweeter spectrum
    fractionalSupport : float
        Fraction of Nyquist frequency over which woofer spectrum is non-zero
    """
    wooferPixelSize = nfftTweeter / (2 * frequencyOverlap)
    f0 = 1 / (2 * wooferPixelSize) * fractionalSupport
    wooferSpectrum, tweeterSpectrum = NestedSpectra(spectrum, f0)
    innerWindows = SlidingWindows(SplineTiles(FftScreen(wooferSpectrum, (nfftWoofer, nfftWoofer), wooferPixelSize)),
                                  dx=dx / wooferPixelSize,
                                  shape=windowShape, origins=windowOrigins,
                                  pixelSize=pixelSize / wooferPixelSize, theta=theta)
    outerWindows = [SlidingWindows(SplineTiles(FftScreen(tweeterSpectrum, (nfftTweeter, nfftTweeter), pixelSize)), dx=dx,
                                   shape=windowShape, origins=[origin],
                                   pixelSize=pixelSize, theta=theta)
                    for origin in windowOrigins]
    while 1:
        inner = next(innerWindows)
        outer = np.squeeze(np.array([next(o) for o in outerWindows]))
        if debug:
            yield inner, outer, inner + outer
        else:
            yield inner + outer


def MegaScreen(r0, L0, windowShape, dx,
               windowOrigins=((0.0, 0.0),), pixelSize=1.0, theta=0.0,
               nfftWoofer=256, nfftTweeter=256,
               frequencyOverlap=4.0, fractionalSupport=1.0,
               debug=False):
    """

    Parameters
    ----------
    r0 : Fried parameter for phase screen
    L0 : Outer scale of turbulence

    Other Parameters and Return value: as for NestedScreen()

    """
    spectrum = functools.partial(VonKarmanSpectrum, r0=r0, L0=L0)
    return NestedScreen(spectrum, windowShape, dx,
                        windowOrigins=windowOrigins, pixelSize=pixelSize, theta=theta,
                        nfftWoofer=nfftWoofer, nfftTweeter=nfftTweeter,
                        frequencyOverlap=frequencyOverlap, fractionalSupport=fractionalSupport,
                        debug=debug)
