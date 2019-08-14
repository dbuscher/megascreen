# Functions to calculate Zernike polynomials
#

import numpy as np
from numpy import sqrt
import functools
import sys


def NumZernike(m):
    """Return the number of polynomials up to and including radial order m"""
    return (m + 1) * (m + 2) // 2


def jtonm(j):
    """Convert Noll index j (1-based) to radial index n and azimuthal index m.  

    The method is described by V. N. Mahajan, Chapter 13, 
    Optical Shop Testing, 3rd. Ed., D. Malacara. ed., Wiley,2007"""
    n = int(sqrt(2 * j - 1) + 0.5) - 1
    m = (1 - 2 * (j % 2)) * (
        2 * int((2 * j + 1 + (n % 2) - n * (n + 1)) / 4.0) - (n % 2)
    )
    return (n, m)


@functools.lru_cache()
def ZernikeGrid(gridSize, maxRadial, diameter=None, orthoganalise=True):
    """Return Zernike polynomials sampled on a square grid.

    Returns all the polynomials for radial orders up to and including maxRadial.
    If diameter is None, then the unit circle just touches the outer edges of the 
    square (with the convention that the square extends +-0.5 pixels beyond the grid).
    Uses a recursive method to evaluate the radial polynomial.
    Returns the polynomials in the order given by Noll (1976), but zero-based.
    """
    if diameter == None:
        diameter = gridSize
    radius = diameter / 2.0
    y, x = np.mgrid[0:gridSize, 0:gridSize]
    x = (x - (gridSize - 1.0) / 2.0) / radius
    y = (y - (gridSize - 1.0) / 2.0) / radius
    # Convert to polar co-ordinates
    temp = x + 1j * y
    r = np.abs(temp)
    R = KinterRadial(r, maxRadial)
    # Create entheta[n] = exp(i*n*theta) recursively from exp(i*theta)
    eitheta = np.where(r == 0.0, 1.0, temp / r)
    entheta = np.concatenate(
        (
            np.ones((1,) + eitheta.shape, dtype=np.complex),
            np.cumprod(np.broadcast_to(eitheta, (maxRadial,) + eitheta.shape), axis=0),
        )
    )
    cntheta = entheta.real
    sntheta = entheta.imag
    jmax = NumZernike(maxRadial)
    Z = np.empty((jmax,) + r.shape)
    for j in range(jmax):
        n, m = jtonm(j + 1)
        const = sqrt((2 - (m == 0)) * (n + 1))
        if m < 0:
            Z[j] = R[n, -m] * sntheta[-m] * const
        else:
            Z[j] = R[n, m] * cntheta[m] * const
    # Make zernike zero outside unit circle (useful for dot product)
    Z = Z * np.less_equal(r, 1.0)
    if orthoganalise:
        return Orthoganalise(Z)
    else:
        return Z


def KinterRadial(r, maxRadial):
    """Compute Zernike radial functions using the modified Kintner method.

    Discussion may be found in Chong, Chee-Way, P. Raveendran, and R. Mukundan. 
    'A Comparative Analysis of Algorithms for Fast Computation of Zernike Moments'.
    Pattern Recognition 36, no. 3 (March 2003): 731-42. 
    doi:10.1016/S0031-3203(02)00091-2."""
    # Storage array for result
    # Creates array 4 * larger than necessary, but simpler to code
    R = np.empty((maxRadial + 1, maxRadial + 1) + r.shape)
    # Calc diagonal elements recursively
    R[0, 0] = np.ones_like(r)
    q = np.arange(1, maxRadial + 1)
    R[q, q] = np.cumprod(np.broadcast_to(r, (maxRadial,) + r.shape), axis=0)
    # Calc second batch of elements according to Fig 3 in Chong et al.
    for q in range(maxRadial - 1):
        R[q + 2, q] = (q + 2) * R[q + 2, q + 2] - (q + 1) * R[q, q]
    # Calc remaining elements
    r2 = r ** 2
    for q in range(maxRadial - 3):
        for p in range(q + 4, maxRadial + 1, 2):
            k1 = (p + q) * (p - q) * (p - 2) / 2.0
            k2 = 2 * p * (p - 1) * (p - 2)
            k3 = -q * q * (p - 1) - p * (p - 1) * (p - 2)
            k4 = -p * (p + q - 2) * (p - q - 2) / 2.0
            R[p, q] = ((k2 / k1) * r2 + (k3 / k1)) * R[p - 2, q] + (k4 / k1) * R[
                p - 4, q
            ]
    return R


def ZernikeHardCoded(gridSize, maxRadial, diameter=None):
    """
    Hardcoded Zernike routine: works only up to 5th radial order
    """
    assert maxRadial <= 5
    if diameter == None:
        diameter = gridSize
    radius = diameter / 2.0
    y, x = np.mgrid[0:gridSize, 0:gridSize]
    x = (x - (gridSize - 1.0) / 2.0) / radius
    y = (y - (gridSize - 1.0) / 2.0) / radius
    # Derive radius and exp(i*theta)
    temp = x + 1j * y
    r1 = np.abs(temp)
    eitheta = np.where(r1 == 0.0, 1.0, temp / r1)
    # Generate powers of r recursively
    r2 = r1 * r1
    r3 = r2 * r1
    r4 = r3 * r1
    r5 = r4 * r1
    # Generate cos and sin terms recursively from exp(i*theta)
    e2 = eitheta * eitheta
    e3 = e2 * eitheta
    e4 = e3 * eitheta
    e5 = e4 * eitheta
    ctheta = eitheta.real
    stheta = eitheta.imag
    c2theta = e2.real
    s2theta = e2.imag
    c3theta = e3.real
    s3theta = e3.imag
    c4theta = e4.real
    s4theta = e4.imag
    c5theta = e5.real
    s5theta = e5.imag
    # Generate all the zernikes
    zernike = np.zeros((21,) + x.shape)
    zernike[0] = 1.0
    zernike[1] = 2.0 * r1 * ctheta
    zernike[2] = 2.0 * r1 * stheta
    zernike[3] = sqrt(3.0) * (2.0 * r2 - 1.0)
    zernike[4] = sqrt(6.0) * r2 * s2theta
    zernike[5] = sqrt(6.0) * r2 * c2theta
    zernike[6] = sqrt(8.0) * (3.0 * r3 - 2.0 * r1) * stheta
    zernike[7] = sqrt(8.0) * (3.0 * r3 - 2.0 * r1) * ctheta
    zernike[8] = sqrt(8.0) * r3 * s3theta
    zernike[9] = sqrt(8.0) * r3 * c3theta
    zernike[10] = sqrt(5.0) * (6.0 * r4 - 6.0 * r2 + 1.0)
    zernike[11] = sqrt(10.0) * (4.0 * r4 - 3.0 * r2) * c2theta
    zernike[12] = sqrt(10.0) * (4.0 * r4 - 3.0 * r2) * s2theta
    zernike[13] = sqrt(10.0) * r4 * c4theta
    zernike[14] = sqrt(10.0) * r4 * s4theta
    zernike[15] = sqrt(12.0) * (10 * r5 - 12 * r3 + 3 * r1) * ctheta
    zernike[16] = sqrt(12.0) * (10 * r5 - 12 * r3 + 3 * r1) * stheta
    zernike[17] = sqrt(12.0) * (5 * r5 - 4 * r3) * c3theta
    zernike[18] = sqrt(12.0) * (5 * r5 - 4 * r3) * s3theta
    zernike[19] = sqrt(12.0) * r5 * c5theta
    zernike[20] = sqrt(12.0) * r5 * s5theta
    # Make zernike zero outside unit circle (useful for dot product)
    zernike = zernike * np.less_equal(r1, 1.0)
    return zernike[: NumZernike(maxRadial)]


def Orthoganalise(Zernikes):
    """Orthoganalise Zernikes to mitigate effects of under-sampling.
    """
    colZernikes = np.transpose(np.reshape(Zernikes, (len(Zernikes), -1)))
    q, r = np.linalg.qr(colZernikes, mode="reduced")
    orthoZernikes = np.reshape(np.transpose(q), Zernikes.shape)
    zernikeNorm = r[0,0]
    orthoZernikes *= zernikeNorm
    return orthoZernikes


def test_vs_hardcoded(gridSize=100, maxRadial=5):
    Z1 = ZernikeGrid(gridSize, maxRadial, orthoganalise=False)
    Z2 = ZernikeHardCoded(gridSize, maxRadial)
    diff = Z1 - Z2
    rms = np.sqrt(np.sum(diff ** 2, axis=(-1, -2)))
    print(rms)
    print(np.amax(rms))


