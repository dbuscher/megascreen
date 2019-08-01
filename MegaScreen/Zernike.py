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
    """Converts Noll index j to indices n and m.  The method is described by V. N. Mahajan, Chapter 13, Optical Shop Testing, 3rd. Ed., D. Malacara. ed., Wiley,2007"""
    n = int(sqrt(2 * j - 1) + 0.5) - 1
    if n % 2 == 0:
        m = 2 * int((2 * j + 1 - n * (n + 1)) / 4.0)
    else:
        m = 2 * int((2 * (j + 1) - n * (n + 1)) / 4.0) - 1
    if j % 2 != 0:
        m = -m
    return (n, m)

#@functools.lru_cache()
def ZernikeGrid(gridSize, maxRadial, diameter=None):
    if diameter == None:
        diameter = gridSize
    radius = diameter / 2.0
    x, y = np.mgrid[0:gridSize, 0:gridSize]
    x = (x - (gridSize - 1.0) / 2.0) / radius
    y = (y - (gridSize - 1.0) / 2.0) / radius
    return ZernikeKintner(x, y, maxRadial)

def ZernikeKintner(x, y, maxRadial):
    """Return Zernike values at a given set of x,y coords up to some maximum Noll index.  Uses modified Kintner method.  If unitNorm=True the polynomials will all have unit normalisation, otherwise they will have standard normalisation.  j is ordered first by n, with lower j corresponding to lower n, and then by m, with lower m coming first.  Discussion may be found in 
  'Chong, Chee-Way, P. Raveendran, and R. Mukundan. 'A Comparative Analysis of Algorithms for Fast Computation of Zernike Moments'.  Pattern Recognition 36, no. 3 (March 2003): 731-42. doi:10.1016/S0031-3203(02)00091-2.'"""

    # Set bounds on indices j, n, m, p, q
    jmax = NumZernike(maxRadial)
    nmmax = jtonm(jmax)
    nmax = nmmax[0]
    mmax = nmmax[1]
    pmax = nmax
    qmax = int(np.absolute(mmax))
    # Dummy variable for holding max value of p or q in specific case
    dmax = 0
    # Storage array for radial and Zernike functions
    # Creates array far larger than necessary, 200*200 for ~ 200 elements but numpy and linux employ lazy memory allocation so only the memory that is used is allocated.  This comes at a cost, if you use more memory than is available the program will crash at runtime.
    R = np.zeros((pmax + 1, pmax + 1) + x.shape)
    Z = np.zeros((jmax,) + x.shape)
    # Polar co-ordinates
    r = sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    # Calc radial part of Zernikes
    # Calc diagonal elements first
    if qmax == pmax - 2 or qmax == pmax:
        dmax = pmax + 1
    else:
        dmax = pmax
    for q in range(dmax):
        R[q, q] = np.power(r, q)
    # Calc second batch of elements according to diagram in source paper
    for q in range(dmax - 2):
        R[q + 2, q] = (q + 2) * R[q + 2, q + 2] - (q + 1) * R[q, q]
    # Calc remaining elements
    for q in range(pmax - 2):
        if q > qmax:
            dmax = pmax
        else:
            dmax = 1 + pmax
        for p in range(q + 4, dmax, 2):
            k1 = (p + q) * (p - q) * (p - 2) / 2.0
            k2 = 2 * p * (p - 1) * (p - 2)
            k3 = -q * q * (p - 1) - p * (p - 1) * (p - 2)
            k4 = -p * (p + q - 2) * (p - q - 2) / 2.0
            R[p, q] = ((k2 * r ** 2 + k3) * R[p - 2, q] + k4 * R[p - 4, q]) / k1

    # Fill Zernike array
    for j in range(0, jmax):
        nm = jtonm(j + 1)
        n = nm[0]
        m = nm[1]
        p = n
        q = int(np.absolute(m))
        if m < 0:
            Z[j] = R[p, q] * np.sin(q * phi)
        else:
            Z[j] = R[p, q] * np.cos(q * phi)
        # Normalise
        Z[j] *= sqrt((2 - (m == 0)) * (n + 1))
    # Make zernike zero outside unit circle (useful for dot product)
    Z = Z * np.less_equal(r, 1.0)
    return Z

def ZernikeHardCoded(gridSize, maxRadial, diameter=None):
    """
    Hardcoded Zernike routine: works only up to 5th radial order
    """
    assert maxRadial<=5
    if diameter == None:
        diameter = gridSize
    radius = diameter / 2.0
    x, y = np.mgrid[0:gridSize, 0:gridSize]
    x = (x - (gridSize - 1.0) / 2.0) / radius
    y = (y - (gridSize - 1.0) / 2.0) / radius
    # Derive radius and exp(i*theta)
    temp = x + 1j * y
    r1 = np.abs(temp)
    e1 = np.where(r1==0.0,1.0,temp / r1)
    # Generate powers of r recursively
    r2 = r1 * r1
    r3 = r2 * r1
    r4 = r3 * r1
    r5 = r4 * r1
    # Generate cos and sin terms recursively from exp(i*theta)
    e2 = e1 * e1
    e3 = e2 * e1
    e4 = e3 * e1
    e5 = e4 * e1
    ctheta = e1.real
    stheta = e1.imag
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


def test_vs_hardcoded(gridSize=100,maxRadial=5):
    Z1=ZernikeGrid(gridSize,maxRadial)
    Z2=ZernikeHardCoded(gridSize,maxRadial)
    diff = Z1 - Z2
    print(np.sqrt(np.sum(diff**2,axis=(-1,-2))))



def GramSchmidt(Zernikes, unitNorm=False):
    """Gram-Schmidt orthoganalise a set of vectors.

    Uses QR factorisation in numpy.
 
    unitNorm=True => ouput array will be orthonormal, otherwise the normalisation will be the same as for the input array"""

    numZernikes = len(Zernikes)
    diameter = len(Zernikes[0])
    if diameter ** 2 < numZernikes:
        raise ValueError(
            "Orthogonalisation only implemented for diameter^2 > numZernikes"
        )
    # Reshape Zernike array into column vectors
    colZernikes = np.reshape(
        np.ravel(Zernikes, order="C"), (diameter ** 2, numZernikes), order="F"
    )
    q, r = np.linalg.qr(colZernikes, mode="reduced")
    # Reshape column vectors back to array of meshgrids
    orthoZernikes = np.reshape(
        np.ravel(q, order="F"), (numZernikes, diameter, diameter), order="C"
    )
    if not unitNorm:
        for i in np.arange(numZernikes):
            zernikeNorm = sqrt(Sum2d(np.square(Zernikes[i])))
            orthoZernikes[i] = np.multiply(orthoZernikes[i], zernikeNorm)
    return orthoZernikes


if __name__ == "__main__":
    test()
