========================
 Basic parameter tuning
========================

Spatial sampling
^^^^^^^^^^^^^^^^
The `MegaScreen` generator yields one or more rectangular "windows" on to an infinite turbulent phase screen. The simplest case is a single window, useful for simulating a conventional single-aperture observation. The data represents the wavefront perturbations across the window in radians, sampled on a grid of square "pixels".

The physical dimensions of this grid are affected by the physical scale of the perturbations being simulated and by three parameters: `r0`, `windowShape`, and `pixelSize`. By default, these parameters have the values 7.0, `(100,100)` and `1.0` respectively, and this corresponds to sampling the perturbations at 7 samples per Fried parameter (:math:`r_0`) inside a square window of 100x100 samples. A circular telescope aperture just filling this window would have a diameter :math:`D` corresponding to :math:`D/r_0=100/7\approx 14`. At optical wavelengths, where typical values of  :math:`r_0` might be around 10cm, this corresponds to a 1.4m aperture diameter, whereas at
an infrared wavelength with :math:`r_0=60\,cm` this would correspond to a physical aperture size of 8.6m.


The `r0` parameter in the call to `MegaScreen()` sets the "tweeter" screen sampling frequency in terms of samples per :math:`r_0` in both the "x" and "y" coordinates and this affects the accuracy of the simulated fluctuations. Setting `r0=7.0` means that the mean squared error due to not simulating wavefront "corrugations" at spatial frequencies above the Nyquist frequency of the tweeter FFT is about
0.01 radians squared. Setting `r0` to values of at least 5 should give adequate accuracy for most applications. Higher accuracy may be needed for simulating ultra-high-contrast experiments using extreme AO, where the wavefront residuals are much smaller than a radian. In such cases significantly larger values of the `r0` parameter may be required.

The value of `pixelSize` is
the size of the interpolated window pixels in units of the tweeter
screen pixel size, and this should usually be less than or equal to unity, otherwise there will be aliasing of the high-spatial-frequency wavefront perturbations to lower frequencies. The value you
choose is application-dependent, but leaving this value at unity is usually acceptable.

When simulating 
a circular aperture, it is best to choose the combination of `r0` and `pixelSize` so that the aperture is sampled with at least  30 samples across a diameter. Qualitatively, using aperture sampling of this order leads to 
the pixelated version of the circular aperture looking reasonably
circular, rather than like a square with the corners cut off, and the diffraction pattern shows low levels of square as opposed to circular symmetry. Quantitatively, good sampling of the
circular pupil is necessary so that the Zernike polynomials up to some desired order are orthogonal to one another to acceptable accuracy. When in doubt, it pays to check what the effects are on the results from any given simulation of, for example, doubling the number of pixels across a given aperture

For computational efficiency, the value of `windowShape` should typically be chosen to just enclose the required aperture.
The dimensions of the window must fit into a single tweeter "strip" which is an infinite strip `nfftTweeter` pixels wide. In fitting into this restriction, `windowShape`, `pixelSize` and the rotation of the window set by `theta` are all important, since the diagonal of the rectangular window is larger than either of the edge dimensions. If a larger window is required, `nfftTweeter` should be increased. Guidance on choosing `nfftTweeter` is given in a later section. 




Temporal sampling
^^^^^^^^^^^^^^^^^^

To simulate temporal evolution of the wavefront, a "frozen turbulence" hypothesis is assumed, namely that the wind is blowing a non-evolving turbulent screen across the aperture. Each call to the generator yields a new snapshot of this screen, with the screen having translated by a distance `dx` between snapshots, where `dx` is expressed in tweeter pixel units.

This can be converted into a time interval via the windspeed, for example if the the parameter `r0` is 7 corresponding to a value of :math:`r_0` of 60cm, then each tweeter pixel is 8.6cm across. Assuming a windspeed of 10m/s then the default value of 3.5 for `dx` corresponds to a snapshot interval of approximately 30 milliseconds. Alternatively, using the conventional definition of the temporal coherence time of the turbulence :math:`t_0\equiv 0.314 r_0/V` where :math:`V` is the windspeed, then the interval between snapshots is given by :math:`3.185(dx/r0) t_0` or approximately :math:`1.6t_0` in this example.

This value for `dx` corresponds to a typical exposure time for a fringe pattern in an interferometer or a speckle pattern in a single telescope, but to simulate the "smearing" of a speckle or fringe pattern during an exposure time, or simulating temporal effects in an adaptive optics system, values of `dx` of order 10 times smaller may be required.

Setting the value of `dx` to be larger than the dimension of the window gives a set of relatively uncorrelated samples of the phase screen, which can be an efficient way to determine average properties of turbulence-affected quantities, e.g. mean image Strehl ratio, mean squared Zernike coefficients etc.

