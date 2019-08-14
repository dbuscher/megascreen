.. MegaScreen documentation master file, created by
   sphinx-quickstart on Fri Aug  9 15:47:48 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MegaScreen's documentation!
======================================

MegaScreen implements the algorithm described in
[Buscher2016]_ for simulating perturbations to an optical wavefront corresponding to the Von-Karman-Kolmogrov-Tatarski model for propagation through turbulence. From the abstract:

   We describe an algorithm for simulating atmospheric wavefront perturbations over ranges of spatial and temporal scales spanning more than 4 orders of magnitude.  An open-source implementation of the algorithm written in Python can simulate the evolution of the perturbations more than an order-of-magnitude faster than real time. Testing of the implementation using metrics appropriate to adaptive optics systems and long-baseline interferometers show accuracies at the few percent level or better.

Installation
^^^^^^^^^^^^

On unix-like systems do

::

    pip3 install MegaScreen

Note that we use `pip3` because this is a Python3 package. Alternative installation methods are given on the `github site`_.
    
Simple usage
^^^^^^^^^^^^


The main interface to the package is the MegaScreen function, which is a Python generator. The following code will generate 10000 snapshots of "frozen turbulence" blowing past a single rectangular window, and process them using a user-defined function `process()`:

.. code:: python

    import MegaScreen
    
    for phaseScreen in MegaScreen.MegaScreen(numIter=10000):
	process(phaseScreen)



More complex example code using this function is given in the `tests` directory on the `github site`_.

The function has a large number of arguments as described in the docstring below. All the parameters have sensible defaults, but these can be changed for different use cases as described in the sections below.

.. autofunction:: MegaScreen.MegaScreen
		  

Basic parameters
^^^^^^^^^^^^^^^^
The generator yields one or more rectangular "windows" on to an infinite turbulent phase screen. The simplest case is a single window, corresponding to a conventional single-aperture observation. The data represents the wavefront perturbations across the window in radians, sampled on a square grid of "pixels".

The physical dimensions of this grid is affected by three parameters: `r0`, `windowShape`, and `pixelSize`. By default, these have the values 7.0, `(100,100)` and `1.0` respectively, and this corresponds to sampling the perturbations at 7 samples per Fried parameter (:math:`r_0`) inside a square window of 100x100 samples. A circular telescope aperture just filling this window would have a diameter :math:`D` corresponding to :math:`D/r_0=100/7\approx 14`. At an infrared wavelength with :math:`r_0=60\,cm` this would correspond to a physical aperture size of :math:`D=100/7\times60\approx 8.6` m, whereas at optical wavelengths where typical values of  :math:`r_0` might be around 10cm, then this corresponds to a 1.4m aperture diameter. 

The value of `r0` in the call to MegaScreen
is the value of the number of "tweeter" screen samples per :math:`r_0` distance. It is recommended to use values greater than or equal to 5. For r0=7
then the missing power in the simulated perturbations due to not simulating wavefront "corrugations" at spatial frequencies above the Nyquist frequency of the tweeter FFT is about
0.01 radians squared, which should give adequate accuracy for most applications. Higher accuracy may be needed for simulating ultra-high-contrast experiments using extreme AO, in which case larger values of `r0` may be needed.

The value of `pixelSize` is
the size of the interpolated window pixels in units of the tweeter
screen pixel size, and this should usually be less than or equal to unity, otherwise there will be aliasing of the high-spatial-frequency wavefront perturbations to lower frequencies. The value you
choose is application-dependent, but leaving this value at unity is usually acceptable.

When sampling
a circular aperture, it is best to choose the combination of `r_0` and `pixelSize` to make the aperture larger than about  30
pixels across. Qualitatively, using aperture sampling of this order leads to 
the pixelated version of the circular aperture looking reasonably
circular, rather than like a square with the corners cut off, and the diffraction pattern shows low levels of square as opposed to circular symmetry. Quantitatively, good sampling of the
circular pupil is necessary so that the Zernike polynomials up to some desired order are orthogonal to one another to acceptable accuracy. When in doubt, it pays to check what the effects of, for example, doubling the number of pixels and halving the pixel size are on the results from any given simulation.

For computational efficiency, the value of `windowShape` should typically be chosen to just enclose the required aperture. One restriction is that the window should fit into a single tweeter "strip" which is an infinite strip `nfftTweeter` pixels wide. In fitting into this restriction, `windowShape`, `pixelSize` and the rotation of the window set by `theta` are important: remember, the diagonal of the rectangular window is larger than either of the edge dimensions. If a larger window is required, `nfftTweeter` should be increased as discussed below. 




Temporal evolution
^^^^^^^^^^^^^^^^^^

To simulate temporal evolution of the wavefront, a "frozen turbulence" hypothesis is assumed, namely that the wind is blowing a fixed turbulent screen across the aperture. Each call to the `MegaScreen` iterator yields a new snapshot of this screen, with the screen having translated by a distance `dx` between snapshots, where `dx` is expressed in tweeter pixel units. This is can be converted into a time interval via the windspeed, for example if the the parameter `r0` is 7 corresponding to a value of :math:`r_0` of 60cm, then each tweeter pixel is 8.6cm across. Assuming a windspeed of 10m/s then the default value of 3.5 for `dx` corresponds to a snapshot interval of :math:`3.5\times 8.6/1000\approx 30` milliseconds. Alternatively, using the conventional definition of the temporal coherence time of the turbulence :math:`t_0\equiv 0.314 r_0/V` where :math:`V` is the windspeed, then the interval between snapshots is given by :math:`3.185(dx/r0) t_0` or approximately :math:`1.6t_0` in this example.

This value for `dx` corresponds approximately to an exposure time for a speckle imaging system, but to simulate the "smearing" of a speckle or fringe pattern during an exposure time, or simulating temporal effects in an adaptive optics system, considerably smaller values of `dx` may be required. Setting the value of `dx` to be larger than the dimension of the window gives a set of relatively uncorrelated samples of the phase screen, which can be an efficient way to determine average properties of turbulence-affected quantities, e.g. mean squared image size, mean squared Zernike coefficients etc.

Simulating multiple turbulent layers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For simple applications, a single phase screen is sufficient. To add additional realism to the space-time correlation of the wavefonts, multiple screens with different wind speeds and directions are required to simulate the effect of multiple turbulent layers at different heights.
The
experimental evidence indicates that at least 3 different screens of approximately the same turbulent strength are needed to give the right level of decorrelation.

Each layer should be generated from an independent `MegaScreen` generator. If the field of view being considered is small enough that anisoplanatism effects are
negligible, then the phase screens so generated can simply be added together. If not, the phase screens need to be shifted with respect to one another by amounts dependending on the heights of the corresponding turbulent layers and the location within the field of view. Where scintillation effects are important, it may be necessary to simulate Fresnel propagation between different screens.

Each layer should have its own value for `theta`, `dx`, `r0` and possibly `L0`. The parameter `theta` sets the wind direction, where it should be noted that, as explained in the docstring, this is at 90 degrees to what might at first be expected, due to a FORTRAN-like convention for the array indexing. This should not usually be a problem, since what matters most is the difference in wind directions.

To derive appropriate `r_0` values for the different layers, we need to consider the effect of multiple turbulent layers on the wavefront propagation. The
variance of a given point on each phase screen is proportional to an integral of the turbulent strength
:math:`C_n^2` over some (usually thin) turbulent layer. By summing together several
statistically independent phase screens the variances add and so we get the
integral over the whole optical path. The variance is  proportional to :math:`r_0^{-5/3}` so
to get a given :math:`C_n^2` we set :math:`r_0` for each layer proportional to :math:`[\int  C_n^2 \,dz]^{-3/5}` for that layer. The effective  value of :math:`r_0` for the sum of all the screens is given by:

.. math::

   r_0 = \left[\sum_i (r_0)_i^{-5/3}\right]^{-3/5}

where :math:`(r_0)_i` is the Fried parameter of the :math:`i^{\text{th}}` layer.

Using multiple windows
^^^^^^^^^^^^^^^^^^^^^^

The
multi-window scheme is aimed at the long-baseline interferometer case, where
it is implicitly assumed that the windows are sufficiently far apart that they can be simulated with statistically independent
tweeter screens. With this assumption the tweeter screens for different windows are
generated independently of each other, but the windows share a common woofer screen to provide the large-scale correlations between apertures, so for example constraining the maximum fringe phase excursion as a function of interferometer baseline. The code does
not check if the tweeter screens are overlapping and simply 
generates uncorrelated tweeter screens in every  case.


Note that the code will also be inaccurate if two windows are
upwind/downwind of one another, as the “frozen” flow will not be correct:
the tweeter screen over one aperture will not pass over the other aperture.
However, this is less of a concern as the evidence for frozen flow over
large timescales is limited, and there are theoretical reasons
to believe that the frozen flow hypothesis is only valid over short
timescales for short-scale-length turbulent eddies, i.e. the components which are simulated by the tweeter screen.

If you want to simulate the effects of two close-together apertures, it may be best in some cases to situate them both inside  a single rectangular window.


Advanced parameters
^^^^^^^^^^^^^^^^^^^

The value of the outer scale of the optical turbulence is set by the parameter `L0`. This parameter is not well-constrained by measurement, but is thought to be in the range of order 10-100 meters for typical astronomical observatories. The default value of this parameter is 7000, which is 1000 times as large as the default value for `r0`. So for :math:`r_0=60` cm then the outer scale simulated is 600m, which will likely give a conservative estimate of the large-scale fluctuations. At optical wavelengths where the physical size of :math:`r_0` may be 10\,cm, then this set of parameter settings corresponds to an outer scale of 100m which is more realistic. 

The sizes of the FFTs for the woofer and tweeter screens can be adjusted for different situations if required. The largest spatial scale simulated is given by `nfftTweeter*nfftWoofer/(2*frequencyOverlap)` tweeter pixels. The default values of `nfftTweeter`=256, `nfftWoofer=256` and `frequencyOverlap=4`, give 8192 pixels as the maximum scale that can be represented in one dimension - the other dimension is semi-infinite. The value of `L0` could be up to :math:`2\pi` larger than this with relatively little effect because of the way :math:`L_0` is defined in the Von Karman model. Keeping the FFTs less than 1024 pixels on a side will mean that effects of typical CPU cache sizes on processing speed are minimised. Setting both `nfftWoofer` and `nfftTweeter` to 1024 gives a maximum spatial scale size of 131,072 pixels, comfortably exceeding the 4 orders of magnitude claimed in the abstract.

.. _`github site`: https://github.com/dbuscher/megascreen
.. [Buscher2016] Buscher, D.F. “Simulating Large Atmospheric Phase Screens Using a Woofer-Tweeter Algorithm.” Optics Express 24, no. 20 (October 3, 2016): 23566–71. https://doi.org/10.1364/OE.24.023566.
		 
.. toctree::
   :maxdepth: 3
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
