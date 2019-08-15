================
 Advanced usage
================

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

If you want to simulate the effects of two close-together apertures, it may be best in some cases to situate them both inside  a single rectangular window of appropriate shape and size.

Polychromatic simulations
^^^^^^^^^^^^^^^^^^^^^^^^^
For simulations covering a not too large range of wavelengths, it is a good approximation to assume that the refractive index perturbations are independent of wavelength, so the phase perturbations at one wavelength can be derived from the phase perturbations at another wavelength by scaling the phase by the ratio of wavelengths. Note that in diffraction calculations, for example speckle image formation, the spatial scale of the diffraction pattern is wavelength-dependent, so this must also be taken into account, for example by magnifying or demagnifying the images as a function of wavelength.

For accurate simulations at infrared wavelengths, it may be necessary to generate independent phase screens representing the water-vapour fluctuations and the "dry" fluctuations, and combine these in different proportions depending on the relative variations of the refractive index of water and air with wavelength. 

Outer scale and woofer screen size
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The value of the outer scale of the optical turbulence is set by the function parameter `L0`. The corresponding physical quantity :math:`L_0` is not well-constrained experimentally, but is thought to be in the range of order 10-100 meters for atmospheric conditions typical over astronomical observatories. The default value of `L0` is 7000, which is 1000 times as large as the default value for `r0`. In the infrared case corresponding to :math:`r_0=60` cm then the outer scale simulated is 600m, and so this is likely to give an over-estimate of the wavefront fluctuations on the largest physical scales. At optical wavelengths where the physical size of :math:`r_0` may be :math:`10\,\text{cm}`, then this set of parameter settings corresponds to an outer scale of 100m which is more realistic. 

The size of a woofer "tile" is given by `nfftTweeter*nfftWoofer/(2*frequencyOverlap)` tweeter pixels. The default values of `nfftTweeter=256`, `nfftWoofer=256` and `frequencyOverlap=4`, gives woofer tiles which are  8192 tweeter pixels across.
The value of `L0` can be perhaps 2-3 times larger than this value with relatively little loss in accuracy because of the way :math:`L_0` is defined in the Von Karman model: the spectrum of fluctuations flattens out for wavefront corrugations whose wavelength is larger than :math:`L_0/2\pi`.

The `MegaScreen()` function generates the woofer screen as a semi-infinite "strip" by splicing together woofer tiles. All the windows for a long-baseline interferometer must fit inside this "strip". If larger spacings than this are needed, the values of `nfftWoofer` and/or `nfftTweeter` can be increased.

For efficiency, `nfftTweeter` and `nfftWoofer` should be powers of 2 and should be larger than 64.
Keeping the FFTs less than 1024 pixels on a side will mean that effects of typical CPU cache sizes on processing speed are minimised. Setting both `nfftWoofer` and `nfftTweeter` to 1024 gives a maximum spatial scale size of 131,072 tweeter pixels, comfortably exceeding the 4-orders-of-magnitude range of scales  claimed in the paper abstract. With these settings, and using `r0=7.0` to simulate conditions where :math:`r_0=10\,\text{cm}` then interferometer baselines up to 1.8\,km can be simulated. 
