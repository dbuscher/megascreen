================
 Advanced usage
================

Simulating multiple turbulent layers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For simple applications, a single phase screen is sufficient. To add additional realism to the space-time correlation of the wavefonts, multiple screens with different wind speeds and directions are required to simulate the effect of multiple turbulent layers at different heights.
The
experimental evidence indicates that at least 3 different screens of approximately the same turbulent strength are needed to give the right level of decorrelation.

Each layer should be generated using an independent `MegaScreen` generator. In most cases, the phase screens so generated can simply be added together to yield a summed phase screen.


Each layer should have a different value for `theta`, the wind direction, and can have different values for `dx`, `r0` and `L0` as appropriate to the atmospheric model being simulated.
To derive appropriate :math:`r_0` values for the different layers, we need to consider the effect of multiple turbulent layers on the wavefront. The structure function of the wavefront perturbations depends on the integral of the turbulent strength
:math:`C_n^2` over the propagation path as
.. math::
   D_\Phi(r) = \left [ 2.91 \, k^2 \,   \int_{\text{path}} C_n^2(z') \, dz' \right ]r^{5/3},
where :math:`k=2\pi/\lambda` is the optical wavenumber of the radiation and it is assumed that :math:`r\llt L_0`. 
Summing together the perturbations due to several
statistically independent phase screens gives a wavefront perturbation whose structure function is the sum of the structure functions of the individual screens. This means we can break the integral into sub-integrals corresponding to traversal of each of the individual layers. 
From the definition of :math:`r_0` we get 
.. math::
   r_0 = \left [ \frac{2.91}{6.88} \, k^2 \,   \int_{\text{path}} C_n^2(z') \, dz' \right ]^{-3/5}.
Thus we set :math:`r_0` for each layer proportional to :math:`[\int  C_n^2 \,dz]^{-3/5}` for that layer, and
the effective  value of :math:`r_0` for the sum of all the screens is given by:
.. math::

   r_0 = \left[\sum_i (r_0)_i^{-5/3}\right]^{-3/5}
where :math:`(r_0)_i` is the Fried parameter of the :math:`i^{\text{th}}` layer.

If the field of view being considered is large, so that anisoplanatism effects are
non-negligible, then the phase screens need to be shifted with respect to one another by different amounts for different  locations within the field of view. Additionally, in situations where the "near-field" approximation breaks down, i.e. when scintillation effects are important, it may be necessary to simulate Fresnel propagation between different layers.


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
Values of `L0` up to this value will be modelled with reasonable accuracy because the spectrum of fluctuations flattens out at the lowest frequencies sampled by the woofer spectrum.

The `MegaScreen()` function generates the woofer screen as a semi-infinite "strip" by splicing together woofer tiles. All the windows for a long-baseline interferometer must fit inside this "strip". If larger spacings than this are needed, the values of `nfftWoofer` and/or `nfftTweeter` can be increased.

For efficiency, `nfftTweeter` and `nfftWoofer` should be powers of 2 and should be larger than 64.
Keeping the FFTs less than 1024 pixels on a side will mean that effects of typical CPU cache sizes on processing speed are minimised. Setting both `nfftWoofer` and `nfftTweeter` to 1024 gives a maximum spatial scale size of 131,072 tweeter pixels, comfortably exceeding the 4-orders-of-magnitude range of scales  claimed in the paper abstract. With these settings, and using `r0=7.0` to simulate conditions where :math:`r_0=10\,\text{cm}` then interferometer baselines up to 1.8\,km can be simulated. 
