Atmospheric phase screen generator using a "woofer-tweeter" algorithm
======================================================================

This Python 3 module implements the algorithm described in my paper “Simulating Large Atmospheric Phase Screens Using a Woofer-Tweeter Algorithm.” Optics Express 24, no. 20 (October 3, 2016): 23566–71. https://doi.org/10.1364/OE.24.023566.

The algorithm superposes two phase screens containing high-frequency ("tweeter") and low-frequency ("woofer") perturbations. The result can be seen in the figure below, where the x- and y-direction structure functions of the individual screens are compared to the theoretical structure function, and their sum is shown on the right. The sum clearly matches theoretical predictions to a high degree of accuracy.

.. image:: tests/component_sf.png
