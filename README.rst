Atmospheric phase screen generator using a "woofer-tweeter" algorithm
======================================================================

This Python 3 module implements the algorithm described in the paper “Simulating Large Atmospheric Phase Screens Using a Woofer-Tweeter Algorithm.” Optics Express 24, no. 20 (October 3, 2016): 23566–71. https://doi.org/10.1364/OE.24.023566.

The algorithm superposes two phase screens containing high-frequency ("tweeter") and low-frequency ("woofer") perturbations. An example result can be seen in the figure below, where the x- and y-direction structure functions of the individual screens are compared to the theoretical structure function, and their sum is shown on the right. The sum clearly matches theoretical predictions to a high degree of accuracy.

.. image:: tests/component_sf.png

Requirements
------------

The module runs under Python3 and requires ``numpy`` and ``scipy``.

For running the test and example code, the ``astropy`` and ``joblib`` libraries are used. There may be some residual dependencies on the ``pois`` library also, but this dependency is being eliminated by including the neccessary functions (e.g. Zernike polynomials) as part of this library.  


Installation
------------

On unix-like systems do

::

    pip3 install MegaScreen

or if that does not work because of file permission errors, then
::

    sudo pip3 install MegaScreen

 
Alternatively download and unpack a copy of this repository, change the working directory to this directory and then use

::

    pip3 install -e .


This should install the package into the Python path.

Basic usage
-----------

The main interface to the package is the MegaScreen function, which is a Python generator. The following code will generate 10000 snapshots of "frozen turbulence" blowing past a single rectangular window, and process them using a user-defined function called "process":

.. code:: python

    import MegaScreen
    
    for phaseScreen in MegaScreen.MegaScreen(numIter=10000):
	process(phaseScreen)

See the documentation for more detail.

Licencing
---------

The code is licenced under the Mozilla Public License Version 2.0
 (see `LICENCE`_).

.. _LICENCE: LICENCE
	   
