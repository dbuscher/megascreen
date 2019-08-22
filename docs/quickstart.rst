=============
 Quick start
=============
Installation
^^^^^^^^^^^^

On unix-like systems do

::

    pip3 install MegaScreen

Note that we use `pip3` because this is a Python3 package. Alternative installation methods are given on the `github site`_.
    
Simple usage
^^^^^^^^^^^^

The main interface to the package is the `MegaScreen()` function, which is a Python generator: a function returning an iterator which can be used in a `for` loop. The following code will generate 10000 snapshots of "frozen turbulence" blowing past a single rectangular window, process them using a user-defined function called `process()`, and yield a list of all the results for later processing, e.g. evaluating the mean and standard deviation:

.. code:: python

    import MegaScreen

    results=[]
    for phaseScreen in MegaScreen.MegaScreen(numIter=10000):
	results.append(process(phaseScreen))


Example code using the `MegaScreen` function (and other internal functions in the package) is available under the `tests` directory on the `github site`_.

.. _`github site`: https://github.com/dbuscher/megascreen

