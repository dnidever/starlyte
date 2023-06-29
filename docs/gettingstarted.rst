***************
Getting Started
***************



Introduction
============

This is a package to coadd images.  Currently, there is no fast python package to resample image.  Therefore, |coadd| relies on the `swarp` softare to perform the resampling.  This will be replaced in the future once a fast python resampling package is written.

There are five main modules:

 - :mod:`~gaussdecomp.driver`:  Decomposes all of the spectra in a datacube.
 - :mod:`~gaussdecomp.fitter`:  Does the actual Gaussian Decomposition.
 - :mod:`~gaussdecomp.cube`:  Contains the :class:`~gaussdecomp.cube.Cube` class for a data cube.
 - :mod:`~gaussdecomp.spectrum`:  Contains the :class:`~gaussdecomp.spectrum.Spectrum` class for a single spectrum.
 - :mod:`~gaussdecomp.utils`:  Various utility functions.

There is a class for data cubes called :class:`~gaussdecomp.cube.Cube` and a class for spectra called :class:`~gaussdecomp.spectrum.Spectrum`.

To fit a single spectrum you first need to create the Spectrum object.

.. code-block:: python

	from gaussdecomp import spectrum,fitter
	sp = spectrum.Spectrum(flux,vel)   # flux and velocity arrays
	out = fitter.gaussfit(sp)          # do the fitting

You can make a nice plot using :func:`~gaussdecomp.utils.gplot`.

.. code-block:: python

	from gaussdecomp import utils
	utils.gplot(vel,flux,par)
	
.. |gaussfitfig| image:: gaussfit.png
  :width: 800
  :alt: Gaussian Fit to Spectrum

|gaussfitfig|


	
To fit an entire datacube, you can either give the driver code a datacube object you have already created or give it a FITS filename.

.. code-block:: python

	from gaussdecomp import cube,driver
	# Load the cube first
	datacube = cube.Cube.read('mycube.fits')
	gstruc = driver.driver(datacube)

	# Give it the FITS filename
	gstruc = driver.driver('mycube.fits')


See the :doc:`examples` page for some examples of how Python |gaussdecomp| runs.
	
