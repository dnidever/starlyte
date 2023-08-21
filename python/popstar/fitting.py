import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from doppler import spec1d


def fitsingle(spec,grid):
    """
    Fit a single SSP spectrum to an observed spectrum.

    Parameters
    ----------
    spec : Spec1D object
       An observed spectrum.
    grid : SpecGrid object
       A SSP spectral grid.

    Returns
    -------
    out : table
      Table with the best-fit results.
    model : Spec1D
      Best-fit SSP model spectrum.

    Example
    -------

    out,model = fitginsle(spec,grid)

    """

    pass


def fitmulti(spec,grid):
    """
    Fit multiple SSP spectral components to an observed spectrum.

    Parameters
    ----------
    spec : Spec1D object
       An observed spectrum.
    grid : SpecGrid object
       A SSP spectral grid.

    Returns
    -------
    out : table
      Table with the best-fit results.
    model : Spec1D
      Best-fit SSP model spectrum.

    Example
    -------

    out,model = fitginsle(spec,grid)

    """

    pass
