import numpy as np
from dlnpyutils import utils as dln
from astropy.table import Table
from scipy.interpolate import interp1d

def ssp(age,metal,alliso=None):
    """ Make SSP for a given age and metal."""

    print('Age = {:.2f} Gyr'.format(age))
    print('[M/H] = {:.2f}'.format(metal))
    
    # --- Isochrones ---
    if alliso is None:
        alliso = Table.read('/Users/nidever/isochrone/parsec_gaiaedr3_2mass/parsec_gaiaedr3_2mass.fits')
    umetal = np.unique(alliso['MH'])
    uage = np.unique(10**alliso['LOGAGE']/1e9)
    bestmetal,bestmetalind = dln.closest(umetal,metal)
    bestage,bestageind = dln.closest(uage,age)

    print('Closest isochrone values')
    print('Age = {:.2f} Gyr'.format(bestage))
    print('[M/H] = {:.2f}'.format(bestmetal))

    # Get the isochrone
    isoind, = np.where((alliso['MH']==bestmetal) & (10**alliso['LOGAGE']/1e9==bestage))
    iso = alliso[isoind]
    iso['AGE'] = 10**iso['LOGAGE']
    
    # --- Create synthetic photometry ---
    bands = ['GMAG','G_BPMAG','G_RPMAG']
    stab = synth(iso,bands,minlabel=1,maxlabel=9,nstars=100000)
    
    
    # Use FERRE to interpolate in the grid to the
    # synthetic photometric Teff, logg, [Fe/H], and [alpha/Fe]

    # Combine all of the spectra with appropriate luminosity weighting

    return stab
    
    
def synth(iso,bands,nstars=None,totmass=None,minlabel=1,maxlabel=8,minmass=0,maxmass=1000,
          columns=['AGE','METAL','MINI','MASS','LOGTE','LOGG','LOGL','LABEL']):
    """ Create synthetic population."""
    
    # By default us 1000 stars
    if nstars is None and totmass is None:
        nstars = 1000

    lab = ((iso['LABEL']>=minlabel) & (iso['LABEL']<=maxlabel))
    data = iso[lab]

    # AGE and METAL columns
    if 'AGE' not in data.colnames and 'LOGAGE' in data.colnames:
        data['AGE'] = 10**data['LOGAGE']
    if 'METAL' not in data.colnames and 'MH' in data.colnames:
        data['METAL'] = data['MH']
    
    ## Mass cuts
    massind, = np.where((data['MINI'] >= minmass) & (data['MINI'] <= maxmass))
    data = data[massind]
    ndata = len(data)
    if ndata==0:
        raise ValueError('No isochrone points left after mass cut')

    # Total mass input, figure out the number of stars we expect
    # for our stellar mass range
    if nstars is None and totmass is not None:
        nstars = np.ceil((np.max(data['INT_IMF'])-np.min(data['INT_IMF']))*totmass)
            
    # Initialize the output catalog
    out = Table()
    out['AGE'] = np.zeros(int(nstars),float)
    for c in columns:
        if c != 'AGE':
            out[c] = 0.0
    for n in bands:  # bands to interpolate
        out[n] = 0.0        

    # PDF, probability distribution function
    # int_IMF, which is the integral of the IMF under consideration (as selected in the form, in number of stars,
    # and normalised to a total mass of 1 Msun) from 0 up to the current M_ini. Differences between 2 values of
    # int_IMF give the absolute number of stars occupying that isochrone section per unit mass of stellar
    # population initially born, as expected for the selected IMF.
    pdf = np.maximum(np.diff(data['INT_IMF'].data),1e-8)
    pdf = np.hstack((0.0, pdf))
        
    # Create normalized index array (from 0-1)
    indx = np.arange(ndata).astype(float)/(ndata-1)
    # Create cumulative distribution function
    cdf = np.cumsum(pdf)
    cdf /= np.max(cdf)

    # Get the indices in the structure from the CDF
    #interp,cdf,indx,randomu(seed,nstars),newindx
    newindx = interp1d(cdf,indx)(np.random.rand(int(nstars)))
        
    # Interpolate all of the relevant columns
    for n in out.colnames:
        if n != 'INT_IMF':
            newval = interp1d(indx,data[n])(newindx)
            out[n] = newval
    out['TEFF'] = 10**out['LOGTE']
            
    return out
