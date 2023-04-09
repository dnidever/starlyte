
## Create simple stellar population (SSP) synthetic spectra

import os
import numpy as np
from glob import glob
from dlnpyutils import utils as dln
from astropy.table import Table
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic_2d
from . import ferre,utils

# Get FERRE grid info
headerfiles = glob(utils.datadir()+'ssp*grid*.hdr')  # Get grid header files
dt = [('filename',str,200),('teffrange',float,2),('loggrange',float,2)]
GRIDINFO = np.zeros(len(headerfiles),dtype=np.dtype(dt))
# Get header Teff/logg ranges
for i in range(len(headerfiles)):
    hinfo = ferre.gridinfo(headerfiles[i])
    GRIDINFO['filename'][i] = headerfiles[i]
    GRIDINFO['teffrange'][i] = [np.min(hinfo['TEFF']),np.max(hinfo['TEFF'])]
    GRIDINFO['loggrange'][i] = [np.min(hinfo['LOGG']),np.max(hinfo['LOGG'])]
# Sort them by temperature
si = np.argsort(GRIDINFO['teffrange'][:,0])
GRIDINFO = GRIDINFO[si]

def ferre_interp(pars):
    """ Interpolated spectra in grid using FERRE."""
    # pars = [teff,logg,metal,alpha]
    # can be 2D
    if pars.ndim==1:
        pars = np.atleast_2d(pars)
    npars = pars.shape[0]
    
    # cool grid 1:     3500 <= Teff <=  6000, 0.0 <= logg <= 5.0
    # cool grid 2:     6000 <= Teff <=  8000, 1.0 <= logg <= 5.0
    # medium grid 1:   8000 <= Teff <= 11750, 2.0 <= logg <= 5.0
    # medium grid 2:  11500 <= Teff <= 19500, 2.5 <= logg <= 5.0
    # hot grid 1:     19000 <= Teff <= 26000, 3.0 <= logg <= 5.0
    # hot grid 2:     26000 <= Teff <= 31000, 3.5 <= logg <= 5.0
    # hot grid 3:     31000 <= Teff <= 39000, 4.0 <= logg <= 5.0
    # hot grid 4:     39000 <= Teff <= 49000, 4.5 <= logg <= 5.0        
    
    #teffranges = np.array([[3500,6000],[6000,8000],[8000,11750],[11500,19500],
    #                       [19000,26000],[26000,31000],[31000,39000],[39000,49000]]).astype(float)
    #loggranges = np.array([[0.0,5.0],[1.0,5.0],[2.0,5.0],[2.5,5.0],[3.0,5.0],
    #                       [3.5,5.0],[4.0,5.0],[4.5,5.0]])
    #gridfiles = ['coolgrid1','coolgrid2','mediumgrid1','mediumgrid2',
    #             'hotgrid1','hotgrid2','hotgrid3','hotgrid4']
    ngrids = len(GRIDINFO)
    
    # Loop over each star and assign it to a grid based on teff/logg
    gridindex = np.zeros(npars,int)
    newpars = pars.copy()
    for i in range(npars):
        pars1 = pars[i,:]  # teff, logg, metal, alpha
        teff1 = pars1[0]
        logg1 = pars1[1]
        if teff1 < np.min(GRIDINFO['teffrange']):
            gridindex[i] = 0
        elif teff1 > np.max(GRIDINFO['teffrange']):
            gridindex[i] = ngrids-1
        else:
            ind, = np.where((teff1 >= GRIDINFO['teffrange'][:,0]) &
                            (teff1 <= GRIDINFO['teffrange'][:,1]))
            gridindex[i] = ind[0]
                            
        # Deal with out of bounds teff/logg values
        #  set to the boundary
        teffr = GRIDINFO['teffrange'][gridindex[i],:]
        loggr = GRIDINFO['loggrange'][gridindex[i],:]        
        if teff1 < teffr[0]:
            newpars[i,0] = teffr[0]+10
        elif teff1 > teffr[1]:
            newpars[i,0] = teffr[1]-10
        if logg1 < loggr[0]:
            newpars[i,1] = loggr[0]+0.02
        elif logg1 > loggr[1]:
            newpars[i,1] = loggr[1]-0.02


    print('Interpolating '+str(npars)+' spectra with FERRE')
            
    # Loop over gridfiles and get all of the spectra
    # that fall in its teff range
    count = 0
    flux = None
    for i in range(ngrids):
        ind, = np.where(gridindex == i)
        nind = len(ind)
        if nind==0:
            continue
        pars2 = pars[ind,:]
        newpars2 = newpars[ind,:]        

        gfile = os.path.basename(GRIDINFO['filename'][i])
        print(i+1,len(ind),gfile)
        
        # Interpolate in the FERRE grid
        fout = ferre.interp(newpars2,grid=gfile)

        # Wavelength solution
        if flux is None:
            wave = fout['wave']
            npix = len(wave)
            flux = np.zeros((npars,npix),float) 

        # Fill in the spectra
        flux[count:count+nind,:] = fout['flux']
        totflux = np.sum(flux,axis=1)

        # some spectra are bad, all -4093520 values
        
        count += nind
        
    # Rescale the blackbody if outside the temperature range
    #bd, = np.where(pars[:,0] != newpars[:,0])
    #nbd = len(bd)
    #for i in range(nbd):
    #    ind = bd[i]
    #    teff0 = pars[ind,0]
    #    teff1 = newpars[ind,0]
    #    #bb0 = planck(wave,teff0)
    #    #bb1 = planck(wave,teff1)        
        
        
    # Put everything together
    out = {'wave':wave,'flux':flux}
    return out


def sspgrid(ages,metals,alphas,outdir='./',clobber=False):
    """  Run a grid of SSP spectra."""

    nage = len(ages)
    nmetal = len(metals)
    nalpha = len(alphas)
    nspectra = nage * nmetal * nalpha
    print('Creating {:d} SSP spectra'.format(nspectra))
    
    # age, metallicity, alpha abundance
    dt = [('age',float),('metal',float),('alpha',float),('luminosity',float)]
    pars = np.zeros((nage,nmetal,nalpha),dtype=np.dtype(dt))
    count = 0
    for i in range(nage):
        age = ages[i]
        for j in range(nmetal):
            metal = metals[j]
            for k in range(nalpha):
                alpha = alphas[k]
                print('{:d} (age,[M/H],[alpha/M])=({:.4f},{:.4f},{:.4f})'.format(count+1,age,metal,alpha))

                outfile = 'ssp_a{:.3f}m{:+.2f}a{:+.2f}.fits'.format(11,metal,alpha)
                if os.path.exists(outdir+outfile) and clobber==False:
                    print(outfile+' already exists and clobber not set')
                    spectrum,hd = fits.getdata(outdir+outfile,header=True)
                    wave = np.arange(hd['naxis1'])*hd['cdelt1']+hd['crval1']
                else:  
                    wave,spectrum = ssp(age,metal,alpha)
                    
                # Start grid
                if count == 0:
                    nwave = len(wave)
                    spectra = np.zeros((nage,nmetal,nalpha,nwave),float)

                # Plug in the spectrum
                spectra[i,j,k,:] = spectrum
                dw = wave[1]-wave[0]
                luminosity = np.sum(spectrum*dw)
                pars['age'][i,j,k] = age
                pars['metal'][i,j,k] = metal
                pars['alpha'][i,j,k] = alpha
                pars['luminosity'][i,j,k] = luminosity

                # Write to a file in case it crashes at some point
                if os.path.exists(outdir+outfile)==False or clobber:
                    hdu = fits.PrimaryHDU(spectrum)
                    hdu.header['CRVAL1'] = wave[0]
                    hdu.header['CRPIX1'] = 1
                    hdu.header['CDELT1'] = np.round(wave[1]-wave[0],6)
                    hdu.header['DC-FLAG'] = 0
                    hdu.writeto(outdir+outfile)
                
                print(' ')
                
                count += 1

    return wave,spectra,pars

                
def ssp(age,metal,alpha,alliso=None):
    """ Make SSP for a given age and metal."""

    print('Age = {:.4f} Gyr'.format(age))
    print('[M/H] = {:.2f}'.format(metal))
    print('[alpha/Fe] = {:.2f}'.format(alpha))    

    # Salaris correction
    metal_salaris = metal + np.log10(0.659*(10**alpha)+0.341)
    print('[M/H]_Salaris = {:.2f}'.format(metal_salaris))    
    
    # --- Isochrones ---
    if alliso is None:
        alliso = Table.read(utils.datadir()+'ssp_isochrones.fits.gz')
        for c in alliso.colnames: alliso[c].name = c.upper()
        #alliso = Table.read('/Users/nidever/isochrone/parsec_gaiaedr3_2mass/parsec_gaiaedr3_2mass.fits')
    umetal = np.unique(alliso['METAL'])
    uage = np.unique(alliso['AGE']/1e9)
    bestmetal,bestmetalind = dln.closest(umetal,metal_salaris)
    bestage,bestageind = dln.closest(uage,age)

    print('Closest isochrone values')
    print('Age = {:.4f} Gyr'.format(bestage))
    print('[M/H] = {:.2f}'.format(bestmetal))
    print('[alpha/Fe] = {:.2f}'.format(alpha))    
    
    # Get the isochrone
    isoind, = np.where((alliso['METAL']==bestmetal) & (alliso['AGE']/1e9==bestage))
    iso = alliso[isoind]
    #iso['AGE'] = 10**iso['LOGAGE']
    
    # --- Create synthetic photometry ---
    #bands = ['GMAG','G_BPMAG','G_RPMAG']
    stab = synth(iso,[],minlabel=1,maxlabel=9,nstars=100000)
    
    # Make 2-D bins so we don't have to interpolate so many spectra
    data = stab['LOGTE'],stab['LOGG']
    # Teff : 2000 to 60,000 K
    # logg : -2.0 to 6.0
    bins = np.arange(3.30,4.78,0.01),np.arange(-2.0,6.0,0.10)    
    binout = binned_statistic_2d(stab['LOGTE'],stab['LOGG'],stab['LOGG'],
                                 statistic='count',bins=bins,expand_binnumbers=True)
    result, teff_edge, logg_edge, binnumber = binout
    teff_center = 10**((teff_edge[0:-1]+teff_edge[1:])*0.5)
    logg_center = (logg_edge[0:-1]+logg_edge[1:])*0.5
    binid = np.char.array(binnumber[0,:]).astype(str) + '-' + np.char.array(binnumber[1,:]).astype(str)
    bin_index = dln.create_index(binid)
    npars = len(bin_index['value'])
    totluminosity = np.zeros(npars,float)  # solar luminosities
    pars = np.zeros((npars,4),float)
    pars[:,2] = metal
    pars[:,3] = alpha
    for i,ubinid in enumerate(bin_index['value']):
        ind = bin_index['index'][bin_index['lo'][i]:bin_index['hi'][i]+1]
        nind = len(ind)
        totluminosity[i] = np.sum(10**stab['LOGL'][ind])
        ubin2d = np.array(ubinid.split('-')).astype(int)
        # the indices give the upper value of the "edge"
        # so we need to subtract 1 to get the correct "center"
        pars[i,0] = teff_center[ubin2d[0]-1]  
        pars[i,1] = logg_center[ubin2d[1]-1]    

        
    #from dlnpyutils import plotting as pl
    #import pdb; pdb.set_trace()

        
    # Plot the binned synthetic stars
    #from dlnpyutils import plotting as pl
    #pl.display(result.T,teff_edge,logg_edge,log=True,xflip=True,yflip=True)
    
    # Use FERRE to interpolate in the grid to the
    # synthetic photometric Teff, logg, [Fe/H], and [alpha/Fe]
    fout = ferre_interp(pars)
    wave = fout['wave']   # 1-D wavelength array
    flux = fout['flux']   # 2-D flux array [Nspec,Nwavelength]
    npix = len(wave)
    
    # Convert each spectrum to a total flux of 1 Lsun
    dw = 0.1
    lsun = 3.846e33   # erg/s
    totflux = np.sum(flux*dw,axis=1)
    scalefactor = lsun / totflux
    scaleflux = flux * scalefactor.reshape(-1,1)
    
    # Combine all of the spectra with appropriate luminosity weighting
    spectrum = np.sum(scaleflux * totluminosity.reshape(-1,1),axis=0)

    # Units
    # synspec returns fluxed spectra in erg/cm2/s/A
    # PARSEC luminosities are in solar luminosity units (3.846 x 10^33 erg/s)
    # each pixel is 0.1 A

    # np.sum(spectrum*dw) / lsun
    print('Final spectrum luminosity = {:8.2f} Lsun'.format(np.sum(spectrum*dw) / lsun))
    print('Total synthetic photometry luminosity = {:8.2f} Lsun'.format(np.sum(totluminosity)))
    
    return wave,spectrum
    
    
def synth(iso,bands,nstars=None,totmass=None,minlabel=1,maxlabel=8,minmass=0,maxmass=1000,
          columns=['AGE','METAL','MINI','LOGTE','LOGG','LOGL','LABEL']):
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
