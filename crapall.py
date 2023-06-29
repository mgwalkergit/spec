import astropy
from astropy import units
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import CCDData
import astropy.units as u
from astropy.modeling import models
import os
import mycode
import specutils
from specutils import Spectrum1D
import scipy
from astropy.nddata import StdDevUncertainty
from astropy.visualization import quantity_support
import numpy as np
from specutils.fitting import find_lines_threshold
from specutils.fitting import find_lines_derivative
from specutils.fitting import fit_lines
from specutils.analysis import centroid
from specutils.analysis import fwhm
from specutils.analysis import line_flux
from specutils.analysis import equivalent_width
from specutils import SpectralRegion
from specutils.manipulation import extract_region
matplotlib.use('pdf')
matplotlib.use('TkAgg')

n_apertures=128
min_separation=5.#pixels, minimum separation between aperture centers
slices=4
window=5#pixels, width of aperture window for fitting (gaussian) aperture profiles (perpendicular to spectral axis)
stddev=2.#pixels, initialization of standard deviation parameter for estimating stddev of aperture profile (perpendicular to spectral axis)
amplitude=1000.#electrons, initialization of amplitude for estimating gaussian fit to aperture profile (perpendicular to spectral axis)

data=astropy.nddata.CCDData.read('/nfs/nas-0-9/mgwalker.proj/m2fs/nov18/m2fs.astro.lsa.umich.edu/data/NovDec2018/ut20181126/b0760_stitched.fits')#,hdu_uncertainty=StdDevUncertainty)

n_cols=np.shape(data)[1]
slice=np.linspace(0,n_cols,slices,dtype='int')

meanvals=[]
stdvals=[]
eqwvals=[]

plt.imshow(data,vmin=0,vmax=1000,cmap='gray')
for i in slice:
    column=data.data[:,i]
    column_uncertainty=data.uncertainty[:,i]
    column_mask=data.mask[:,i]

    lamb=(np.arange(len(column),dtype='float'))*u.AA#unit is pixels, but specutils apparently can't handle that, so we lie and say Angs.
    spec1d=0.
    aperture_centers=0.
    apertures=0.
    spec1d=Spectrum1D(spectral_axis=lamb,flux=column*u.electron,uncertainty=column_uncertainty,mask=column_mask)

    apertures=find_lines_threshold(spec1d,noise_factor=5)

###identify false-positive apertures as those that violate minimum aperture separation
    aperture_centers=apertures['line_center'].value
    plt.scatter(np.zeros(len(aperture_centers))+np.float(i),aperture_centers,s=3,color='r')
    dcenters=[]
    for k in range(0,len(aperture_centers)-1):
        dcenters.append(aperture_centers[k+1]-aperture_centers[k])
    dcenters=np.array(dcenters)
    keep_apertures=np.concatenate(((dcenters>min_separation),np.array([True])))
    apertures=apertures[keep_apertures]
    keep_apertures=np.array(keep_apertures)

    print('found '+str(len(apertures))+' apertures in column '+str(i))

#    gs=plt.GridSpec(10,10) # define multi-panel plot
#    gs.update(wspace=0,hspace=0) # specify inter-panel spacing
#    fig=plt.figure(figsize=(6,6)) # define plot size
#    ax1=fig.add_subplot(gs[0:4,0:10])
#    ax2=fig.add_subplot(gs[6:10,0:4])
#    ax3=fig.add_subplot(gs[6:10,6:10])

#    ax1.plot(spec1d.spectral_axis,spec1d.flux,lw=0.1,color='k')
    eqw=[]
    std=[]
    mean=[]
    for j in range(0,len(apertures)):
        g_init=models.Gaussian1D(amplitude=amplitude*u.electron,mean=apertures[j][0].value*u.AA,stddev=stddev*u.AA)
        g_fit=fit_lines(spec1d,g_init,window=window*u.AA)
        ax1.axvline(x=apertures[j][0].value,color='b',lw=0.1,alpha=0.7)
        val1=apertures[j][0].value-window/2.
        val2=apertures[j][0].value+window/2.
        sub_region=SpectralRegion(val1*u.AA,val2*u.AA)
        sub_spectrum=extract_region(spec1d,sub_region)
        y_fit=g_fit(sub_spectrum.spectral_axis)
        std.append(g_fit.stddev.value)
        mean.append(g_fit.mean.value)
        eqw.append(equivalent_width(spec1d,regions=sub_region).value)
        ax1.plot(sub_spectrum.spectral_axis,y_fit,color='r',lw=0.05,alpha=0.7)
        diff=g_fit.mean.value-apertures[j][0].value
    mean=np.array(mean)
    std=np.array(std)
    eqw=np.array(eqw)

    meanvals.append(mean)
    stdvals.append(std)
    eqwvals.append(eqw)

#    if len(apertures)==129:
#        np.pause()
#    ax2.scatter(np.arange(len(apertures)),eqw,s=5,color='k')
#    ax3.scatter(np.arange(len(apertures)),std,s=5,color='k')
#    plt.savefig('test'+str(i)+'.pdf',dpi=200)
#    plt.show()
#    plt.close()

meanvals=np.array(meanvals)
stdvals=np.array(stdvals)
eqwvals=np.array(eqwvals)

shite=np.transpose(meanvals)

#plt.imshow(data,vmin=0,vmax=1000,cmap='gray')
#for i in range(0,len(shite)):
#for i in range(0,16):
#    plt.plot(slice,shite[i],color='r',lw=0.2)
#plt.savefig('test.pdf',dpi=200)
#plt.show()
#plt.close()
plt.show()
plt.close()
