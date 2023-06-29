import astropy
import dill as pickle
from astropy import units
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import CCDData
from astropy.nddata import Cutout2D
import astropy.units as u
from astropy.modeling import models
import crapm2fs_process as m2fs
import os
import mycode
import specutils
from specutils import Spectrum1D
import scipy
from astropy.nddata import StdDevUncertainty
from astropy.visualization import quantity_support
import numpy as np
from specutils.fitting import fit_lines
from specutils.analysis import centroid
from specutils.analysis import fwhm
from specutils.analysis import line_flux
from specutils.analysis import equivalent_width
from specutils import SpectralRegion
from specutils.manipulation import extract_region
from specutils.fitting import fit_generic_continuum
from astropy.modeling import models,fitting
from specutils.fitting import estimate_line_parameters
from ccdproc import Combiner
#matplotlib.use('pdf')
matplotlib.use('TkAgg')

initialize=False
find=False
trace=True

threshold_factor=25.#multiple of continuum residual rms threshold to impose for aperture detection in find_lines_derivative
n_lines=20#columns to combine when scanning across rows to identify apertures (as 'emission lines')
continuum_rejection_iterations=10#number of iterations of outlier rejection for fitting "continuum"
window=10#pixels, width of aperture window for fitting (gaussian) aperture profiles (perpendicular to spectral axis)
trace_step=n_lines#tracing step
trace_nlost=3#number of consecutive failures to trace before quitting
trace_shift_max=3.
trace_order=4
trace_rejection_iterations=10
trace_rejection_sigma=3.#largest rms deviation to accept in fit to aperture trace
trace_pix_min=250.#conservative limits for trace fit
trace_pix_max=1500.#conservative limits for trace fit

directory='/nfs/nas-0-9/mgwalker.proj/m2fs/'
m2fsrun='nov18'
datadir='/nfs/nas-0-9/mgwalker.proj/m2fs/m2fs.astro.lsa.umich.edu/data/NovDec2018/'

for ccd in ('b','r'):
    root=directory+m2fsrun+'_'+ccd+'_master_fibermap1'
    master_fibermap=root+'.fits'
    data=astropy.nddata.CCDData.read(master_fibermap)#format is such that data.data[:,0] has column-0 value in all rows

    if initialize:#find apertures in each column stack and save to pickle files
        columnspec_array0=m2fs.get_columnspec(data,trace_step,n_lines,continuum_rejection_iterations,threshold_factor,window)
        '''
        columnspec_array is array of columnspec objects.
        columnspec object contains: 
        columns: numbers of columns in original data frame that are stacked
        spec1d: 'spectrum' across stacked column, where each aperture appears as an 'emission line'
        pixel: value of pixel across stacked column 'spectrum', artificially given units of AA in order to comply with specutils requirements for fitting spectra
        continuum: parameters of 1dGaussian continuum fit to stacked column 'spectrum'
        rms: rms residuals around continuum fit (using only regions between fiber bundles)
        apertures_initial: initial aperture centers returned by specutils find_lines_derivative
        apertures_profile: contains parameters of fits to apertures detected in apertures_initial
        '''
        pickle.dump(columnspec_array0,open(root+'_columnspec_array.pickle','wb'))#save pickle to file

    if find:
        columnspec_array=pickle.load(open(root+'_columnspec_array.pickle','rb'))

        '''
        Make initial array of more precise aperture location fits in middle column stack
        '''
        middle_column=np.long(len(columnspec_array)/2)
        apertures_profile_middle=columnspec_array[middle_column].apertures_profile

        '''
        apertures_profile_middle object contains: 
        fit: parameters of fit to cross-aperture profile ('line' profile) *for middle column*
        subregion: pixel range of cross-aperture 'spectrum' used for fit *for middle column*
        realvirtual: value of 1 means aperture corresponds to real spectrum, value of 0 means virtual aperture used as place-holder (due to bad or unplugged fiber)
        '''

        '''
        display apertures initially identified in central column stack and allow user to add/delete apertures with screen input
        '''
        apertures_profile_middle=m2fs.fiddle_apertures(columnspec_array,middle_column,window,apertures_profile_middle)
        pickle.dump([apertures_profile_middle,middle_column],open(root+'_apertures_profile_middle.pickle','wb'))#save pickle to file

    if trace:
        columnspec_array=pickle.load(open(root+'_columnspec_array.pickle','rb'))
        apertures_profile_middle,middle_column=pickle.load(open(root+'_apertures_profile_middle.pickle','rb'))

        aperture_trace_array=[]
        for j in range(0,len(apertures_profile_middle.fit)):
            if apertures_profile_middle.realvirtual[j]==1:
                mean0=apertures_profile_middle.fit[j].mean.value
                trace_x=[]
                trace_y=[]

                nlost=0
                for k in reversed(range(0,middle_column)):
                    print(k)
                    means=np.array([columnspec_array[k].apertures_profile.fit[q].mean.value for q in range(0,len(columnspec_array[k].apertures_profile.fit))])
                    means=means[means==means]#remove any NaNs due to poor fits
                    dist=(mean0-means)**2
                    best=np.where(dist==np.min(dist))[0][0]
                    if dist[best]<trace_shift_max:
                        nlost=0
                        mean0=means[best]
                        trace_x.append(np.median(columnspec_array[k].columns))
                        trace_y.append(means[best])
                    else:
                        nlost=+1

                nlost=0
                mean0=apertures_profile_middle.fit[j].mean.value
                for k in range(middle_column,len(columnspec_array)):
                    print(k)
                    means=np.array([columnspec_array[k].apertures_profile.fit[q].mean.value for q in range(0,len(columnspec_array[k].apertures_profile.fit))])
                    means=means[means==means]#remove any NaNs due to poor fits
                    dist=(mean0-means)**2
                    best=np.where(dist==np.min(dist))[0][0]
                    if dist[best]<trace_shift_max:
                        nlost=0
                        mean0=means[best]
                        trace_x.append(np.median(columnspec_array[k].columns))
                        trace_y.append(means[best])
                    else:
                        nlost=+1
                trace_x=np.array(trace_x)
                trace_y=np.array(trace_y)
                order=np.argsort(trace_x)
                trace_x=trace_x[order]
                trace_y=trace_y[order]
            #DO ANOTHER FIDDLE TO INTERACTIVELY CHOOSE MIN,MAX FOR FIT
                aperture_trace_func,aperture_trace_order,aperture_trace_pixel_min,aperture_trace_pixel_max=m2fs.get_aperture_trace(trace_x,trace_y,trace_order,trace_rejection_sigma,trace_pix_min,trace_pix_max)
                
                print(aperture_trace_func,aperture_trace_order)
#                aperture=aperture(xxx)###write aperture class above that holds trace_func, minimum and maximum pixels of fit,  function for width of aperture vs column, etc.)

    np.pause()
