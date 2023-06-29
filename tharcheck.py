import astropy
import dill as pickle
from astropy import units
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import CCDData
from astropy.nddata import Cutout2D
from astropy.nddata import NDDataRef
import astropy.units as u
from astropy.modeling import models
import m2fs_process as m2fs
import os
from os import path
import mycode
import specutils
from specutils.spectra import Spectrum1D
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
from scipy import interpolate
from astropy import time, coordinates as coord, units as u
from astropy.coordinates import SkyCoord, EarthLocation
from mpl_toolkits.axes_grid1 import make_axes_locatable
#matplotlib.use('pdf')
matplotlib.use('TkAgg')

shite=False
trim=False#define boundaries of useful spectral region
initialize=False
find=False#find apertures in image
trace_all=False#trace all apertures
trace_edit=False#edit traces of individual apertures
make_image=False#generate PDF of data frame with aperture traces overlaid
apmask=False#mask apertures to create 2d data frame containing only extra-aperture light
apflat=False
apflatcorr=False
scatteredlightcorr=False#fit 2d function to extra-aperture light and subtract from data frame
extract1d_flat=False#extract 1d spectra for flat frames
extract1d_thar=False#extract 1d spectra for thar frames
extract1d_sci=False#extract 1d spectra for science frames
id_lines_template=False#identify lines in thar template and fit wavelength solution
id_lines_translate=True
id_lines_check=False
tharcheck=False
wavcal=False
cr_reject=False
stack_twilight=False
throughputcorr=False#perform throughput correction (wavelength-dependent)
plugmap=False
skysubtract=False
stack_frames=False
sky_target_check=False
writefits=False
overwrite=False#overwrite previous results

#linelist_file='/nfs/nas-0-9/mgwalker.proj/m2fs/m2fs_config1b_thar_list'
lco=coord.EarthLocation.from_geodetic(lon=-70.6919444*u.degree,lat=-29.0158333*u.degree,height=2380.*u.meter)
threshold_factor=25.#multiple of continuum residual rms threshold to impose for aperture detection in find_lines_derivative
n_lines=20#columns to combine when scanning across rows to identify apertures (as 'emission lines')
columnspec_continuum_rejection_low=-5.
columnspec_continuum_rejection_high=1.
columnspec_continuum_rejection_iterations=10#number of iterations of outlier rejection for fitting "continuum"
columnspec_continuum_rejection_order=10
profile_rejection_iterations=10#number of iterations of outlier rejection for fitting profile amplitude and sigma
profile_nsample=50#number of points along spectral- (x-) axis at which to measure profile amplitude and sigma before performing fit of amplitude(x) and sigma(x)
profile_order=4#order of polynomial used to fit profile amplitude and sigma as functions of pixel along dispersion direction
window=10#pixels, width of aperture window for fitting (gaussian) aperture profiles (perpendicular to spectral axis)
trace_step=n_lines#tracing step
trace_nlost_max=2
trace_shift_max=1.5
trace_order=4
trace_rejection_iterations=10
trace_rejection_sigma=3.#largest rms deviation to accept in fit to aperture trace
trace_rejection_iterations=10
id_lines_continuum_rejection_low=-5.
id_lines_continuum_rejection_high=1.
id_lines_continuum_rejection_sigma=3.
id_lines_continuum_rejection_iterations=10#number of iterations of outlier rejection for fitting "continuum"
id_lines_continuum_rejection_order=10
id_lines_threshold_factor=10.#multiple of continuum residual rms threshold to impose for aperture detection in find_lines_derivative
id_lines_window=5.#pixels, width of aperture window for fitting (gaussian) line profiles in arc spectra
id_lines_order=5#order of wavelength solution
id_lines_tol_angs=0.05#tolerance for finding new lines to add from linelist (Angstroms)
id_lines_tol_pix=2.#tolerance for matching lines between template and new spectrum (pixels)
id_lines_minlines_hires=25#mininum number of ID'd lines for acceptable wavelength solution (less than this, and cannot fit reliable throughput correction and beyond)
id_lines_minlines_medres=15#mininum number of ID'd lines for acceptable wavelength solution (less than this, and cannot fit reliable throughput correction and beyond)
resolution_order=1
resolution_rejection_iterations=10
scatteredlightcorr_order=4
scatteredlightcorr_rejection_iterations=10
scatteredlightcorr_rejection_sigma=3.
id_lines_translate_add_lines_iterations=5
extract1d_aperture_width=3.#maximum (half-)width of aperture for extraction (too large and we get weird edge effects)
throughputcorr_continuum_rejection_low=-1.
throughputcorr_continuum_rejection_high=3.
throughputcorr_continuum_rejection_iterations=5#number of iterations of outlier rejection for fitting "continuum"
throughputcorr_continuum_rejection_order=4
cr_rejection_low=-2.
cr_rejection_high=3.
cr_rejection_order=4
cr_rejection_iterations=5
cr_rejection_tol=5.#multiple of rms residual above fit to flag as CR
cr_rejection_collateral=2#number of pixels adjacent to CR-flagged pixel to mask
hires_exptime=29.
medres_exptime=10.
use_flat=True

directory='/nfs/nas-0-9/mgwalker.proj/m2fs/'
m2fsrun='may18'
datadir=m2fs.get_datadir(m2fsrun)

edit_header_filename=[]
edit_header_keyword=[]
edit_header_value=[]
with open(directory+'edit_headers') as f:
    data=f.readlines()[0:]
for line in data:
    p=line.split()
    edit_header_filename.append(p[0])
    edit_header_keyword.append(p[1])
    edit_header_value.append(p[2])
edit_header_filename=np.array(edit_header_filename)
edit_header_keyword=np.array(edit_header_keyword)
edit_header_value=np.array(edit_header_value)
for i in range(0,len(edit_header_filename)):
    print('editing header for ',edit_header_filename[i])
    fits.setval(edit_header_filename[i]+'_stitched.fits',edit_header_keyword[i],value=edit_header_value[i])

utdate=[]
file1=[]
file2=[]
flatfile=[]
tharfile=[]
field_name=[]
scifile=[]
fibermap_file=[]
fiber_changes=[]

with open(directory+m2fsrun+'_science_raw') as f:
    data=f.readlines()[0:]
for line in data:
    p=line.split()
    if p[0]!='none':
        utdate.append(str(p[0]))
        file1.append(int(p[1]))
        file2.append(int(p[2]))
        flatfile.append(p[3])
        tharfile.append(p[4])
        field_name.append(p[5])
        scifile.append(p[6])
        fibermap_file.append(p[7])
        fiber_changes.append(p[8])
utdate=np.array(utdate)
file1=np.array(file1)
file2=np.array(file2)
flatfile=np.array(flatfile)
tharfile=np.array(tharfile)
field_name=np.array(field_name)
scifile=np.array(scifile)
fibermap_file=np.array(fibermap_file)
fiber_changes=np.array(fiber_changes)

flatfile0=[]
tharfile0=[]
scifile0=[]
allfile0=[]
fiber_changes0=[]
for i in range(0,len(tharfile)):
    flatfile0.append(flatfile[i].split('-'))
    tharfile0.append(tharfile[i].split('-'))
    scifile0.append(scifile[i].split('-'))
    allfile0.append(flatfile[i].split('-')+tharfile[i].split('-')+scifile[i].split('-'))
    fiber_changes0.append(fiber_changes[i].split(','))
flatfile0=np.array(flatfile0,dtype='object')
tharfile0=np.array(tharfile0,dtype='object')
scifile0=np.array(scifile0,dtype='object')
allfile0=np.array(allfile0,dtype='object')
#fiber_changes0=np.array(fiber_changes0,dtype='str')

hires_npoints=[]
medres_npoints=[]
hires_rms=[]
medres_rms=[]
for i in range(0,len(utdate)):
    for ccd in ('r','b'):
        for k in range(0,len(tharfile0[i])):
            root=datadir+utdate[i]+'/'+ccd+str(flatfile[i]).zfill(4)
            root0=datadir+utdate[i]+'/'+ccd+str(tharfile0[i][k]).zfill(4)
            id_lines_array_file=root0+'_id_lines_array.pickle'
            columnspec_array_file=root+'_columnspec_array.pickle'
            apertures_profile_middle_file=root+'_apertures_profile_middle.pickle'
            aperture_array_file=root+'_aperture_array.pickle'
            image_file=root+'_apertures2d.pdf'
            find_apertures_file=root+'_find_apertures.pdf'
            apmask_file=root+'_apmask.pickle'
            apflat_file=root+'_apflat.pickle'
            apflat_residual_file=root+'_apflat_residual.pickle'
            extract1d_array_file=root0+'_extract1d_array.pickle'
            throughput_continuum_file=root+'_throughput_continuum.pickle'

            columnspec_array_exists=path.exists(columnspec_array_file)
            apertures_profile_middle_exists=path.exists(apertures_profile_middle_file)
            aperture_array_exists=path.exists(aperture_array_file)
            apmask_exists=path.exists(apmask_file)
            apflat_exists=path.exists(apflat_file)
            apflat_residual_exists=path.exists(apflat_residual_file)
            id_lines_array_exists=path.exists(id_lines_array_file)

            data_file=root0+'_stitched.fits'
            data=astropy.nddata.CCDData.read(data_file)#format is such that data.data[:,0] has column-0 value in all rows

            filtername=data.header['FILTER']
            print(filtername)
            if filtername=='Mgb_Rev2':
                filtername='Mgb_HiRes'
            if filtername=='Mgb_HiRes':
                id_lines_minlines=id_lines_minlines_hires
            if filtername=='Mgb_MedRes':
                id_lines_minlines=id_lines_minlines_medres
#        columnspec_array=pickle.load(open(columnspec_array_file,'rb'))
#        aperture_array=pickle.load(open(aperture_array_file,'rb'))
    
            print(id_lines_array_file,id_lines_array_exists)
            if id_lines_array_exists:
                id_lines_array=pickle.load(open(id_lines_array_file,'rb'))
                apertures_profile_middle,crap=pickle.load(open(apertures_profile_middle_file,'rb'))
                extract1d_array=pickle.load(open(extract1d_array_file,'rb'))
                for j in range(0,len(apertures_profile_middle.realvirtual)):
                    if apertures_profile_middle.realvirtual[j]:
                        shite=np.where(np.array([id_lines_array[q].aperture for q in range(0,len(id_lines_array))])==j+1)[0]
                        if len(shite)>0:
                            keep=np.where(extract1d_array[shite[0]].spec1d_mask==False)[0]
                            print(j,len(keep),id_lines_array[shite[0]].npoints,id_lines_array[shite[0]].rms)
                            if filtername=='Mgb_HiRes':
                                hires_npoints.append(id_lines_array[shite[0]].npoints)
                                hires_rms.append(id_lines_array[shite[0]].rms)
                            if filtername=='Mgb_MedRes':
                                medres_npoints.append(id_lines_array[shite[0]].npoints)
                                medres_rms.append(id_lines_array[shite[0]].rms)

hires_npoints=np.array(hires_npoints)
hires_rms=np.array(hires_rms)
medres_npoints=np.array(medres_npoints)
medres_rms=np.array(medres_rms)
gs=plt.GridSpec(7,7) # define multi-panel plot
gs.update(wspace=0,hspace=0) # specify inter-panel spacing
fig=plt.figure(figsize=(6,6)) # define plot size
ax1=fig.add_subplot(gs[1:6,0:5])
ax2=fig.add_subplot(gs[0:1,0:5])
ax3=fig.add_subplot(gs[1:6,5:6])

ax1.scatter(hires_npoints,np.log10(hires_rms),s=3,color='r')
ax1.scatter(medres_npoints,np.log10(medres_rms),s=3,color='b')
ax2.hist(hires_npoints,bins=100,range=[0,50],color='r')
ax2.hist(medres_npoints,bins=100,range=[0,50],color='b')
ax3.hist(np.log10(hires_rms),bins=100,range=[-3,1],color='r',orientation='horizontal')
ax3.hist(np.log10(medres_rms),bins=100,range=[-3,1],color='b',orientation='horizontal')
#plt.hist(medres_npoints,bins=100,range=[0,50],color='b')
ax1.set_xlim([0,50])
ax2.set_xlim([0,50])
ax1.set_ylim([-3,1])
ax3.set_ylim([-3,1])

plt.show()
plt.close()

