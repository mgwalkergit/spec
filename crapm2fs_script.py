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

directory='/nfs/nas-0-9/mgwalker.proj/m2fs/'

with open(directory+'m2fs_apertures0.py') as f:
    apertures_text=f.readlines()[0:]

m2fsrun0=['feb14']

for i in range(0,len(m2fsrun0)):
    b_apertures_out='m2fs_b_apertures_'+m2fsrun0[i]+'.py'
    b_apertures_slurm_out='/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_b_apertures_'+m2fsrun0[i]+'.slurm'
    r_apertures_out='m2fs_r_apertures_'+m2fsrun0[i]+'.py'
    r_apertures_slurm_out='/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_r_apertures_'+m2fsrun0[i]+'.slurm'

    g1b=open(b_apertures_out,'w')
    for line in apertures_text:
        if 'enter_run_here' in line:
            g1b.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
        elif 'for ccd in ' in line:
            g1b.write('    for ccd in (\''+'b'+'\'): \n')
        else:
            g1b.write(line)
    g1b.close()

    g1r=open(r_apertures_out,'w')
    for line in apertures_text:
        if 'enter_run_here' in line:
            g1r.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
        elif 'for ccd in ' in line:
            g1r.write('    for ccd in (\''+'r'+'\'): \n')
        else:
            g1r.write(line)
    g1r.close()

    g1b=open(b_apertures_slurm_out,'w')
    g1b.write('#!/bin/bash \n')
    g1b.write('#SBATCH -c 1 \n')
    g1b.write('#SBATCH -t 0-07:59 \n')
    g1b.write('#SBATCH --partition=short \n')
    g1b.write('#SBATCH --mem-per-cpu=100000 \n')
    g1b.write('#SBATCH -o m2fs_b_apertures_'+m2fsrun0[i]+'.o \n')
    g1b.write('#SBATCH -e m2fs_b_apertures_'+m2fsrun0[i]+'.err \n')
    g1b.write(' \n')
    g1b.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/'+b_apertures_out+' \n')
    g1b.close()

    g1r=open(r_apertures_slurm_out,'w')
    g1r.write('#!/bin/bash \n')
    g1r.write('#SBATCH -c 1 \n')
    g1r.write('#SBATCH -t 0-07:59 \n')
    g1r.write('#SBATCH --partition=short \n')
    g1r.write('#SBATCH --mem-per-cpu=100000 \n')
    g1r.write('#SBATCH -o m2fs_r_apertures_'+m2fsrun0[i]+'.o \n')
    g1r.write('#SBATCH -e m2fs_r_apertures_'+m2fsrun0[i]+'.err \n')
    g1r.write(' \n')
    g1r.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/'+r_apertures_out+' \n')
    g1r.close()

    utdate=[]
    file1=[]
    file2=[]
    flatfile=[]
    tharfile=[]
    field_name=[]
    scifile=[]
    fibermap_file=[]
    filter_correction=[]

    with open(directory+m2fsrun0[i]+'_science_raw') as f:
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
        if 'filter_correction' in line:
            p2=line.split('filter_correction:')
            filter_correction.append(p2[1])
        else:
            filter_correction.append('none')
    utdate=np.array(utdate)
    file1=np.array(file1)
    file2=np.array(file2)
    flatfile=np.array(flatfile)
    tharfile=np.array(tharfile)
    field_name=np.array(field_name)
    scifile=np.array(scifile)
    fibermap_file=np.array(fibermap_file)
    
    for j in range(0,len(utdate)):
        science_raw_out=m2fsrun0[i]+'_'+str(j+1)+'_science_raw'
        b_apertures2_out='m2fs_b_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.py'
        b_apertures2_slurm_out='/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_b_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.slurm'
        r_apertures2_out='m2fs_r_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.py'
        r_apertures2_slurm_out='/nfs/nas-0-9/mgwalker.proj/scripts/m2fs_r_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.slurm'
        g1=open(science_raw_out,'w')
        g2b=open(b_apertures2_out,'w')
        g3b=open(b_apertures2_slurm_out,'w')
        g2r=open(r_apertures2_out,'w')
        g3r=open(r_apertures2_slurm_out,'w')
        g1.write(utdate[j]+' '+str(file1[j]).zfill(4)+' '+str(file2[j]).zfill(4)+' '+flatfile[j]+' '+tharfile[j]+' '+field_name[j]+' '+scifile[j]+' '+fibermap_file[j]+' \n')
        g1.close()
        
        for line in apertures_text:
            if 'enter_run_here' in line:
                g2b.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
            elif 'with open(directory+m2fsrun+\'_science_raw\') as f:' in line:
                g2b.write('with open(directory+m2fsrun+\''+'_'+str(j+1)+'_science_raw\') as f: \n')
            elif 'for ccd in ' in line:
                g2b.write('    for ccd in (\''+'b'+'\'): \n')
            else:
                g2b.write(line)
        g2b.close()

        for line in apertures_text:
            if 'enter_run_here' in line:
                g2r.write('m2fsrun=\''+m2fsrun0[i]+'\' \n')
            elif 'with open(directory+m2fsrun+\'_science_raw\') as f:' in line:
                g2r.write('with open(directory+m2fsrun+\''+'_'+str(j+1)+'_science_raw\') as f: \n')
            elif 'for ccd in ' in line:
                g2r.write('    for ccd in (\''+'r'+'\'): \n')
            else:
                g2r.write(line)
        g2r.close()

        g3b.write('#!/bin/bash \n')
        g3b.write('#SBATCH -c 1 \n')
        g3b.write('#SBATCH -t 0-07:59 \n')
        g3b.write('#SBATCH --partition=short \n')
        g3b.write('#SBATCH --mem-per-cpu=100000 \n')
        g3b.write('#SBATCH -o m2fs_b_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.o \n')
        g3b.write('#SBATCH -e m2fs_b_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.err \n')
        g3b.write(' \n')
        g3b.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/m2fs_b_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.py \n')
        g3b.close()

        g3r.write('#!/bin/bash \n')
        g3r.write('#SBATCH -c 1 \n')
        g3r.write('#SBATCH -t 0-07:59 \n')
        g3r.write('#SBATCH --partition=short \n')
        g3r.write('#SBATCH --mem-per-cpu=100000 \n')
        g3r.write('#SBATCH -o m2fs_r_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.o \n')
        g3r.write('#SBATCH -e m2fs_r_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.err \n')
        g3r.write(' \n')
        g3r.write('python /nfs/nas-0-9/mgwalker.proj/m2fs/m2fs_r_apertures_'+m2fsrun0[i]+'_'+str(j+1)+'.py \n')
        g3r.close()

