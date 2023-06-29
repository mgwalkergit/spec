import numpy as np
import astropy
from astropy import units
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import CCDData
import ccdproc
import astropy.units as u
from astropy.modeling import models
from ccdproc import Combiner
import os
import mycode

directory='/nfs/nas-0-9/mgwalker.proj/m2fs/'
m2fsrun='nov18'
datadir='/nfs/nas-0-9/mgwalker.proj/m2fs/nov18/m2fs.astro.lsa.umich.edu/data/NovDec2018/'

utdate=[]
file1=[]
file2=[]
#with open(directory+m2fsrun+'_fibermap_raw') as f:
#    data=f.readlines()[0:]
#for line in data:
#    p=line.split()
#    utdate.append(str(p[0]))
#    file1.append(long(p[1]))
#    file2.append(long(p[2]))
#with open(directory+m2fsrun+'_twilight_raw') as f:
#    data=f.readlines()[0:]
#for line in data:
#    p=line.split()
#    utdate.append(str(p[0]))
#    file1.append(long(p[1]))
#    file2.append(long(p[2]))
with open(directory+m2fsrun+'_crapscience_raw') as f:
    data=f.readlines()[0:]
for line in data:
    p=line.split()
    utdate.append(str(p[0]))
    file1.append(long(p[1]))
    file2.append(long(p[2]))

utdate=np.array(utdate)
file1=np.array(file1)
file2=np.array(file2)

for ccd in ('b','r'):
    for chip in ('c1','c2','c3','c4'):
        master_bias=astropy.nddata.CCDData.read(directory+m2fsrun+'_'+ccd+'_'+chip+'_master_bias.fits')
        master_dark=astropy.nddata.CCDData.read(directory+m2fsrun+'_'+ccd+'_'+chip+'_master_dark.fits')
        for i in range(0,len(utdate)):
            for j in range(file1[i],file2[i]+1):
                filename=datadir+utdate[i]+'/'+ccd+str(j).zfill(4)+chip+'.fits'
                out=datadir+utdate[i]+'/'+ccd+str(j).zfill(4)+chip+'_debias_dedark.fits'
                data=astropy.nddata.CCDData.read(filename,unit=u.adu)#header is in data.meta
                print(filename,data.header['object'],data.header['binning'])
                data_with_deviation=ccdproc.create_deviation(data,gain=data.meta['egain']*u.electron/u.adu,readnoise=data.meta['enoise']*u.electron)
                gain_corrected=ccdproc.gain_correct(data_with_deviation,data_with_deviation.meta['egain']*u.electron/u.adu,add_keyword={'gain_corr':'Done'})
                cr_cleaned=ccdproc.cosmicray_lacosmic(gain_corrected,sigclip=5)
                oscan_subtracted=ccdproc.subtract_overscan(cr_cleaned,overscan=cr_cleaned[:,1024:],overscan_axis=1,model=models.Polynomial1D(3),add_keyword={'oscan_corr':'Done'})
                trimmed1=ccdproc.trim_image(oscan_subtracted[:,:1024],add_keyword={'trim1':'Done'})
                trimmed2=ccdproc.trim_image(trimmed1[:1028,:1024],add_keyword={'trim2':'Done'})
                debiased0=ccdproc.subtract_bias(trimmed2,master_bias)
                dedark0=ccdproc.subtract_dark(debiased0,master_dark,exposure_time='exptime',exposure_unit=u.second,scale=True,add_keyword={'dark_corr':'Done'})
                os.popen('rm '+out)
                dedark0.write(out)

                if chip=='c1':
                    c1_process=dedark0
                if chip=='c2':
                    c2_process=dedark0
                if chip=='c3':
                    c3_process=dedark0
                if chip=='c4':
                    c4_process=dedark0
    left=astropy.nddata.CCDData(np.concatenate((c1_process,np.flipud(c4_process)),axis=0),unit=u.electron)
    right=astropy.nddata.CCDData(np.concatenate((np.fliplr(c2_process),np.fliplr(np.flipud(c3_process))),axis=0),unit=u.electron)
    stitched=astropy.nddata.CCDData(np.concatenate((left,right),axis=1),unit=u.electron)
    np.pause()
                
