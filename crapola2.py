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

with open('m2fshires_jul15.dat') as f:
    data=f.readlines()[1:]
snratio=[]
rmag=[]
for line in data:
    p=line.split()
    snratio.append(float(p[29]))
    rmag.append(float(p[24]))
snratio=np.array(snratio)
rmag=np.array(rmag)

keep=np.where(rmag>10.)[0]
plt.scatter(rmag[keep],np.log10(snratio[keep]),s=5,color='r',label='HiRes',marker='o')
plt.legend(loc=1)
plt.xlabel('r mag')
plt.ylabel(r'median $\log_{10}[S/N/\mathrm{pixel}]$')
plt.savefig('m2fs_sn.pdf',dpi=200)
plt.show()
plt.close()
