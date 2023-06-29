import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits 
from scipy import stats
import random
import mycode2
np.set_printoptions(threshold='nan')

#plt.rc('grid',alpha=0.7) # modify rcparams
#plt.rc('grid',color='white')
#plt.rc('grid',linestyle='-')
plt.rc('legend',frameon='False')
plt.rc('xtick.major',size=2)
plt.rc('ytick.major',size=2)
plt.rc('axes',axisbelow='True')
plt.rc('axes',grid='False')
plt.rc('axes',facecolor='white')
plt.rc('axes',linewidth=0.5)
plt.rc('xtick.major',width=0.5)
plt.rc('ytick.major',width=0.5)
plt.rc('xtick',direction='in')
plt.rc('ytick',direction='in')
plt.rc('xtick',labelsize='9')
plt.rc('ytick',labelsize='9')
plt.rc('text',usetex='True')
plt.rc('font',family='serif')
plt.rc('font',style='normal')
plt.rc('font',serif='Century')
plt.rc('savefig',format='eps')
plt.rc('lines',markeredgewidth=0)
plt.rc('lines',markersize=2)
            
tuc2_dmodulus=18.8
gru1_dmodulus=20.4

tuc2_data1='plates/tuc2/tuc2_ascii.txt'
tuc2_data2='age12fehm250.iso'
tuc2_data3='craptuc2.dat'
tuc2_data4='tuc2_ellipsehalf.res'

with open(tuc2_data2) as g: # read data file
    data=g.readlines()

tuc2_iso_g=[]
tuc2_iso_r=[]
tuc2_iso_i=[]

for line in data: # fill arrays
    p=line.split()
    tuc2_iso_g.append(float(p[6]))
    tuc2_iso_r.append(float(p[7]))
    tuc2_iso_i.append(float(p[8]))

tuc2_iso_g=np.array(tuc2_iso_g)
tuc2_iso_r=np.array(tuc2_iso_r)
tuc2_iso_i=np.array(tuc2_iso_i)

with open(tuc2_data4) as g: # read data file
    data=g.readlines()

tuc2_ellipse_x=[]
tuc2_ellipse_y=[]

for line in data: # fill arrays
    p=line.split()
    tuc2_ellipse_x.append(float(p[0]))
    tuc2_ellipse_y.append(float(p[1]))

tuc2_ellipse_x=np.array(tuc2_ellipse_x)
tuc2_ellipse_y=np.array(tuc2_ellipse_y)

with open(tuc2_data1) as f: # read data file
    data=f.readlines()

tuc2_v_radeg=[]
tuc2_v_decdeg=[]
tuc2_v_gmag=[]
tuc2_v_rmag=[]
tuc2_v_extg=[]
tuc2_v_extr=[]
tuc2_v_spreadg=[]
tuc2_v_spreadr=[]
tuc2_v_spreadi=[]

for line in data: # fill arrays
    p=line.split()
    tuc2_v_radeg.append(float(p[0]))
    tuc2_v_decdeg.append(float(p[1]))
    tuc2_v_gmag.append(float(p[2]))
    tuc2_v_rmag.append(float(p[3]))
    tuc2_v_extg.append(float(p[4]))
    tuc2_v_extr.append(float(p[5]))
    tuc2_v_spreadg.append(float(p[7]))
    tuc2_v_spreadr.append(float(p[8]))

tuc2_v_radeg=np.array(tuc2_v_radeg)
tuc2_v_decdeg=np.array(tuc2_v_decdeg)
tuc2_v_gmag=np.array(tuc2_v_gmag)
tuc2_v_rmag=np.array(tuc2_v_rmag)
tuc2_v_extg=np.array(tuc2_v_extg)
tuc2_v_extr=np.array(tuc2_v_extr)
tuc2_v_spreadg=np.array(tuc2_v_spreadg)
tuc2_v_spreadr=np.array(tuc2_v_spreadr)
tuc2_v_gmag=tuc2_v_gmag-tuc2_v_extg
tuc2_v_rmag=tuc2_v_rmag-tuc2_v_extr

tuc2_v_ra=tuc2_v_radeg*np.pi/180.
tuc2_v_dec=tuc2_v_decdeg*np.pi/180.
tuc2_v_radegcenter=0.73
tuc2_v_decdegcenter=-60.85
tuc2_v_racenter=np.zeros(len(tuc2_v_ra))+tuc2_v_radegcenter*np.pi/180.
tuc2_v_deccenter=np.zeros(len(tuc2_v_dec))+tuc2_v_decdegcenter*np.pi/180.
tuc2_v_xi,tuc2_v_eta=mycode2.etaxiarr(tuc2_v_ra,tuc2_v_dec,tuc2_v_racenter,tuc2_v_deccenter)
tuc2_v_r=np.sqrt(tuc2_v_xi**2+tuc2_v_eta**2)

tuc2_v_cand=np.zeros(len(tuc2_v_ra),dtype='int')
for i in range(0,len(tuc2_v_cand)):
    tuc2_dist=np.sqrt((tuc2_v_rmag[i]-(tuc2_iso_r+tuc2_dmodulus))**2+((tuc2_v_gmag[i]-tuc2_v_rmag[i])-(tuc2_iso_g-tuc2_iso_r))**2)
    if (min(tuc2_dist) < 0.2) & (tuc2_v_rmag[i] < 22.):
        tuc2_v_cand[i]=1

tuc2_v_cands=np.where(tuc2_v_cand == 1)

tuc2_data1='tuc2.dat'

with open(tuc2_data3) as f: # read data file
    data=f.readlines()[1:]
    
tuc2_radeg=[]
tuc2_decdeg=[]
tuc2_xi=[]
tuc2_eta=[]
tuc2_r=[]
tuc2_v=[]
tuc2_sigv=[]
tuc2_teff=[]
tuc2_sigteff=[]
tuc2_logg=[]
tuc2_siglogg=[]
tuc2_feh=[]
tuc2_sigfeh=[]
tuc2_gmag=[]
tuc2_rmag=[]
for line in data: # fill arrays
    p=line.split()
    tuc2_radeg.append(float(p[0]))
    tuc2_decdeg.append(float(p[1]))
    tuc2_xi.append(float(p[2]))
    tuc2_eta.append(float(p[3]))
    tuc2_r.append(float(p[4]))
    tuc2_v.append(float(p[5]))
    tuc2_sigv.append(float(p[6]))
    tuc2_teff.append(float(p[7]))
    tuc2_sigteff.append(float(p[8]))
    tuc2_logg.append(float(p[9]))
    tuc2_siglogg.append(float(p[10]))
    tuc2_feh.append(float(p[11]))
    tuc2_sigfeh.append(float(p[12]))
    tuc2_gmag.append(float(p[13]))
    tuc2_rmag.append(float(p[14]))
tuc2_radeg=np.array(tuc2_radeg)
tuc2_decdeg=np.array(tuc2_decdeg)
tuc2_xi=np.array(tuc2_xi)
tuc2_eta=np.array(tuc2_eta)
tuc2_r=np.array(tuc2_r)
tuc2_v=np.array(tuc2_v)
tuc2_sigv=np.array(tuc2_sigv)
tuc2_teff=np.array(tuc2_teff)
tuc2_sigteff=np.array(tuc2_sigteff)
tuc2_logg=np.array(tuc2_logg)
tuc2_siglogg=np.array(tuc2_siglogg)
tuc2_feh=np.array(tuc2_feh)
tuc2_sigfeh=np.array(tuc2_sigfeh)
tuc2_sigr=tuc2_r-tuc2_r
tuc2_gmag=np.array(tuc2_gmag)
tuc2_rmag=np.array(tuc2_rmag)
tuc2_pmem=tuc2_r-tuc2_r
#tuc2_keep=np.where((tuc2_sigv < 15.) & (tuc2_skewv > -1.) & (tuc2_skewv < 1.) & (tuc2_kurtv > -1.1) & (tuc2_kurtv < 1.1))
tuc2_keep=np.where((tuc2_sigv < 15.))# & (tuc2_skewv > -1.) & (tuc2_skewv < 1.) & (tuc2_kurtv > -1.1) & (tuc2_kurtv < 1.1))
tuc2_mem=np.where((tuc2_sigv < 15.) & (tuc2_pmem > 0.5))
tuc2_goodmem=np.where((tuc2_sigv < 15.) & (tuc2_pmem > 0.9))
tuc2_nonmem=np.where((tuc2_sigv < 15.) & (tuc2_pmem < 0.5))
tuc2_goodnonmem=np.where((tuc2_sigv < 15.) & (tuc2_pmem < 0.5))
tuc2_bad=np.where(tuc2_sigv >= 15.)

tuc2_vmem=[-150.,-110.]
tuc2_loggmem=[0.,4.0]
tuc2_fehmem=[-5.,0.]
tuc2_teffmem=[4000,8000]
tuc2_rmem=[0.,12.]
tuc2_magmem=[21,16]
tuc2_colmem=[-0.2,1]

tuc2_center=np.where((tuc2_v_r < 15.))

gs=plt.GridSpec(10,10) # define multi-panel plot
gs.update(wspace=0,hspace=0) # specify inter-panel spacing
fig=plt.figure(figsize=(6,6)) # define plot size
#fig.set_
ax0_0=fig.add_subplot(gs[0:5,0:4])
ax0_1=fig.add_subplot(gs[0:5,4:8])

ax0_0.set_xlabel(r'$g-r$ [mag]',fontsize=12,rotation=0,labelpad=5)
ax0_0.set_ylabel(r'$r$ [mag]',fontsize=12,rotation=90,labelpad=5)
ax0_0.set_xlim([-0.5,1.5])
ax0_0.set_ylim([22.5,15.5])
ax0_0.set_xscale(u'linear')
ax0_0.set_yscale(u'linear')
ax0_0.set_yticks([22,21,20,19,18,17,16])
#ax0_0.set_xticklabels(vticks,rotation=0)
ax0_0.set_xticks([-0.5,0,0.5,1])
#ax0_0.set_yticklabels([4000]+teffticks,rotation=0)
#ax0_0.errorbar(v[tuc2_keep],teff[tuc2_keep],xerr=sigv[tuc2_keep],yerr=sigteff[tuc2_keep],elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
#ax0_0.scatter(tuc2_v_gmag[tuc2_center]-tuc2_v_rmag[tuc2_center],tuc2_v_rmag[tuc2_center],s=2,lw=0,edgecolor='none',alpha=0.99,marker='x',color='0.5',rasterized=True)
ax0_0.scatter(tuc2_v_gmag[tuc2_center]-tuc2_v_rmag[tuc2_center],tuc2_v_rmag[tuc2_center],s=2,lw=0,edgecolor='none',alpha=0.99,marker='o',color='0.5',rasterized=True)
ax0_0.scatter(tuc2_gmag[tuc2_bad]-tuc2_rmag[tuc2_bad],tuc2_rmag[tuc2_bad],s=10,lw=1,edgecolor='k',alpha=0.99,marker='x',color='k',rasterized=True)
ax0_0.scatter(tuc2_gmag[tuc2_nonmem]-tuc2_rmag[tuc2_nonmem],tuc2_rmag[tuc2_nonmem],s=15,lw=1,edgecolor='k',alpha=0.99,marker='o',facecolor='none',rasterized=True)
ax0_0.scatter(tuc2_gmag[tuc2_goodnonmem]-tuc2_rmag[tuc2_goodnonmem],tuc2_rmag[tuc2_goodnonmem],s=15,lw=1,edgecolor='k',alpha=0.99,marker='o',color='k',rasterized=True)
ax0_0.scatter(tuc2_gmag[tuc2_mem]-tuc2_rmag[tuc2_mem],tuc2_rmag[tuc2_mem],s=15,lw=1,edgecolor='r',alpha=0.99,marker='o',facecolor='none',rasterized=True)
ax0_0.scatter(tuc2_gmag[tuc2_goodmem]-tuc2_rmag[tuc2_goodmem],tuc2_rmag[tuc2_goodmem],s=15,lw=1,edgecolor='r',alpha=0.99,marker='o',color='r',rasterized=True)
ax0_0.text(-0.4,16.25,'Tuc 4',fontsize=12)


ax0_0.plot(tuc2_iso_g-tuc2_iso_r,tuc2_iso_r+tuc2_dmodulus,lw=1,color='r')
#ax0_0.plot(vteffx2,vteffy2,color='b',linewidth=0.25)


ax0_1.yaxis.set_major_formatter(plt.NullFormatter())
ax0_1.set_xlabel(r'$g-r$ [mag]',fontsize=12,rotation=0,labelpad=5)
#ax0_1.set_ylabel(r'$r$ [mag]',fontsize=8,rotation=90,labelpad=5)
ax0_1.set_xlim([-0.5,1.5])
ax0_1.set_ylim([22.5,15.5])
ax0_1.set_xscale(u'linear')
ax0_1.set_yscale(u'linear')
ax0_1.set_yticks([22,21,20,19,18,17,16])
ax0_1.set_yticklabels([])
#ax0_1.set_xticklabels(vticks,rotation=0)
ax0_1.set_xticks([0,0.5,1,1.5])
#ax0_1.set_yticklabels([4000]+teffticks,rotation=0)
#ax0_1.errorbar(v[gru1_keep],teff[gru1_keep],xerr=sigv[gru1_keep],yerr=sigteff[gru1_keep],elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
#ax0_1.plot(vteffx2,vteffy2,color='b',linewidth=0.25)


plotfilename='craptuc2_cmd.pdf'
plt.savefig(plotfilename,dpi=400)
plt.show()
plt.close()
