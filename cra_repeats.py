import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits 
from scipy import stats
import random
import matplotlib.mlab as mlab
#import mycode2
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
plt.rc('xtick',labelsize='6')
plt.rc('ytick',labelsize='6')
plt.rc('text',usetex='True')
plt.rc('font',family='serif')
plt.rc('font',style='normal')
plt.rc('font',serif='Century')
plt.rc('savefig',format='eps')
plt.rc('lines',markeredgewidth=0)
plt.rc('lines',markersize=2)
            
data1='cra_repeats.dat'
data2='cra_repeats2.dat'
data3='cra_repeats3.dat'
data4='cra_repeats4.dat'

with open(data1) as f: # read data file
    data=f.readlines()

v=[]
sigv=[]
vmean=[]
sigvmean=[]
teff=[]
sigteff=[]
teffmean=[]
sigteffmean=[]
logg=[]
siglogg=[]
loggmean=[]
sigloggmean=[]
feh=[]
sigfeh=[]
fehmean=[]
sigfehmean=[]

for line in data: # fill arrays
    p=line.split()
    v.append(float(p[0]))
    sigv.append(float(p[1]))
    vmean.append(float(p[2]))
    sigvmean.append(float(p[3]))
    teff.append(float(p[4]))
    sigteff.append(float(p[5]))
    teffmean.append(float(p[6]))
    sigteffmean.append(float(p[7]))
    logg.append(float(p[8]))
    siglogg.append(float(p[9]))
    loggmean.append(float(p[10]))
    sigloggmean.append(float(p[11]))
    feh.append(float(p[12]))
    sigfeh.append(float(p[13]))
    fehmean.append(float(p[14]))
    sigfehmean.append(float(p[15]))

v=np.array(v)
sigv=np.array(sigv)
vmean=np.array(vmean)
sigvmean=np.array(sigvmean)
teff=np.array(teff)
sigteff=np.array(sigteff)
teffmean=np.array(teffmean)
sigteffmean=np.array(sigteffmean)
logg=np.array(logg)
siglogg=np.array(siglogg)
loggmean=np.array(loggmean)
sigloggmean=np.array(sigloggmean)
feh=np.array(feh)
sigfeh=np.array(sigfeh)
fehmean=np.array(fehmean)
sigfehmean=np.array(sigfehmean)

with open(data2) as f: # read data file
    data=f.readlines()

v1=[]
sigv1=[]
v2=[]
sigv2=[]
teff1=[]
sigteff1=[]
teff2=[]
sigteff2=[]
logg1=[]
siglogg1=[]
logg2=[]
siglogg2=[]
feh1=[]
sigfeh1=[]
feh2=[]
sigfeh2=[]

for line in data: # fill arrays
    p=line.split()
    v1.append(float(p[0]))
    sigv1.append(float(p[1]))
    v2.append(float(p[2]))
    sigv2.append(float(p[3]))
    teff1.append(float(p[4]))
    sigteff1.append(float(p[5]))
    teff2.append(float(p[6]))
    sigteff2.append(float(p[7]))
    logg1.append(float(p[8]))
    siglogg1.append(float(p[9]))
    logg2.append(float(p[10]))
    siglogg2.append(float(p[11]))
    feh1.append(float(p[12]))
    sigfeh1.append(float(p[13]))
    feh2.append(float(p[14]))
    sigfeh2.append(float(p[15]))

v1=np.array(v1)
sigv1=np.array(sigv1)
v2=np.array(v2)
sigv2=np.array(sigv2)
teff1=np.array(teff1)
sigteff1=np.array(sigteff1)
teff2=np.array(teff2)
sigteff2=np.array(sigteff2)
logg1=np.array(logg1)
siglogg1=np.array(siglogg1)
logg2=np.array(logg2)
siglogg2=np.array(siglogg2)
feh1=np.array(feh1)
sigfeh1=np.array(sigfeh1)
feh2=np.array(feh2)
sigfeh2=np.array(sigfeh2)

vdiff=(v1-v2)/np.sqrt(sigv1**2+sigv2**2)#(v-vmean)/sigv
teffdiff=(teff1-teff2)/np.sqrt(sigteff1**2+sigteff2**2)#(teff-teffmean)/sigteff
loggdiff=(logg1-logg2)/np.sqrt(siglogg1**2+siglogg2**2)#(logg-loggmean)/siglogg
fehdiff=(feh1-feh2)/np.sqrt(sigfeh1**2+sigfeh2**2)

with open(data4) as f4: # read data file
    data=f4.readlines()

fake=[]
for line in data: # fill arrays
    p=line.split()
    fake.append(float(p[0]))

fake=np.array(fake)

gs=plt.GridSpec(4,4) # define multi-panel plot
gs.update(wspace=0.5,hspace=0.5) # specify inter-panel spacing
fig=plt.figure(figsize=(6,6)) # define plot size
#fig.set_

ax0_0=fig.add_subplot(gs[0,0])
ax0_1=fig.add_subplot(gs[0,1])
ax0_2=fig.add_subplot(gs[0,2])
ax0_3=fig.add_subplot(gs[0,3])
ax1_0=fig.add_subplot(gs[1,0])
ax1_1=fig.add_subplot(gs[1,1])
ax1_2=fig.add_subplot(gs[1,2])
ax1_3=fig.add_subplot(gs[1,3])

   
vlim=[-25.,199.]
tefflim=[4001.,7999.] 
revtefflim=[8000.,4000.] 
logglim=[0.01,4.99]
revlogglim=[5,0]
fehlim=[-4.99,0.99]
rlim=[0.,20.]

vlabelx1=[r'$v_{\rm los,first}$ [km/s]']
vlabely1=[r'$v_{\rm los,later}$ [km/s]']
tefflabelx1=[r'$T_{\rm eff,1}$ [$10^3 \mathrm{K}$]']
tefflabely1=[r'$T_{\rm eff,2}-T_{\rm eff,1}$ [$10^3 \mathrm{K}$]']
logglabelx1=[r'$\log g_1$']
logglabely1=[r'$\log g_2 - \log g_1$']
fehlabelx1=[r'$[\mathrm{Fe/H}]_1$']
fehlabely1=[r'$[\mathrm{Fe/H}]_2$-[Fe/H]$_1$']

ax0_0.set_xlabel(r"$v_{\rm los,1}$ [km/s]",fontsize=8,rotation=0)
ax0_0.set_ylabel(r"$v_{\rm los,1}-v_{\rm los,2}$ [km/s]",fontsize=8,rotation=90,labelpad=5)
ax0_0.set_xlim([-115,400])
ax0_0.set_ylim([-10,10])
ax0_0.set_xscale(u'linear')
ax0_0.set_yscale(u'linear')
#ax0_0.set_xticks([0,100,200,300])
#ax0_0.set_yticks([0,100,200,300])
ax0_0.set_xticks([0,200,400])
ax0_0.set_yticks([-10,-5,0,5,10])
ax0_0.plot([-1000,1000],[000,000],linestyle=':',color='k',linewidth=0.25)
ax0_0.errorbar(v1,v1-v2,xerr=sigv1,yerr=np.sqrt(sigv1**2+sigv2**2),elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)

ax0_1.set_xlabel(tefflabelx1[0],fontsize=8,rotation=0)
ax0_1.set_ylabel(tefflabely1[0],fontsize=8,rotation=90,labelpad=1)
ax0_1.set_xlim([4,8])
ax0_1.set_ylim([-2,2])
ax0_1.set_xscale(u'linear')
ax0_1.set_yscale(u'linear')
ax0_1.set_xticks([4,5,6,7,8])
ax0_1.set_yticks([-2,-1,0,1,2])
ax0_1.plot([-1000,1000],[-000,000],linestyle=':',color='k',linewidth=0.25)
ax0_1.errorbar(teff1/1000.,teff2/1000.-teff1/1000.,xerr=sigteff1/1000.,yerr=np.sqrt((sigteff2/1000.)**2+(sigteff1/1000.)**2),elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)

ax0_2.set_xlabel(logglabelx1[0],fontsize=8,rotation=0)
ax0_2.set_ylabel(logglabely1[0],fontsize=8,rotation=90,labelpad=1)
ax0_2.set_xlim([0,5])
ax0_2.set_ylim([-2,2])
ax0_2.set_xscale(u'linear')
ax0_2.set_yscale(u'linear')
ax0_2.set_xticks([0,1,2,3,4,5])
ax0_2.set_yticks([-2,-1,0,1,2])
ax0_2.plot([-1000,1000],[-000,000],linestyle=':',color='k',linewidth=0.25)
ax0_2.errorbar(logg1,logg2-logg1,xerr=siglogg1,yerr=np.sqrt(siglogg2**2+siglogg1**2),elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)


ax0_3.set_xlabel(fehlabelx1[0],fontsize=8,rotation=0)
ax0_3.set_ylabel(fehlabely1[0],fontsize=8,rotation=90,labelpad=1)
ax0_3.set_xlim([-5,1])
ax0_3.set_ylim([-2,2])
ax0_3.set_xscale(u'linear')
ax0_3.set_yscale(u'linear')
ax0_3.set_xticks([-5,-3,-1,1])
ax0_3.set_yticks([-2,-1,0,1,2])
ax0_3.plot([-1000,1000],[-000,000],linestyle=':',color='k',linewidth=0.25)
ax0_3.errorbar(feh1,feh2-feh1,xerr=sigfeh1,yerr=np.sqrt(sigfeh2**2+sigfeh1**2),elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)


x=np.linspace(-5,5,100)
mean=0.
sigma=1.

ax1_0.set_xlabel(r'$\frac{v_{\rm los,1}-v_{\rm los,2}}{(\delta^2_{v_{\rm los,1}}+\delta^2_{v_{\rm los,2}})^{1/2}}$')
ax1_0.set_ylabel('N',rotation=90,labelpad=10)
ax1_0.set_xlim([-5,5])
#ax1_0.set_ylim([0,10])
ax1_0.set_xscale(u'linear')
ax1_0.set_yscale(u'linear')
ax1_0.set_xticks([-4,-2,0,2,4])
#ax1_0.yaxis.set_major_formatter(plt.NullFormatter())
#ax1_0.hist(fake,bins=15,range=[-5,5],normed=False,histtype='stepfilled',align='mid',rwidth=0,orientation='vertical',color='0.5',linewidth=0,alpha=0.5)
ax1_0.plot(x,mlab.normpdf(x,mean,sigma)*float(len(vdiff))*0.67,color='0.15',linewidth=0.65,zorder=0,linestyle='-',alpha=0.5)
ax1_0.hist(vdiff,bins=15,range=[-5,5],normed=False,histtype='step',align='mid',rwidth=0,orientation='vertical',color='k',linewidth=0.5)

ax1_1.set_xlabel(r'$\frac{T_{\rm eff,1}-T_{\rm eff,2}}{(\delta^2_{T_{\rm eff,1}}+\delta^2_{\rm eff,2})^{1/2}}$')
ax1_1.set_xlim([-5,5])
#ax1_1.set_ylim([0,10])
ax1_1.set_xscale(u'linear')
ax1_1.set_yscale(u'linear')
ax1_1.set_xticks([-4,-2,0,2,4])
#ax1_1.yaxis.set_major_formatter(plt.NullFormatter())
ax1_1.plot(x,mlab.normpdf(x,mean,sigma)*float(len(teffdiff))*0.67,color='0.15',linewidth=0.65,zorder=0,linestyle='-',alpha=0.5)
#ax1_1.hist(fake,bins=15,range=[-5,5],normed=False,histtype='stepfilled',align='mid',rwidth=0,orientation='vertical',color='0.5',linewidth=0,alpha=0.5)
ax1_1.hist(teffdiff,bins=15,range=[-5,5],normed=False,histtype='step',align='mid',rwidth=0,orientation='vertical',color='k',linewidth=0.5)

ax1_2.set_xlabel(r'$\frac{\log g_1-\log g_2}{(\delta^2_{\log g,1}+\delta^2_{\log g,2})^{1/2}}$')
ax1_2.set_xlim([-5,5])
#ax1_2.set_ylim([0,10])
ax1_2.set_xscale(u'linear')
ax1_2.set_yscale(u'linear')
ax1_2.set_xticks([-4,-2,0,2,4])
#ax1_2.yaxis.set_major_formatter(plt.NullFormatter())
ax1_2.plot(x,mlab.normpdf(x,mean,sigma)*float(len(loggdiff))*0.67,color='0.15',linewidth=0.65,zorder=0,linestyle='-',alpha=0.5)
#ax1_2.hist(fake,bins=15,range=[-5,5],normed=False,histtype='stepfilled',align='mid',rwidth=0,orientation='vertical',color='0.5',linewidth=0,alpha=0.5)
ax1_2.hist(loggdiff,bins=15,range=[-5,5],normed=False,histtype='step',align='mid',rwidth=0,orientation='vertical',color='k',linewidth=0.5)

ax1_3.set_xlabel(r'$\frac{[\mathrm{Fe/H}]-\langle [\mathrm{Fe/H}]\rangle}{\sigma_{[\mathrm{Fe/H}]}}$')
ax1_3.set_xlim([-5,5])
#ax1_3.set_ylim([0,10])
ax1_3.set_xscale(u'linear')
ax1_3.set_yscale(u'linear')
ax1_3.set_xticks([-4,-2,0,2,4])
#ax1_3.yaxis.set_major_formatter(plt.NullFormatter())
ax1_3.plot(x,mlab.normpdf(x,mean,sigma)*float(len(fehdiff))*0.67,color='0.15',linewidth=0.65,zorder=0,linestyle='-',alpha=0.5)
#ax1_3.hist(fake,bins=15,range=[-5,5],normed=False,histtype='stepfilled',align='mid',rwidth=0,orientation='vertical',color='0.5',linewidth=0,alpha=0.5)
ax1_3.hist(fehdiff,bins=15,range=[-5,5],normed=False,histtype='step',align='mid',rwidth=0,orientation='vertical',color='k',linewidth=0.5)

plotfilename='cra_repeats.pdf'
plt.savefig(plotfilename,dpi=400)
plt.show()
plt.close()
