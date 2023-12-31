import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits 
from scipy import stats
import random
#import mycode2
np.set_printoptions(threshold='nan')

#plt.rc('grid',alpha=0.7) # modify rcparams
#plt.rc('grid',color='white')
#plt.rc('grid',linestyle='-')
plt.rc('legend',frameon='False')
plt.rc('xtick.major',size=2)
plt.rc('ytick.major',size=2)
plt.rc('xtick.major',width=0.25)
plt.rc('ytick.major',width=0.25)
plt.rc('xtick.minor',size=2)
plt.rc('ytick.minor',size=2)
plt.rc('xtick.minor',width=0.25)
plt.rc('ytick.minor',width=0.25)
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
            
data1='m2fsmedres_jul15.dat'
with open(data1) as f: # read data file
    data=f.readlines()
    
radeg=[]
decdeg=[]
xi=[]
eta=[]
r=[]
hjd=[]
v=[]
sigv=[]
skewv=[]
kurtv=[]
teff=[]
sigteff=[]
skewteff=[]
kurtteff=[]
logg=[]
siglogg=[]
skewlogg=[]
kurtlogg=[]
feh=[]
sigfeh=[]
skewfeh=[]
kurtfeh=[]
pmem=[]
snratio=[]
rmag=[]

for line in data: # fill arrays
    p=line.split()
    radeg.append(float(p[0]))
    decdeg.append(float(p[1]))
    xi.append(float(p[2]))
    eta.append(float(p[3]))
    r.append(float(p[4]))
    hjd.append(float(p[5]))
    v.append(float(p[6]))
    sigv.append(float(p[7]))
    skewv.append(float(p[8]))
    kurtv.append(float(p[9]))
    teff.append(float(p[10]))
    sigteff.append(float(p[11]))
    skewteff.append(float(p[12]))
    kurtteff.append(float(p[13]))
    logg.append(float(p[14]))
    siglogg.append(float(p[15]))
    skewlogg.append(float(p[16]))
    kurtlogg.append(float(p[17]))
    feh.append(float(p[18]))
    sigfeh.append(float(p[19]))
    skewfeh.append(float(p[20]))
    kurtfeh.append(float(p[21]))
    pmem.append(float(p[22]))
    snratio.append(float(p[29]))
    rmag.append(float(p[24]))

radeg=np.array(radeg)
decdeg=np.array(decdeg)
xi=np.array(xi)
eta=np.array(eta)
r=np.array(r)
hjd=np.array(hjd)
v=np.array(v)
sigv=np.array(sigv)
skewv=np.array(skewv)
kurtv=np.array(kurtv)
teff=np.array(teff)
sigteff=np.array(sigteff)
skewteff=np.array(skewteff)
kurtteff=np.array(kurtteff)
logg=np.array(logg)
siglogg=np.array(siglogg)
skewlogg=np.array(skewlogg)
kurtlogg=np.array(kurtlogg)
feh=np.array(feh)
sigfeh=np.array(sigfeh)
skewfeh=np.array(skewfeh)
kurtfeh=np.array(kurtfeh)
pmem=np.array(pmem)
snratio=np.array(snratio)
sigr=r-r
rmag=np.array(rmag)

keep=np.where((sigv < 10.) & (skewv >= -1.) & (skewv <= 1.) & (kurtv >= -1.) & (kurtv <= 1.))
mem=np.where((sigv < 10.) & (skewv > -1.) & (skewv < 1.) & (kurtv > -1.) & (kurtv < 1.) & (pmem > 0.9))
mem2=np.where((sigv < 10.) & (skewv > -1.) & (skewv < 1.) & (kurtv > -1.) & (kurtv < 1.) & (pmem > 0.5) & (pmem < 0.8))
mem3=np.where((sigv < 10.) & (skewv > -1.) & (skewv < 1.) & (kurtv > -1.) & (kurtv < 1.) & (pmem > 0.5))

vmem=[55.,80.]
vmem2=[50.,85.]
loggmem=[0.25,4.5]
loggmem2=[0.05,4.7]
fehmem=[-4,-1.1]
fehmem2=[-4.2,-1.]
teffmem=[4250,7000]
teffmem2=[4200,7050]
rmem=[0.,12.]
rmem2=[0.,15]


gs=plt.GridSpec(20,20) # define multi-panel plot
gs2=plt.GridSpec(5,5) # define multi-panel plot
gs.update(wspace=0,hspace=0) # specify inter-panel spacing
gs2.update(wspace=1,hspace=1) # specify inter-panel spacing
fig=plt.figure(figsize=(6,6)) # define plot size
#fig.set_

ax1_1=fig.add_subplot(gs[7:9,0:2])
ax1_2=fig.add_subplot(gs[5:7,0:2])
ax1_3=fig.add_subplot(gs[3:5,0:2])
ax1_4=fig.add_subplot(gs[1:3,0:2])
ax11_1=fig.add_subplot(gs[7:9,2:3])   
ax11_2=fig.add_subplot(gs[5:7,2:3])   
ax11_3=fig.add_subplot(gs[3:5,2:3])   

ax2_1=fig.add_subplot(gs[7:9,4:6])
ax2_2=fig.add_subplot(gs[5:7,4:6])
ax2_3=fig.add_subplot(gs[3:5,4:6])
#ax2_4=fig.add_subplot(gs[1:3,4:6])
ax22_1=fig.add_subplot(gs[7:9,6:7])   
ax22_2=fig.add_subplot(gs[5:7,6:7])   
ax22_3=fig.add_subplot(gs[3:5,6:7])   

ax3_1=fig.add_subplot(gs[7:9,8:10])
ax3_2=fig.add_subplot(gs[5:7,8:10])
ax3_3=fig.add_subplot(gs[3:5,8:10])
#ax3_4=fig.add_subplot(gs[1:3,8:10])
ax33_1=fig.add_subplot(gs[7:9,10:11])   
ax33_2=fig.add_subplot(gs[5:7,10:11])   
ax33_3=fig.add_subplot(gs[3:5,10:11])   

ax4_4=fig.add_subplot(gs[1:3,12:14])

ax4_1=fig.add_subplot(gs[7:9,12:14])
ax4_2=fig.add_subplot(gs[5:7,12:14])
ax4_3=fig.add_subplot(gs[3:5,12:14])
#ax4_4=fig.add_subplot(gs[1:3,12:14])
ax44_1=fig.add_subplot(gs[7:9,14:15])   
ax44_2=fig.add_subplot(gs[5:7,14:15])   
ax44_3=fig.add_subplot(gs[3:5,14:15])   


ax1_1.set_xlabel('S/N',fontsize=8,rotation=0)
ax1_1.set_ylabel(r'variance$^{1/2}$',fontsize=8,rotation=90,labelpad=5)
ax1_1.set_xlim([0.1,200])
ax1_1.set_ylim([0.1,900])
ax1_1.set_xscale(u'log')
ax1_1.set_yscale(u'log')
ax1_1.set_xticks([0.1,1,10,100])
ax1_1.set_xticklabels([0.1,1,10,100],rotation=90)
ax1_1.set_yticks([0.1,1,10,100])
ax1_1.set_yticklabels([0.1,1,10,100],rotation=0)
ax1_1.scatter(snratio,sigv,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax1_1.scatter(snratio[keep],sigv[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax1_1.text(0.15,0.25,'[km/s]',fontsize=4.,color='k',alpha=1)

#ax1_2.set_xlabel('S/N',fontsize=8,rotation=0)
ax1_2.set_ylabel(r'skewness',fontsize=8,rotation=90,labelpad=5)
ax1_2.set_xlim([0.1,200])
ax1_2.set_ylim([-9,9])
ax1_2.set_xscale(u'log')
ax1_2.set_yscale(u'linear')
ax1_2.set_xticks([1,10,100])
#ax1_2.set_xticklabels([1,10,100],rotation=0)
ax1_2.xaxis.set_major_formatter(plt.NullFormatter())
ax1_2.set_yticks([-5,0,5])
ax1_2.set_yticklabels([-5,0,5],rotation=0)
ax1_2.scatter(snratio,skewv,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax1_2.scatter(snratio[keep],skewv[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)

#ax1_3.set_xlabel('S/N',fontsize=8,rotation=0)
ax1_3.set_ylabel(r'kurtosis',fontsize=8,rotation=90,labelpad=5)
ax1_3.set_xlim([0.1,200])
ax1_3.set_ylim([1,200])
ax1_3.set_xscale(u'log')
ax1_3.set_yscale(u'log')
ax1_3.set_xticks([1,10,100])
#ax1_3.set_xticklabels([1,10,100],rotation=0)
ax1_3.xaxis.set_major_formatter(plt.NullFormatter())
ax1_3.set_yticks([1,10,100])
ax1_3.set_yticklabels([1,10,100],rotation=0)
ax1_3.scatter(snratio,3.+kurtv,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax1_3.scatter(snratio[keep],3.+kurtv[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax1_3.text(100,40,r'$v_{\rm los}$',fontsize=7,color='k',alpha=1,horizontalalignment='right')

#ax1_4.set[keep]_xlabel[keep]('S/N',fontsize=8,rotation=0)
ax1_4.set_ylabel(r'r [mag]',fontsize=8,rotation=90,labelpad=5)
ax1_4.set_xlim([0.1,200])
ax1_4.set_ylim([23,19])
ax1_4.set_xscale(u'log')
ax1_4.set_yscale(u'linear')
#ax1_4.set_xticks([1,10,100])
#ax1_4.set_xticklabels([1,10,100],rotation=0)
ax1_4.xaxis.set_major_formatter(plt.NullFormatter())
ax1_4.set_yticks([22,21,20,19])
#ax1_4.set_yticklabels([10,20,30,40,50],rotation=0)
ax1_4.scatter(snratio,rmag,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax1_4.scatter(snratio[keep],rmag[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)


#ax4_4.set[keep]_xlabel[keep]('S/N',fontsize=8,rotation=0)
ax4_4.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
ax4_4.set_xlim(np.log10([0.1,200]))
ax4_4.set_ylim([0,50])
ax4_4.set_xscale(u'linear')
ax4_4.set_yscale(u'linear')
#ax4_4.set_xticks([1,10,100])
#ax4_4.set_xticklabels([1,10,100],rotation=0)
ax4_4.xaxis.set_major_formatter(plt.NullFormatter())
ax4_4.set_yticks([10,20,30,40,50])
#ax4_4.set_yticklabels([10,20,30,40,50],rotation=0)
ax4_4.hist(np.log10(snratio),bins=25,range=np.log10([0.1,200]),normed=False,align='mid',histtype='step',rwidth=0,orientation='vertical',color='k',linewidth=0.25,alpha=1)

#ax11_1.set_xlabel('N',fontsize=8,rotation=0)
#ax11_1.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax11_1.set_xlim([0,100])
ax11_1.set_ylim(np.log10([0.1,900]))
ax11_1.set_xscale(u'linear')
ax11_1.set_yscale(u'linear')
#ax11_1.set_xticks([30,60,90])
#ax11_1.set_xticklabels([30,60,90],rotation=90)
ax11_1.yaxis.set_major_formatter(plt.NullFormatter())
ax11_1.xaxis.set_major_formatter(plt.NullFormatter())
ax11_1.set_yticks(np.log10([0.1,1,10,100]))
#ax11_1.set_yticklabels([25,50,75],rotation=0)
ax11_1.hist(np.log10(sigv),bins=25,range=np.log10([0.1,900]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax11_2.set_xlabel('S/N',fontsize=8,rotation=0)
#ax11_2.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax11_2.set_xlim([0,150])
ax11_2.set_ylim([-9,9])
ax11_2.set_xscale(u'linear')
ax11_2.set_yscale(u'linear')
#ax11_2.set_xticks([1,10,100])
#ax11_2.set_xticklabels([1,10,100],rotation=0)
ax11_2.yaxis.set_major_formatter(plt.NullFormatter())
ax11_2.xaxis.set_major_formatter(plt.NullFormatter())
ax11_2.set_yticks([-5,0,5])
#ax11_2.set_yticklabels([25,50,75],rotation=0)
ax11_2.hist(skewv,bins=25,range=[-9,9],normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax11_3.set_xlabel('S/N',fontsize=8,rotation=0)
#ax11_3.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax11_3.set_xlim([0,150])
ax11_3.set_ylim(np.log10([1.,200.]))
ax11_3.set_xscale(u'linear')
ax11_3.set_yscale(u'linear')
#ax11_3.set_xticks([1,10,100])
#ax11_3.set_xticklabels([1,10,100],rotation=0)
ax11_3.yaxis.set_major_formatter(plt.NullFormatter())
ax11_3.xaxis.set_major_formatter(plt.NullFormatter())
#ax11_3.set_yticks(np.log10([0.1,1,10,100]))
#ax11_3.set_yticklabels([25,50,75],rotation=0)
ax11_3.hist(np.log10(3.+kurtv),bins=15,range=np.log10([1.,200.]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)




ax2_1.set_xlabel('S/N',fontsize=8,rotation=0)
#ax2_1.set_ylabel(r'variance$^{1/2}$',fontsize=8,rotation=90,labelpad=5)
ax2_1.set_xlim([0.1,200])
ax2_1.set_ylim([10,2000])
ax2_1.set_xscale(u'log')
ax2_1.set_yscale(u'log')
ax2_1.set_xticks([0.1,1,10,100])
ax2_1.set_xticklabels([0.1,1,10,100],rotation=90)
ax2_1.set_yticks([10,100,1000])
ax2_1.set_yticklabels([10,100,1000],rotation=0,fontsize=4)
ax2_1.scatter(snratio,sigteff,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax2_1.scatter(snratio[keep],sigteff[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax2_1.text(0.15,17,'[K]',fontsize=4.,color='k',alpha=1)

#ax2_2.set_xlabel('S/N',fontsize=8,rotation=0)
#ax2_2.set_ylabel(r'skewness',fontsize=8,rotation=90,labelpad=5)
ax2_2.set_xlim([0.1,200])
ax2_2.set_ylim([-9,9])
ax2_2.set_xscale(u'log')
ax2_2.set_yscale(u'linear')
ax2_2.set_xticks([1,10,100])
#ax2_2.set_xticklabels([1,10,100],rotation=0)
ax2_2.xaxis.set_major_formatter(plt.NullFormatter())
ax2_2.set_yticks([-5,0,5])
ax2_2.set_yticklabels([-5,0,5],rotation=0)
ax2_2.scatter(snratio,skewteff,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax2_2.scatter(snratio[keep],skewteff[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)

#ax2_3.set[keep]_xlabel[keep]('S/N',fontsize=8,rotation=0)
#ax2_3.set_ylabel(r'kurtosis',fontsize=8,rotation=90,labelpad=5)
ax2_3.set_xlim([0.1,200])
ax2_3.set_ylim([1,200])
ax2_3.set_xscale(u'log')
ax2_3.set_yscale(u'log')
ax2_3.set_xticks([1,10,100])
#ax2_3.set_xticklabels([1,10,100],rotation=0)
ax2_3.xaxis.set_major_formatter(plt.NullFormatter())
ax2_3.set_yticks([1,10,100])
ax2_3.set_yticklabels([1,10,100],rotation=0)
ax2_3.scatter(snratio,3.+kurtteff,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax2_3.scatter(snratio[keep],3.+kurtteff[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax2_3.text(100,40,r'$T_{\rm eff}$',fontsize=7,color='k',alpha=1,horizontalalignment='right')

#ax22_1.set_xlabel('N',fontsize=8,rotation=0)
#ax22_1.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax22_1.set_xlim([0,100])
ax22_1.set_ylim(np.log10([10,2000]))
ax22_1.set_xscale(u'linear')
ax22_1.set_yscale(u'linear')
#ax22_1.set_xticks([30,60,90])
#ax22_1.set_xticklabels([30,60,90],rotation=90)
ax22_1.yaxis.set_major_formatter(plt.NullFormatter())
ax22_1.xaxis.set_major_formatter(plt.NullFormatter())
ax22_1.set_yticks(np.log10([10,100,1000]))
#ax22_1.set_yticklabels([25,50,75],rotation=0)
ax22_1.hist(np.log10(sigteff),bins=25,range=np.log10([10,2000]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax22_2.set_xlabel('S/N',fontsize=8,rotation=0)
#ax22_2.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax22_2.set_xlim([0,150])
ax22_2.set_ylim([-9,9])
ax22_2.set_xscale(u'linear')
ax22_2.set_yscale(u'linear')
#ax22_2.set_xticks([1,10,100])
#ax22_2.set_xticklabels([1,10,100],rotation=0)
ax22_2.yaxis.set_major_formatter(plt.NullFormatter())
ax22_2.xaxis.set_major_formatter(plt.NullFormatter())
ax22_2.set_yticks([-5,0,5])
#ax22_2.set_yticklabels([25,50,75],rotation=0)
ax22_2.hist(skewteff,bins=25,range=[-9,9],normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax22_3.set_xlabel('S/N',fontsize=8,rotation=0)
#ax22_3.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax22_3.set_xlim([0,150])
ax22_3.set_ylim(np.log10([1.,200.]))
ax22_3.set_xscale(u'linear')
ax22_3.set_yscale(u'linear')
#ax22_3.set_xticks([1,10,100])
#ax22_3.set_xticklabels([1,10,100],rotation=0)
ax22_3.yaxis.set_major_formatter(plt.NullFormatter())
ax22_3.xaxis.set_major_formatter(plt.NullFormatter())
#ax22_3.set_yticks(np.log10([0.1,1,10,100]))
#ax22_3.set_yticklabels([25,50,75],rotation=0)
ax22_3.hist(np.log10(3.+kurtteff),bins=15,range=np.log10([1.,200.]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)




ax3_1.set_xlabel('S/N',fontsize=8,rotation=0)
#ax3_1.set_ylabel(r'variance$^{1/2}$',fontsize=8,rotation=90,labelpad=5)
ax3_1.set_xlim([0.1,200])
ax3_1.set_ylim([0.05,2])
ax3_1.set_xscale(u'log')
ax3_1.set_yscale(u'log')
ax3_1.set_xticks([0.1,1,10,100])
ax3_1.set_xticklabels([0.1,1,10,100],rotation=90)
ax3_1.set_yticks([0.1,1])
ax3_1.set_yticklabels([0.1,1],rotation=0)
ax3_1.scatter(snratio,siglogg,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax3_1.scatter(snratio,siglogg,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax3_1.scatter(snratio[keep],siglogg[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax3_1.text(0.15,0.07,'[dex]',fontsize=4.,color='k',alpha=1)

#ax3_2.set_xlabel('S/N',fontsize=8,rotation=0)
#ax3_2.set_ylabel(r'skewness',fontsize=8,rotation=90,labelpad=5)
ax3_2.set_xlim([0.1,200])
ax3_2.set_ylim([-9,9])
ax3_2.set_xscale(u'log')
ax3_2.set_yscale(u'linear')
ax3_2.set_xticks([1,10,100])
#ax3_2.set_xticklabels([1,10,100],rotation=0)
ax3_2.xaxis.set_major_formatter(plt.NullFormatter())
ax3_2.set_yticks([-5,0,5])
ax3_2.set_yticklabels([-5,0,5],rotation=0)
ax3_2.scatter(snratio,skewlogg,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax3_2.scatter(snratio[keep],skewlogg[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)

#ax3_3.set[keep]_xlabel[keep]('S/N',fontsize=8,rotation=0)
#ax3_3.set_ylabel(r'kurtosis',fontsize=8,rotation=90,labelpad=5)
ax3_3.set_xlim([0.1,200])
ax3_3.set_ylim([1,200])
ax3_3.set_xscale(u'log')
ax3_3.set_yscale(u'log')
ax3_3.set_xticks([1,10,100])
#ax3_3.set_xticklabels([1,10,100],rotation=0)
ax3_3.xaxis.set_major_formatter(plt.NullFormatter())
ax3_3.set_yticks([1,10,100])
ax3_3.set_yticklabels([1,10,100],rotation=0)
ax3_3.scatter(snratio,3.+kurtlogg,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax3_3.scatter(snratio[keep],3.+kurtlogg[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax3_3.text(100,40,r'$\log g$',fontsize=7,color='k',alpha=1,horizontalalignment='right')

#ax33_1.set_xlabel('N',fontsize=8,rotation=0)
#ax33_1.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax33_1.set_xlim([0,100])
ax33_1.set_ylim(np.log10([0.05,2]))
ax33_1.set_xscale(u'linear')
ax33_1.set_yscale(u'linear')
#ax33_1.set_xticks([30,60,90])
#ax33_1.set_xticklabels([30,60,90],rotation=90)
ax33_1.yaxis.set_major_formatter(plt.NullFormatter())
ax33_1.xaxis.set_major_formatter(plt.NullFormatter())
ax33_1.set_yticks(np.log10([0.1,2]))
#ax33_1.set_yticklabels([25,50,75],rotation=0)
ax33_1.hist(np.log10(siglogg),bins=25,range=np.log10([0.1,2]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax33_2.set_xlabel('S/N',fontsize=8,rotation=0)
#ax33_2.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax33_2.set_xlim([0,150])
ax33_2.set_ylim([-9,9])
ax33_2.set_xscale(u'linear')
ax33_2.set_yscale(u'linear')
#ax33_2.set_xticks([1,10,100])
#ax33_2.set_xticklabels([1,10,100],rotation=0)
ax33_2.yaxis.set_major_formatter(plt.NullFormatter())
ax33_2.xaxis.set_major_formatter(plt.NullFormatter())
ax33_2.set_yticks([-5,0,5])
#ax33_2.set_yticklabels([25,50,75],rotation=0)
ax33_2.hist(skewlogg,bins=25,range=[-9,9],normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax33_3.set_xlabel('S/N',fontsize=8,rotation=0)
#ax33_3.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax33_3.set_xlim([0,150])
ax33_3.set_ylim(np.log10([1.,200.]))
ax33_3.set_xscale(u'linear')
ax33_3.set_yscale(u'linear')
#ax33_3.set_xticks([1,10,100])
#ax33_3.set_xticklabels([1,10,100],rotation=0)
ax33_3.yaxis.set_major_formatter(plt.NullFormatter())
ax33_3.xaxis.set_major_formatter(plt.NullFormatter())
#ax33_3.set_yticks(np.log10([0.1,1,10,100]))
#ax33_3.set_yticklabels([25,50,75],rotation=0)
ax33_3.hist(np.log10(3.+kurtlogg),bins=15,range=np.log10([1.,200.]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)




ax4_1.set_xlabel('S/N',fontsize=8,rotation=0)
#ax4_1.set_ylabel(r'variance$^{1/2}$',fontsize=8,rotation=90,labelpad=5)
ax4_1.set_xlim([0.1,200])
ax4_1.set_ylim([0.05,2])
ax4_1.set_xscale(u'log')
ax4_1.set_yscale(u'log')
ax4_1.set_xticks([0.1,1,10,100])
ax4_1.set_xticklabels([0.1,1,10,100],rotation=90)
ax4_1.set_yticks([0.1,1])
ax4_1.set_yticklabels([0.1,1],rotation=0)
ax4_1.scatter(snratio,sigfeh,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax4_1.scatter(snratio[keep],sigfeh[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax4_1.text(0.15,0.07,'[dex]',fontsize=4.,color='k',alpha=1)

#ax4_2.set_xlabel('S/N',fontsize=8,rotation=0)
#ax4_2.set_ylabel(r'skewness',fontsize=8,rotation=90,labelpad=5)
ax4_2.set_xlim([0.1,200])
ax4_2.set_ylim([-9,9])
ax4_2.set_xscale(u'log')
ax4_2.set_yscale(u'linear')
ax4_2.set_xticks([1,10,100])
#ax4_2.set_xticklabels([1,10,100],rotation=0)
ax4_2.xaxis.set_major_formatter(plt.NullFormatter())
ax4_2.set_yticks([-5,0,5])
ax4_2.set_yticklabels([-5,0,5],rotation=0)
ax4_2.scatter(snratio,skewfeh,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax4_2.scatter(snratio[keep],skewfeh[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)

#ax4_3.set_xlabel[keep]('S[keep]/N',fontsize=8,rotation=0)
#ax4_3.set_ylabel(r'kurtosis',fontsize=8,rotation=90,labelpad=5)
ax4_3.set_xlim([0.1,200])
ax4_3.set_ylim([1,200])
ax4_3.set_xscale(u'log')
ax4_3.set_yscale(u'log')
ax4_3.set_xticks([1,10,100])
#ax4_3.set_xticklabels([1,10,100],rotation=0)
ax4_3.xaxis.set_major_formatter(plt.NullFormatter())
ax4_3.set_yticks([1,10,100])
ax4_3.set_yticklabels([1,10,100],rotation=0)
ax4_3.scatter(snratio,3.+kurtfeh,s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='k',rasterized=True)
ax4_3.scatter(snratio[keep],3.+kurtfeh[keep],s=1,lw=0,edgecolor='none',alpha=0.75,marker='o',color='r',rasterized=True)
ax4_3.text(100,40,r'[Fe/H]',fontsize=7,color='k',alpha=1,horizontalalignment='right')

#ax44_1.set_xlabel('N',fontsize=8,rotation=0)
#ax44_1.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax44_1.set_xlim([0,100])
ax44_1.set_ylim(np.log10([0.05,2]))
ax44_1.set_xscale(u'linear')
ax44_1.set_yscale(u'linear')
#ax44_1.set_xticks([30,60,90])
#ax44_1.set_xticklabels([30,60,90],rotation=90)
ax44_1.yaxis.set_major_formatter(plt.NullFormatter())
ax44_1.xaxis.set_major_formatter(plt.NullFormatter())
ax44_1.set_yticks(np.log10([0.1,2]))
#ax44_1.set_yticklabels([25,50,75],rotation=0)
ax44_1.hist(np.log10(sigfeh),bins=25,range=np.log10([0.1,2]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax44_2.set_xlabel('S/N',fontsize=8,rotation=0)
#ax44_2.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax44_2.set_xlim([0,150])
ax44_2.set_ylim([-9,9])
ax44_2.set_xscale(u'linear')
ax44_2.set_yscale(u'linear')
#ax44_2.set_xticks([1,10,100])
#ax44_2.set_xticklabels([1,10,100],rotation=0)
ax44_2.yaxis.set_major_formatter(plt.NullFormatter())
ax44_2.xaxis.set_major_formatter(plt.NullFormatter())
ax44_2.set_yticks([-5,0,5])
#ax44_2.set_yticklabels([25,50,75],rotation=0)
ax44_2.hist(skewfeh,bins=25,range=[-9,9],normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)

#ax44_3.set_xlabel('S/N',fontsize=8,rotation=0)
#ax44_3.set_ylabel(r'N',fontsize=8,rotation=90,labelpad=5)
#ax44_3.set_xlim([0,150])
ax44_3.set_ylim(np.log10([1.,200.]))
ax44_3.set_xscale(u'linear')
ax44_3.set_yscale(u'linear')
#ax44_3.set_xticks([1,10,100])
#ax44_3.set_xticklabels([1,10,100],rotation=0)
ax44_3.yaxis.set_major_formatter(plt.NullFormatter())
ax44_3.xaxis.set_major_formatter(plt.NullFormatter())
#ax44_3.set_yticks(np.log10([0.1,1,10,100]))
#ax44_3.set_yticklabels([25,50,75],rotation=0)
ax44_3.hist(np.log10(3.+kurtfeh),bins=15,range=np.log10([1.,200.]),normed=True,align='mid',histtype='step',rwidth=0,orientation='horizontal',color='k',linewidth=0.25,alpha=1)


plotfilename='m2fsmedres_jul15_diagnostics.pdf'
plt.savefig(plotfilename,dpi=400)
plt.show()
plt.close()
