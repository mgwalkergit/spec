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
plt.rc('legend',frameon='True')
plt.rc('legend',framealpha=0.5)
plt.rc('legend',borderpad=1)
plt.rc('legend',borderaxespad=1)
plt.rc('legend',scatterpoints=1)
plt.rc('legend',numpoints=1)
plt.rc('xtick.major',size=2)
plt.rc('ytick.major',size=2)
plt.rc('xtick.minor',size=2)
plt.rc('ytick.minor',size=2)
plt.rc('axes',axisbelow='True')
plt.rc('axes',grid='False')
plt.rc('axes',facecolor='white')
plt.rc('axes',linewidth=0.5)
plt.rc('xtick.major',width=0.5)
plt.rc('ytick.major',width=0.5)
plt.rc('xtick.minor',width=0.5)
plt.rc('ytick.minor',width=0.5)
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
            
data1='mcconnachie.dat'
data2='harris2.dat'
data3='harris3.dat'
data4='gru1teffprior_gradientpost_equal_weights.dat'
data5='ret2_gradientpost_equal_weights.dat'
data6='tuc2teffprior_gradientpost_equal_weights.dat'
out1='tuc2gru1teffprior_scaling_newcommands.tex'
g1=open(out1,'w')

with open(data1) as f: # read data file
    data=f.readlines()

parent=[]
dsph=[]
rah=[]
ram=[]
ras=[]
chardecd=[]
decm=[]
decs=[]
glon=[]
glat=[]
distance=[]
sigdistance=[]
vsys=[]
sigvsys=[]
vmag=[]
sigvmag=[]
rhalfarcmin=[]
sigrhalfarcmin=[]
muv=[]
sigmuv=[]
ellip=[]
absvmag=[]
sigabsvmag=[]
rhalf=[]
sigrhalf=[]
vdisp=[]
sigvdisp=[]
feh=[]
sigfeh=[]
method=[]

for line in data: # fill arrays
    p=line.split()
    parent.append(str(p[0]))
    dsph.append(str(p[1]))
    rah.append(float(p[3]))
    ram.append(float(p[4]))
    ras.append(float(p[5]))
    chardecd.append(float(p[6]))
    decm.append(float(p[7]))
    decs.append(str(p[8]))
    glon.append(float(p[9]))
    glat.append(float(p[10]))
    distance.append(float(p[11]))
    sigdistance.append(float(p[12]))
    vsys.append(float(p[14]))
    sigvsys.append(float(p[15]))
    vmag.append(float(p[16]))
    sigvmag.append(float(p[17]))
    rhalfarcmin.append(float(p[19]))
    sigrhalfarcmin.append(float(p[20]))
    muv.append(float(p[22]))
    sigmuv.append(float(p[23]))
    ellip.append(float(p[25]))
    absvmag.append(float(p[26]))
    sigabsvmag.append(float(p[27]))
    rhalf.append(float(p[29]))
    sigrhalf.append(float(p[30]))
    vdisp.append(float(p[32]))
    sigvdisp.append(float(p[33]))
    feh.append(float(p[35]))
    sigfeh.append(float(p[36]))
    method.append(str(p[37]))

parent=np.array(parent)
dsph=np.array(dsph)
rah=np.array(rah)
ram=np.array(ram)
ras=np.array(ras)
chardecd=np.array(chardecd)
decm=np.array(decm)
decs=np.array(decs)
glon=np.array(glon)
glat=np.array(glat)
distance=np.array(distance)
sigdistance=np.array(sigdistance)
vsys=np.array(vsys)
sigvsys=np.array(sigvsys)
vmag=np.array(vmag)
sigvmag=np.array(sigvmag)
rhalfarcmin=np.array(rhalfarcmin)
sigrhalfarcmin=np.array(sigrhalfarcmin)
muv=np.array(muv)
sigmuv=np.array(sigmuv)
ellip=np.array(ellip)
absvmag=np.array(absvmag)
sigabsvmag=np.array(sigabsvmag)
rhalf=np.array(rhalf)
sigrhalf=np.array(sigrhalf)
vdisp=np.array(vdisp)
sigvdisp=np.array(sigvdisp)
feh=np.array(feh)
sigfeh=np.array(sigfeh)
method=np.array(method)

for i in range(0,len(rhalf)):
    if ellip[i] > 0.:
        rhalf[i]=rhalf[i]*np.sqrt(1.-ellip[i])########convert to geometric mean radius
        sigrhalf[i]=sigrhalf[i]*np.sqrt(1.-ellip[i])

with open(data2) as f: # read data file
    data=f.readlines()

gcname=[]
gcfeh=[]
crap=[]
gcred=[]
crap=[]
gcdistmodulus=[]
gcvmag=[]
gcabsvmag=[]

for line in data: # fill arrays
    p=line.split()

    gcname.append(str(p[0]))
    gcfeh.append(float(p[1]))
    gcred.append(float(p[3]))
    crap.append(float(p[4]))
    gcdistmodulus.append(float(p[5]))
    gcvmag.append(float(p[6]))
    gcabsvmag.append(float(p[7]))

gcname=np.array(gcname)
gcfeh=np.array(gcfeh)
gcred=np.array(gcred)
crap=np.array(crap)
gcdistmodulus=np.array(gcdistmodulus)
gcvmag=np.array(gcvmag)
gcabsvmag=np.array(gcabsvmag)

gcdist=10.**((gcdistmodulus+5.)/5.)

with open(data3) as f: # read data file
    data=f.readlines()

gcname2=[]
gcvr=[]
gcsigvr=[]
gcvlsr=[]
gcvdisp=[]
gcsigvdisp=[]
gcc=[]
gcrcore=[]
gcrhalf=[]

for line in data: # fill arrays
    p=line.split()
    gcname2.append(str(p[0]))
    gcvr.append(float(p[1]))
    gcsigvr.append(float(p[2]))
    gcvlsr.append(float(p[3]))
    gcvdisp.append(float(p[4]))
    gcsigvdisp.append(float(p[5]))
    gcc.append(float(p[6]))
    gcrcore.append(float(p[7]))
    gcrhalf.append(float(p[8]))

gcvr=np.array(gcvr)
gcsigvr=np.array(gcsigvr)
gcvlsr=np.array(gcvlsr)
gcvdisp=np.array(gcvdisp)
gcsigvdisp=np.array(gcsigvdisp)
gcc=np.array(gcc)
gcrcore=np.array(gcrcore)
gcrhalf=np.array(gcrhalf)

gcrhalf=gcdist*np.tan(gcrhalf/60.*np.pi/180.)

with open(data4) as f: # read data file
    data=f.readlines()

gru1vmean=[]
gru1vdisp=[]
gru1vgrad=[]
gru1vtheta=[]
gru1fehmean=[]
gru1fehdisp=[]
gru1fehgrad=[]
gru1like=[]

for line in data: # fill arrays
    p=line.split()
    gru1vmean.append(float(p[0]))
    gru1vdisp.append(float(p[1]))
    gru1vgrad.append(float(p[2]))
    gru1vtheta.append(float(p[3]))
    gru1fehmean.append(float(p[4]))
    gru1fehdisp.append(float(p[5]))
    gru1fehgrad.append(float(p[6]))
    gru1like.append(float(p[7]))

gru1vmean=np.array(gru1vmean)
gru1vdisp=np.array(gru1vdisp)
gru1vgrad=np.array(gru1vgrad)
gru1vtheta=np.array(gru1vtheta)
gru1fehmean=np.array(gru1fehmean)
gru1fehdisp=np.array(gru1fehdisp)
gru1fehgrad=np.array(gru1fehgrad)
gru1like=np.array(gru1like)

gru1_vdisp0=np.array([np.median(gru1vdisp)])
gru1_sigvdisp0=np.array([np.std(gru1vdisp)])
gru1_feh=np.array([np.median(gru1fehmean)])
gru1_sigfeh=np.array([np.std(gru1fehmean)])
gru1_rhalf0=62.
gru1_sigrhalf0=30.
gru1_rhalf=np.random.normal(loc=gru1_rhalf0,scale=gru1_sigrhalf0,size=len(gru1vdisp))
#gru1_sigrhalf=np.array([1.])
gru1_absvmag0=-3.4
gru1_sigabsvmag0=0.3
gru1_absvmag=np.random.normal(loc=gru1_absvmag0,scale=gru1_sigabsvmag0,size=len(gru1vdisp))
#gru1_sigabsvmag=np.array([0.1])

with open(data5) as f: # read data file
    data=f.readlines()

ret2vmean=[]
ret2vdisp=[]
ret2vgrad=[]
ret2vtheta=[]
ret2fehmean=[]
ret2fehdisp=[]
ret2fehgrad=[]
ret2like=[]

for line in data: # fill arrays
    p=line.split()
    ret2vmean.append(float(p[0]))
    ret2vdisp.append(float(p[1]))
    ret2vgrad.append(float(p[2]))
    ret2vtheta.append(float(p[3]))
    ret2fehmean.append(float(p[4]))
    ret2fehdisp.append(float(p[5]))
    ret2fehgrad.append(float(p[6]))
    ret2like.append(float(p[7]))

ret2vmean=np.array(ret2vmean)
ret2vdisp=np.array(ret2vdisp)
ret2vgrad=np.array(ret2vgrad)
ret2vtheta=np.array(ret2vtheta)
ret2fehmean=np.array(ret2fehmean)
ret2fehdisp=np.array(ret2fehdisp)
ret2fehgrad=np.array(ret2fehgrad)
ret2like=np.array(ret2like)

ret2_vdisp0=np.array([np.median(ret2vdisp)])
ret2_sigvdisp0=np.array([np.std(ret2vdisp)])
ret2_feh=np.array([np.median(ret2fehmean)])
ret2_sigfeh=np.array([np.std(ret2fehmean)])
ret2_rhalf0=32.
ret2_sigrhalf0=2.
ret2_rhalf=np.random.normal(loc=ret2_rhalf0,scale=ret2_sigrhalf0,size=len(ret2vdisp))
#ret2_sigrhalf=np.array([1.])
ret2_absvmag0=-2.7
ret2_sigabsvmag0=0.1
ret2_absvmag=np.random.normal(loc=ret2_absvmag0,scale=0.2,size=len(ret2vdisp))
#ret2_sigabsvmag=np.array([0.1])


with open(data6) as f: # read data file
    data=f.readlines()

tuc2vmean=[]
tuc2vdisp=[]
tuc2vgrad=[]
tuc2vtheta=[]
tuc2fehmean=[]
tuc2fehdisp=[]
tuc2fehgrad=[]
tuc2like=[]

for line in data: # fill arrays
    p=line.split()
    tuc2vmean.append(float(p[0]))
    tuc2vdisp.append(float(p[1]))
    tuc2vgrad.append(float(p[2]))
    tuc2vtheta.append(float(p[3]))
    tuc2fehmean.append(float(p[4]))
    tuc2fehdisp.append(float(p[5]))
    tuc2fehgrad.append(float(p[6]))
    tuc2like.append(float(p[7]))

tuc2vmean=np.array(tuc2vmean)
tuc2vdisp=np.array(tuc2vdisp)
tuc2vgrad=np.array(tuc2vgrad)
tuc2vtheta=np.array(tuc2vtheta)
tuc2fehmean=np.array(tuc2fehmean)
tuc2fehdisp=np.array(tuc2fehdisp)
tuc2fehgrad=np.array(tuc2fehgrad)
tuc2like=np.array(tuc2like)

tuc2_vdisp0=np.array([np.median(tuc2vdisp)])
tuc2_sigvdisp0=np.array([np.std(tuc2vdisp)])
tuc2_feh=np.array([np.median(tuc2fehmean)])
tuc2_sigfeh=np.array([np.std(tuc2fehmean)])
tuc2_rhalf0=165.
tuc2_sigrhalf0=28.
tuc2_rhalf=np.random.normal(loc=tuc2_rhalf0,scale=tuc2_sigrhalf0,size=len(tuc2vdisp))
#tuc2_sigrhalf=np.array([1.])
tuc2_absvmag0=-3.8
tuc2_sigabsvmag0=0.1
tuc2_absvmag=np.random.normal(loc=tuc2_absvmag0,scale=tuc2_sigabsvmag0,size=len(tuc2vdisp))
#tuc2_sigabsvmag=np.array([0.1])


luminosity=10.**((absvmag-4.83)/(-2.5))
sigluminosity=np.log(10.)/2.5*10.**((absvmag-4.83)/(-2.5))*sigabsvmag
mrhalf=1./0.0043*rhalf*vdisp**2
sigmrhalf=np.sqrt(2.*rhalf/0.0043*vdisp*(sigvdisp**2)+((1./0.0043*(vdisp**2))**2)*sigrhalf**2)
mlratio=mrhalf/(luminosity)
sigmlratio=np.sqrt((sigmrhalf**2)/(luminosity**2)+((mrhalf/luminosity**2)**2)*sigluminosity**2)

gru1_luminosity0=10.**((gru1_absvmag0-4.83)/(-2.5))
gru1_luminosity=10.**((gru1_absvmag-4.83)/(-2.5))
gru1_sigluminosity0=np.log(10.)/2.5*10.**((gru1_absvmag0-4.83)/(-2.5))*gru1_sigabsvmag0
gru1_mrhalf0=1./0.0043*gru1_rhalf0*gru1_vdisp0**2
gru1_mrhalf=1./0.0043*gru1_rhalf*gru1vdisp**2
gru1_sigmrhalf0=np.sqrt(2.*gru1_rhalf0/0.0043*gru1_vdisp0*(gru1_sigvdisp0**2)+((1./0.0043*(gru1_vdisp0**2))**2)*gru1_sigrhalf0**2)
gru1_mlratio=gru1_mrhalf/(gru1_luminosity)
gru1_mlratio0=gru1_mrhalf0/gru1_luminosity0
gru1_sigmlratio0=np.sqrt((gru1_sigmrhalf0**2)/(gru1_luminosity0**2)+((gru1_mrhalf0/gru1_luminosity0**2)**2)*gru1_sigluminosity0**2)

ret2_luminosity0=10.**((ret2_absvmag0-4.83)/(-2.5))
ret2_luminosity=10.**((ret2_absvmag-4.83)/(-2.5))
ret2_sigluminosity0=np.log(10.)/2.5*10.**((ret2_absvmag0-4.83)/(-2.5))*ret2_sigabsvmag0
ret2_mrhalf0=1./0.0043*ret2_rhalf0*ret2_vdisp0**2
ret2_mrhalf=1./0.0043*ret2_rhalf*ret2vdisp**2
ret2_sigmrhalf0=np.sqrt(2.*ret2_rhalf0/0.0043*ret2_vdisp0*(ret2_sigvdisp0**2)+((1./0.0043*(ret2_vdisp0**2))**2)*ret2_sigrhalf0**2)
ret2_mlratio=ret2_mrhalf/(ret2_luminosity)
ret2_mlratio0=ret2_mrhalf0/ret2_luminosity0
ret2_sigmlratio0=np.sqrt((ret2_sigmrhalf0**2)/(ret2_luminosity0**2)+((ret2_mrhalf0/ret2_luminosity0**2)**2)*ret2_sigluminosity0**2)

tuc2_luminosity0=10.**((tuc2_absvmag0-4.83)/(-2.5))
tuc2_luminosity=10.**((tuc2_absvmag-4.83)/(-2.5))
tuc2_sigluminosity0=np.log(10.)/2.5*10.**((tuc2_absvmag0-4.83)/(-2.5))*tuc2_sigabsvmag0
tuc2_mrhalf0=1./0.0043*tuc2_rhalf0*tuc2_vdisp0**2
tuc2_mrhalf=1./0.0043*tuc2_rhalf*tuc2vdisp**2
tuc2_sigmrhalf0=np.sqrt(2.*tuc2_rhalf0/0.0043*tuc2_vdisp0*(tuc2_sigvdisp0**2)+((1./0.0043*(tuc2_vdisp0**2))**2)*tuc2_sigrhalf0**2)
tuc2_mlratio=tuc2_mrhalf/(tuc2_luminosity)
tuc2_mlratio0=tuc2_mrhalf0/tuc2_luminosity0
tuc2_sigmlratio0=np.sqrt((tuc2_sigmrhalf0**2)/(tuc2_luminosity0**2)+((tuc2_mrhalf0/tuc2_luminosity0**2)**2)*tuc2_sigluminosity0**2)

gcluminosity=10.**((gcabsvmag-4.83)/(-2.5))
#gcsigluminosity=np.log(10.)/2.5*10.**((gcabsvmag-4.83)/(-2.5))*gcsigabsvmag
gcmrhalf=1./0.0043*gcrhalf*gcvdisp**2
#gcsigmrhalf=np.sqrt(2.*gcrhalf*gcvdisp*(gcsigvdisp**2)+((1.*(gcvdisp**2))**2)*gcsigrhalf**2)
gcmlratio=gcmrhalf/(gcluminosity/2.)
#gcsigmlratio=np.sqrt((gcsigmrhalf**2)/(gcluminosity**2)+((gcmrhalf/gcluminosity**2)**2)*gcsigluminosity**2)

gs=plt.GridSpec(10,10) # define multi-panel plot
gs.update(wspace=0,hspace=0) # specify inter-panel spacing
fig=plt.figure(figsize=(6,6)) # define plot size

ax0_0=fig.add_subplot(gs[0:3,0:4])
ax1_0=fig.add_subplot(gs[3:6,0:4])
ax2_0=fig.add_subplot(gs[6:9,0:4])


#ax0_0.set_xlabel(r'$M_{\rm V}$ [mag]',fontsize=8,rotation=0)
ax0_0.set_ylabel(r'$R_{\rm h}$ [pc]',fontsize=13,rotation=90,labelpad=5)
ax0_0.set_xlim([0,-15])
ax0_0.set_ylim([1,3000])
ax0_0.set_xscale(u'linear')
ax0_0.set_yscale(u'log')
ax0_0.set_xticks([0,-2.5,-5,-7.5,-10,-12.5,-15])
ax0_0.xaxis.set_major_formatter(plt.NullFormatter())
#ax0_0.set_xticklabels(vticks,rotation=0)
ax0_0.set_yticks([10,100,1000])
ax0_0.set_yticklabels([10,100,1000],rotation=0,fontsize=10)
ax0_0.scatter(gcabsvmag[gcvdisp > 0.],gcrhalf[gcvdisp > 0.],s=2,lw=0,edgecolor='none',alpha=0.99,marker='o',color='0.5',rasterized=True)
#ax0_0.scatter(np.median(gru1_absvmag),gru1_rhalf,s=20,lw=0,edgecolor='none',alpha=0.99,marker='*',color='k',rasterized=True)
#ax0_0.scatter(np.median(ret2_absvmag),ret2_rhalf,s=20,lw=0,edgecolor='none',alpha=0.99,marker='o',color='k',rasterized=True)
ax0_0.scatter(np.median(tuc2_absvmag),np.median(tuc2_rhalf),s=10,lw=0,edgecolor='none',alpha=0.99,marker='s',color='k',rasterized=True)
ax0_0.scatter(np.median(gru1_absvmag),np.median(gru1_rhalf),s=20,lw=0,edgecolor='none',alpha=0.99,marker='^',color='k',rasterized=True)
ax0_0.scatter(np.median(ret2_absvmag),np.median(ret2_rhalf),s=10,lw=0,edgecolor='none',alpha=0.99,marker='D',color='k',rasterized=True)
ax0_0.errorbar(absvmag[(parent == 'MW') & (sigfeh > 0)],rhalf[(parent == 'MW') & (sigfeh > 0)],xerr=sigabsvmag[(parent == 'MW') & (sigfeh > 0)],yerr=sigrhalf[(parent == 'MW') & (sigfeh > 0)],elinewidth=0.25,fmt='o',capsize=0,alpha=1,color='b',rasterized=True)
ax0_0.errorbar(absvmag[(parent == 'M31') & (sigfeh >0)],rhalf[(parent == 'M31') & (sigfeh >0)],xerr=sigabsvmag[(parent == 'M31') & (sigfeh > 0)],yerr=sigrhalf[(parent == 'M31') & (sigfeh > 0)],elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='r',rasterized=True)
#ax0_0.errorbar(np.median(gru1_absvmag),gru1_rhalf,xerr=np.std(gru1_absvmag),yerr=gru1_sigrhalf,elinewidth=1,fmt='o',capsize=0,alpha=1,color='k',rasterized=True)
ax0_0.errorbar(np.median(ret2_absvmag),np.median(ret2_rhalf),xerr=np.std(ret2_absvmag),yerr=ret2_sigrhalf0,elinewidth=0.5,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
ax0_0.errorbar(np.median(tuc2_absvmag),np.median(tuc2_rhalf),xerr=np.std(tuc2_absvmag),yerr=tuc2_sigrhalf0,elinewidth=0.5,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
ax0_0.errorbar(np.median(gru1_absvmag),np.median(gru1_rhalf),xerr=np.std(gru1_absvmag),yerr=gru1_sigrhalf0,elinewidth=0.5,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
#ax0_0.text(np.median(gru1_absvmag)-0.15,gru1_rhalf[0]-0.4,'Gru 1',fontsize=7,color='k')
#ax0_0.text(np.median(ret2_absvmag)-0.2,ret2_rhalf[0]+0.2,'Ret 2',fontsize=7,color='k')
ax0_0.text(np.median(tuc2_absvmag)+2.7,np.median(tuc2_rhalf)-0.2,'Tuc 2',fontsize=7,color='k')
ax0_0.text(np.median(gru1_absvmag)+2.5,np.median(gru1_rhalf)-0.2,'Gru 1',fontsize=7,color='k')




#ax1_0.set_xlabel(r'$M_{\rm V}$ [mag]',fontsize=8,rotation=0)
ax1_0.set_ylabel(r'[Fe/H]',fontsize=13,rotation=90,labelpad=5)
ax1_0.set_xlim([0,-15])
ax1_0.set_ylim([-3.5,0.5])
ax1_0.set_xscale(u'linear')
ax1_0.set_yscale(u'linear')
ax1_0.set_xticks([0,-2.5,-5,-7.5,-10,-12.5,-15])
ax1_0.xaxis.set_major_formatter(plt.NullFormatter())
#ax1_0.set_xticklabels(vticks,rotation=0)
ax1_0.set_yticks([-3,-2,-1,0])
ax1_0.set_yticklabels([-3,-2,-1,0],rotation=0,fontsize=10)
ax1_0.scatter(gcabsvmag[gcvdisp > 0.],gcfeh[gcvdisp > 0.],s=2,lw=0,edgecolor='none',alpha=0.99,marker='o',color='0.5',rasterized=True)
#ax1_0.scatter(np.median(gru1_absvmag),gru1_feh,s=20,lw=0,edgecolor='none',alpha=0.99,marker='*',color='k',rasterized=True)
ax1_0.scatter(np.median(tuc2_absvmag),tuc2_feh,s=10,lw=0,edgecolor='none',alpha=0.99,marker='s',color='k',rasterized=True)
ax1_0.scatter(np.median(gru1_absvmag),gru1_feh,s=20,lw=0,edgecolor='none',alpha=0.99,marker='^',color='k',rasterized=True)
ax1_0.scatter(np.median(ret2_absvmag),ret2_feh,s=10,lw=0,edgecolor='none',alpha=0.99,marker='D',color='k',rasterized=True)
ax1_0.errorbar(absvmag[(parent == 'MW') & (sigfeh > 0)],feh[(parent == 'MW') & (sigfeh > 0)],xerr=sigabsvmag[(parent == 'MW') & (sigfeh > 0)],yerr=sigfeh[(parent == 'MW') & (sigfeh > 0)],elinewidth=0.25,fmt='o',capsize=0,alpha=1,color='b',rasterized=True)

ax1_0.errorbar(absvmag[(parent == 'M31') & (sigfeh >0)],feh[(parent == 'M31') & (sigfeh >0)],xerr=sigabsvmag[(parent == 'M31') & (sigfeh > 0)],yerr=sigfeh[(parent == 'M31') & (sigfeh > 0)],elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='r',rasterized=True)
#ax1_0.errorbar(np.median(gru1_absvmag),gru1_feh,xerr=np.std(gru1_absvmag),yerr=gru1_sigfeh,elinewidth=1,fmt='o',capsize=0,alpha=1,color='k',rasterized=True)
ax1_0.errorbar(np.median(ret2_absvmag),ret2_feh,xerr=np.std(ret2_absvmag),yerr=ret2_sigfeh,elinewidth=0.5,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
ax1_0.errorbar(np.median(tuc2_absvmag),tuc2_feh,xerr=np.std(tuc2_absvmag),yerr=tuc2_sigfeh,elinewidth=0.5,fmt='o',capsize=0,alpha=1,color='k',rasterized=True)
ax1_0.errorbar(np.median(gru1_absvmag),gru1_feh,xerr=np.std(gru1_absvmag),yerr=gru1_sigfeh,elinewidth=0.5,fmt='o',capsize=0,alpha=1,color='k',rasterized=True)
#ax1_0.text(np.median(gru1_absvmag)-0.15,gru1_feh[0]-0.4,'Gru 1',fontsize=7,color='k')
#ax1_0.text(np.median(ret2_absvmag)-0.2,ret2_feh[0]+0.2,'Ret 2',fontsize=7,color='k')
ax1_0.text(np.median(tuc2_absvmag)-0.2,tuc2_feh[0]-0.2,'Tuc 2',fontsize=7,color='k')
ax1_0.text(np.median(gru1_absvmag)-0.3,gru1_feh[0]-0.2,'Gru 1',fontsize=7,color='k')

ax2_0.set_xlabel(r'$M_{\rm V}$ [mag]',fontsize=13,rotation=0)
ax2_0.set_ylabel(r'$\frac{R_{h}\sigma^2_{v_{\rm los}}}{GL_{\rm V}}$  [M/L$_{\rm V}$]$_{\odot}$',fontsize=13,rotation=90,labelpad=3)
ax2_0.set_xlim([0,-15])
ax2_0.set_ylim([0.1,9999])
ax2_0.set_xscale(u'linear')
ax2_0.set_yscale(u'log')
ax2_0.set_xticks([0,-2.5,-5,-7.5,-10,-12.5,-15])
ax2_0.set_yticks([0.1,1,10,100])
ax2_0.set_xticklabels([0,-2.5,-5,-7.5,-10,-12.5,-15],rotation=0,fontsize=10)
ax2_0.set_yticklabels([0.1,1,10,100],rotation=0,fontsize=10)
#ax2_0.set_yticks([1,10,100,1000])
#ax2_0.xaxis.set_minor_formatter(tick.FuncFormatter(showOnlySomeTicks))

#ax2_0.scatter(np.median(gru1_absvmag),np.median(gru1_mlratio),s=20,lw=0,edgecolor='none',alpha=0.99,marker='s',color='b',rasterized=True,label='Ret 2')
#ax2_0.scatter(np.median(gru1_absvmag0),np.median(gru1_mlratio0),s=20,lw=0,edgecolor='none',alpha=0.99,marker='*',color='k',rasterized=True,label='Gru 1')
ax2_0.scatter(np.median(tuc2_absvmag0),np.median(tuc2_mlratio0),s=10,lw=0,edgecolor='none',alpha=0.99,marker='s',color='k',rasterized=True,label='Tuc 2')
ax2_0.text(np.median(tuc2_absvmag)-0.2,tuc2_mlratio0[0]-0.2,'Tuc 2',fontsize=7,color='k')
ax2_0.scatter(np.median(gru1_absvmag0),np.median(gru1_mlratio0),s=20,lw=0,edgecolor='none',alpha=0.99,marker='^',color='k',rasterized=True,label='Gru 1')
ax2_0.text(np.median(gru1_absvmag)-0.2,gru1_mlratio0[0]-0.2,'Gru 1',fontsize=7,color='k')
ax2_0.scatter(np.median(ret2_absvmag0),np.median(ret2_mlratio0),s=10,lw=0,edgecolor='none',alpha=0.99,marker='D',color='k',rasterized=True,label='Ret 2')

ax2_0.errorbar(absvmag[(parent == 'MW') & (sigvdisp > 0)],mlratio[(parent == 'MW') & (sigvdisp > 0)],xerr=sigabsvmag[(parent == 'MW') & (sigvdisp > 0)],yerr=sigmlratio[(parent == 'MW') & (sigvdisp > 0)],elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='b',rasterized=True,label='MW dSph')
ax2_0.errorbar(absvmag[(parent == 'M31') & (sigvdisp > 0)],mlratio[(parent == 'M31') & (sigvdisp > 0)],xerr=sigabsvmag[(parent == 'M31') & (sigvdisp > 0)],yerr=sigmlratio[(parent == 'M31') & (sigvdisp > 0)],elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='r',rasterized=True,label='M31 dSph')
ax2_0.scatter(gcabsvmag[gcvdisp > 0.],gcmlratio[gcvdisp > 0.],s=2,lw=0,edgecolor='none',alpha=0.99,marker='o',color='0.5',rasterized=True,label='GCs')


lower_error=np.array([(gru1_mlratio[gru1like == max(gru1like)])[0]-mycode2.hpd(gru1_mlratio,0.32)[0]])
upper_error=np.array([mycode2.hpd(gru1_mlratio,0.32)[1]-(gru1_mlratio[gru1like == max(gru1like)])[0]])       
#print (gru1_mlratio[gru1like == max(gru1like)])[0]
#print lower_error
#print upper_error

gru1_lower_error=[np.percentile(gru1_mlratio,50)-np.percentile(gru1_mlratio,16)]
gru1_lower_error=[10000]#because velocity dispersion is unresolved!!!!
gru1_upper_error=[np.percentile(gru1_mlratio,84)-np.percentile(gru1_mlratio,50)]
#ax2_0.errorbar(np.median(gru1_absvmag),np.median(gru1_mlratio),xerr=np.std(gru1_absvmag),yerr=[gru1_lower_error,gru1_upper_error],elinewidth=1,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
ret2_lower_error=[np.percentile(ret2_mlratio,50)-np.percentile(ret2_mlratio,16)]
ret2_upper_error=[np.percentile(ret2_mlratio,84)-np.percentile(ret2_mlratio,50)]
tuc2_lower_error=[np.percentile(tuc2_mlratio,50)-np.percentile(tuc2_mlratio,16)]
tuc2_upper_error=[np.percentile(tuc2_mlratio,84)-np.percentile(tuc2_mlratio,50)]
ax2_0.errorbar(np.median(tuc2_absvmag),np.median(tuc2_mlratio),xerr=np.std(tuc2_absvmag),yerr=[tuc2_lower_error,tuc2_upper_error],elinewidth=0.5,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
ax2_0.errorbar(np.median(gru1_absvmag),np.median(gru1_mlratio),xerr=np.std(gru1_absvmag),yerr=[gru1_lower_error,gru1_upper_error],elinewidth=0.5,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
ax2_0.errorbar(np.median(ret2_absvmag),np.median(ret2_mlratio),xerr=np.std(ret2_absvmag),yerr=[ret2_lower_error,ret2_upper_error],elinewidth=0.5,fmt='.',capsize=0,alpha=1,color='k',rasterized=True)
#ax2_0.errorbar([gru1_absvmag0],[gru1_mlratio0],xerr=[gru1_sigabsvmag0],yerr=[gru1_sigmlratio0],elinewidth=0.25,fmt='.',capsize=0,alpha=1,color='g',rasterized=True)
#ax2_0.text(gru1_absvmag[0]-0.15,gru1_mlratio[0]+45,'Ret 2',fontsize=7,color='k')
ax2_0.legend(loc=1,fontsize=5,handlelength=0,numpoints=1,scatterpoints=1,shadow=False)

#maxlikemdyn=1./0.0043*gru1_rhalf0*gru1vdisp[(gru1like == max(gru1like))[0]]**2
#maxlikemrhalf=2.5/0.0043*gru1_rhalf0*gru1vdisp[(gru1like == max(gru1like))[0]]**2
#maxlikemldyn=maxlikemdyn/gru1_luminosity0
#maxlikemlratio=maxlikemrhalf/(gru1_luminosity0/2.)
#lower_error=np.array(maxlikemldyn-mycode2.hpd(gru1_mlratio,0.32)[0])
#upper_error=np.array(mycode2.hpd(gru1_mlratio,0.32)[1]-maxlikemldyn)

g1.write(r'\newcommand{\grumrhalf}{$'+str(round(np.median(2.5*gru1_mrhalf/1.e+6),1))+r'_{'+str(round(np.percentile(2.5*gru1_mrhalf/1.e+6,16)-np.percentile(2.5*gru1_mrhalf/1.e+6,50),1))+r'}^{+'+str(round(np.percentile(2.5*gru1_mrhalf/1.e+6,84)-np.percentile(2.5*gru1_mrhalf/1.e+6,50),1))+r'}\times 10^6$}'+'\n')
g1.write(r'\newcommand{\grumlratio}{$'+str(int(np.median(2.*2.5*gru1_mlratio)))+r'_{'+str(int(np.percentile(2.*2.5*gru1_mlratio,16)-np.percentile(2.*2.5*gru1_mlratio,50)))+r'}^{+'+str(int(np.percentile(2.*2.5*gru1_mlratio,84)-np.percentile(2.*2.5*gru1_mlratio,50)))+r'}$}'+'\n')

g1.write(r'\newcommand{\grumrhalfupperlim}{$'+str(round(np.percentile(2.5*gru1_mrhalf/1.e+6,84),1))+r'\times 10^6$}'+'\n')
g1.write(r'\newcommand{\grumlratioupperlim}{$'+str(int(np.percentile(2.*2.5*gru1_mlratio,84)))+r'$}'+'\n')

g1.write(r'\newcommand{\tucmrhalfupperlim}{$'+str(round(np.percentile(2.5*tuc2_mrhalf/1.e+6,84),1))+r'\times 10^6$}'+'\n')
g1.write(r'\newcommand{\tucmlratioupperlim}{$'+str(int(np.percentile(2.*2.5*tuc2_mlratio,84)))+r'$}'+'\n')

g1.write(r'\newcommand{\tucmrhalf}{$'+str(round(np.median(2.5*tuc2_mrhalf/1.e+6),1))+r'_{'+str(round(np.percentile(2.5*tuc2_mrhalf/1.e+6,16)-np.percentile(2.5*tuc2_mrhalf/1.e+6,50),1))+r'}^{+'+str(round(np.percentile(2.5*tuc2_mrhalf/1.e+6,84)-np.percentile(2.5*tuc2_mrhalf/1.e+6,50),1))+r'}\times 10^6$}'+'\n')
g1.write(r'\newcommand{\tucmlratio}{$'+str(int(np.median(2.*2.5*tuc2_mlratio)))+r'_{'+str(int(np.percentile(2.*2.5*tuc2_mlratio,16)-np.percentile(2.*2.5*tuc2_mlratio,50)))+r'}^{+'+str(int(np.percentile(2.*2.5*tuc2_mlratio,84)-np.percentile(2.*2.5*tuc2_mlratio,50)))+r'}$}'+'\n')
g1.write(r'\newcommand{\tucmlratiotwosigma}{$'+str(int(np.median(2.*2.5*tuc2_mlratio)))+r'_{'+str(int(np.percentile(2.*2.5*tuc2_mlratio,2.5)-np.percentile(2.*2.5*tuc2_mlratio,50)))+r'}^{+'+str(int(np.percentile(2.*2.5*tuc2_mlratio,97.5)-np.percentile(2.*2.5*tuc2_mlratio,50)))+r'}$}'+'\n')

g1.write(r'\newcommand{\retmrhalf}{$'+str(round(np.median(2.5*ret2_mrhalf/1.e+6),1))+r'_{'+str(round(np.percentile(2.5*ret2_mrhalf/1.e+6,16)-np.percentile(2.5*ret2_mrhalf/1.e+6,50),1))+r'}^{+'+str(round(np.percentile(2.5*ret2_mrhalf/1.e+6,84)-np.percentile(2.5*ret2_mrhalf/1.e+6,50),1))+r'}\times 10^6$}'+'\n')
g1.write(r'\newcommand{\retmlratio}{$'+str(int(np.median(2.*2.5*ret2_mlratio)))+r'_{'+str(int(np.percentile(2.*2.5*ret2_mlratio,16)-np.percentile(2.*2.5*ret2_mlratio,50)))+r'}^{+'+str(int(np.percentile(2.*2.5*ret2_mlratio,84)-np.percentile(2.*2.5*ret2_mlratio,50)))+r'}$}'+'\n')

g1.close()

plotfilename='tuc2gru1teffprior_scaling.pdf'
plt.savefig(plotfilename,dpi=400)
plt.show()
plt.close()
