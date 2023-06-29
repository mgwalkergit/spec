import numpy as np
import astropy
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import StdDevUncertainty
from astropy.nddata import CCDData
import dustmaps.sfd
import astropy.units as u
from astropy.modeling import models
import os
from os import path
import scipy
import mycode
import m2fs_process as m2fs
import crossmatcher
#matplotlib.use('TkAgg')
import dill as pickle
from scipy.stats import uniform,norm
from pymultinest.solve import solve
from isochrones.mist import MIST_Isochrone
from isochrones.mist import MIST_EvolutionTrack
from isochrones.mist.bc import MISTBolometricCorrectionGrid
from isochrones import get_ichrone
from isochrones.priors import GaussianPrior,SalpeterPrior,DistancePrior,FlatPrior,ChabrierPrior
from isochrones.populations import StarFormationHistory,StarPopulation
from dynesty import DynamicNestedSampler,NestedSampler
import minimint


err_mag_v,err_err_v,err_std_v=np.load('hst_empirical_errors_v.npy')
err_mag_i,err_err_i,err_std_i=np.load('hst_empirical_errors_i.npy')

n1=100000
n2=10000

mist=get_ichrone('mist')
#imf=SalpeterPrior(bounds=(0.08,10))
imf=ChabrierPrior(bounds=(0.08,10))
fb=0.4
gamma=0.3
sfh=StarFormationHistory(dist=uniform(np.log10(1.e+10),np.log10(1.2e+10)))
sfh2=StarFormationHistory(dist=uniform(np.log10(1.e+9),np.log10(1.2e+9)))
feh=GaussianPrior(-1.,0.5)
feh2=GaussianPrior(-1.,0.01)
distance=FlatPrior(bounds=[138000,142000])
av=FlatPrior(bounds=[0.01,0.05])
pop=StarPopulation(mist,imf=imf,fB=fb,gamma=gamma,sfh=sfh,feh=feh,distance=distance,AV=av)
pop2=StarPopulation(mist,imf=imf,fB=fb,gamma=gamma,sfh=sfh2,feh=feh2,distance=distance,AV=av)
df=pop.generate(n1,exact_N=True)
df2=pop2.generate(n2,exact_N=True)

vmag_err1=np.random.normal(loc=0.,scale=1.,size=n1)
vmag_err2=np.random.normal(loc=0.,scale=1.,size=n2)
imag_err1=np.random.normal(loc=0.,scale=1.,size=n1)
imag_err2=np.random.normal(loc=0.,scale=1.,size=n2)
                           
mbol=df['Mbol_0'].array
eep=df['eep_0'].array
initial_mass=df['initial_mass_0'].array
feh=df['feh_0'].array
teff=df['Teff_0'].array
logg=df['logg_0'].array
mass=df['mass_0'].array
gmag=df['G_mag'].array
bpmag=df['BP_mag'].array
rpmag=df['RP_mag'].array
age=df['age_0'].array
av=df['AV_0'].array
ag=df['A_G'].array
distance=df['distance_0'].array
dmodulus=5.*np.log10(distance)-5.

mbol2=df2['Mbol_0'].array
eep2=df2['eep_0'].array
initial_mass2=df2['initial_mass_0'].array
feh2=df2['feh_0'].array
teff2=df2['Teff_0'].array
logg2=df2['logg_0'].array
mass2=df2['mass_0'].array
gmag2=df2['G_mag'].array
bpmag2=df2['BP_mag'].array
rpmag2=df2['RP_mag'].array
age2=df2['age_0'].array
av2=df2['AV_0'].array
ag2=df2['A_G'].array
distance2=df2['distance_0'].array
dmodulus2=5.*np.log10(distance2)-5.

mist=MIST_Isochrone()
mist_track=MIST_EvolutionTrack()
bc_grid=MISTBolometricCorrectionGrid(['G','BP','RP','DECam_g','DECam_r','DECam_i','DECam_z','DECam_Y','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','PS_g','PS_r','PS_i','PS_z','F150W','F090W','WFC3_UVIS_F555W','WFC3_UVIS_F606W','WFC3_UVIS_F814W'])

piss=bc_grid.interp([teff,logg,feh,av],['WFC3_UVIS_F555W','WFC3_UVIS_F606W','WFC3_UVIS_F814W'])
piss2=bc_grid.interp([teff2,logg2,feh2,av2],['WFC3_UVIS_F555W','WFC3_UVIS_F606W','WFC3_UVIS_F814W'])

bc_f555w=piss.T[0]
bc_f606w=piss.T[1]
bc_f814w=piss.T[2]
f555w=-bc_f555w+mbol+dmodulus
f606w=-bc_f606w+mbol+dmodulus
f814w=-bc_f814w+mbol+dmodulus

f606w_obs=[]
f814w_obs=[]
for i in range(0,len(f606w)):
    mean=np.interp(f606w[i],err_mag_v,err_err_v)
    std=np.interp(f606w[i],err_mag_v,err_std_v)
    f606w_obs.append(f606w[i]+mean*vmag_err1[i]/5)
    mean=np.interp(f814w[i],err_mag_i,err_err_i)
    std=np.interp(f814w[i],err_mag_i,err_std_i)
    f814w_obs.append(f814w[i]+mean*imag_err1[i]/5)
f606w_obs=np.array(f606w_obs)
f814w_obs=np.array(f814w_obs)
    
bc_f555w2=piss2.T[0]
bc_f606w2=piss2.T[1]
bc_f814w2=piss2.T[2]
f555w2=-bc_f555w2+mbol2+dmodulus2
f606w2=-bc_f606w2+mbol2+dmodulus2
f814w2=-bc_f814w2+mbol2+dmodulus2

f606w2_obs=[]
f814w2_obs=[]
for i in range(0,len(f606w2)):
    mean=np.interp(f606w2[i],err_mag_v,err_err_v)
    std=np.interp(f606w2[i],err_mag_v,err_std_v)
    f606w2_obs.append(f606w2[i]+mean*vmag_err2[i]/5)
    mean=np.interp(f814w2[i],err_mag_i,err_err_i)
    std=np.interp(f814w2[i],err_mag_i,err_std_i)
    f814w2_obs.append(f814w2[i]+mean*imag_err2[i]/5)
f606w2_obs=np.array(f606w2_obs)
f814w2_obs=np.array(f814w2_obs)


col=f606w_obs-f814w_obs
mag=f606w
col2=f606w2_obs-f814w2_obs
mag2=f606w2
plt.scatter(col,mag,s=1,alpha=0.2,color='k',rasterized=True)
plt.scatter(col2,mag2,s=1,alpha=0.2,color='r',rasterized=True)
plt.xlim([-1,3])
plt.ylim([28,16])
plt.xlabel('F606W-F814W')
plt.ylabel('F606W')
plt.savefig('fornax_test.pdf',dpi=200)
plt.show()
plt.close()

np.pause()



#ii=minimint.Interpolator(['Gaia_G_EDR3','Gaia_BP_EDR3','Gaia_RP_EDR3','DECam_g','DECam_r','DECam_i','DECam_z','PS_g','PS_r','PS_i','PS_z','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','HST_ACS_HR_F555W'])

imf_alpha1=0.3
imf_alpha2=1.3
imf_alpha3=2.3
imf_break1=0.08
imf_break2=0.5
imf_min_mass=0.08
imf_max_mass=10.#np.inf
imf_part1=(imf_break1**(1.-imf_alpha1)-imf_min_mass**(1.-imf_alpha1))/(1.-imf_alpha1)
imf_part2=imf_break1**(imf_alpha2-imf_alpha1)*(imf_break2**(1.-imf_alpha2)-imf_break1**(1.-imf_alpha2))/(1.-imf_alpha2)
imf_part3=imf_break1**(imf_alpha2-imf_alpha1)*imf_break2**(imf_alpha3-imf_alpha2)*(imf_max_mass**(1.-imf_alpha3)-imf_break2**(1.-imf_alpha3))/(1.-imf_alpha3)
imf_const=imf_part1+imf_part2+imf_part3

ran=np.random.uniform(low=0.,high=1.,size=1000000)
mass=[]
for i in range(0,len(ran)):
    if ran[i]<=imf_part1/imf_const:
        x=((1.-imf_alpha1)*imf_const*ran[i]+imf_min_mass**(1.-imf_alpha1))**(1./(1.-imf_alpha1))
    elif ((ran[i]>imf_part1/imf_const)&(ran[i]<=(imf_part1+imf_part2)/imf_const)):
        x=((1.-imf_alpha2)*(imf_const*ran[i]-imf_part1)*imf_break1**(imf_alpha1-imf_alpha2)+imf_break1**(1.-imf_alpha2))**(1./(1.-imf_alpha2))
    elif ran[i]>(imf_part1+imf_part2)/imf_const:
        x=((1.-imf_alpha3)*(imf_const*ran[i]-imf_part1-imf_part2)*imf_break1**(imf_alpha1-imf_alpha2)*imf_break2**(imf_alpha2-imf_alpha3)+imf_break2**(1.-imf_alpha3))**(1./(1.-imf_alpha3))
    else:
        raise ValueError('something wrong in sampling Kroupa IMF')
    mass.append(x)
mass=np.array(mass)
#matplotlib.use('pdf')

mist=MIST_Isochrone()
mist_track=MIST_EvolutionTrack()
bc_grid=MISTBolometricCorrectionGrid(['G','BP','RP','DECam_g','DECam_r','DECam_i','DECam_z','DECam_Y','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','PS_g','PS_r','PS_i','PS_z','F150W','F090W','WFC3_UVIS_F555W','WFC3_UVIS_F606W','WFC3_UVIS_F814W'])


feh=np.ones(len(mass))*(0.1)
logage=np.ones(len(mass))*np.log10(1.e+10)
av=np.ones(len(mass))*0.

shite=mist_track.generate(mass,logage,feh,accruate=True)[['eep','logg','Teff','Mbol']]
#piss=bc_grid.interp([shite.Teff.array[0],shite.logg.array[0],feh,av],['WFC3_UVIS_F555W'])
piss=bc_grid.interp([shite.Teff.array,shite.logg.array,feh,av],['WFC3_UVIS_F555W','WFC3_UVIS_F606W','WFC3_UVIS_F814W'])

mbol=shite.Mbol.array
bc_f555w=piss.T[0]
bc_f606w=piss.T[1]
bc_f814w=piss.T[2]

f555w=-bc_f555w+mbol
f606w=-bc_f606w+mbol
f814w=-bc_f814w+mbol

dmodulus=20.
col=f606w-f814w
mag=f606w+dmodulus
plt.scatter(col,mag,s=1)
plt.xlim([0,2])
plt.ylim([28,20])
plt.show()
plt.close()
np.pause()

def myprior(cube):
    from scipy.special import erfinv
    from statistics import NormalDist
    prior=[]
    prior.append([0.,0.])#stellar mass
    prior.append([200,808])#EEP
    prior.append([-4.,+0.5])#[Fe/H]
    prior.append([0.,300.])# distance / kpc
    prior.append([0.,0.])# reddening Av
    prior=np.array(prior)
    x=np.array(cube)
    #draw stellar mass from Kroupa IMF
    if cube[0]<=imf_part1/imf_const:
        x[0]=((1.-imf_alpha1)*imf_const*cube[0]+imf_min_mass**(1.-imf_alpha1))**(1./(1.-imf_alpha1))
    elif ((cube[0]>imf_part1/imf_const)&(cube[0]<=(imf_part1+imf_part2)/imf_const)):
        x[0]=((1.-imf_alpha2)*(imf_const*cube[0]-imf_part1)*imf_break1**(imf_alpha1-imf_alpha2)+imf_break1**(1.-imf_alpha2))**(1./(1.-imf_alpha2))
    elif cube[0]>(imf_part1+imf_part2)/imf_const:
        x[0]=((1.-imf_alpha3)*(imf_const*cube[0]-imf_part1-imf_part2)*imf_break1**(imf_alpha1-imf_alpha2)*imf_break2**(imf_alpha2-imf_alpha3)+imf_break2**(1.-imf_alpha3))**(1./(1.-imf_alpha3))
    else:
        raise ValueError('something wrong in sampling Kroupa IMF')
    #Gaussian prior on A_V, centered on dustmap value, with sigma = 15% of dustmap value
    av_mean=ebv[j]*2.742#reddening at R_V=3.1 value from Schlafly & Finkbeiner 2011
    av_sigma=0.15*av_mean
    x[4]=np.max(np.array([NormalDist(mu=av_mean,sigma=av_sigma).inv_cdf(cube[4])]))
    #for the rest of parameters, use flat priors
    for i in range(1,len(x)-1):
        x[i]=prior[i][0]+(prior[i][1]-prior[i][0])*cube[i]
    return x
            
def myloglike(cube):
    mass=cube[0]
    eep=cube[1]
    feh=cube[2]
    distance=cube[3]*1000.
    av=cube[4]
    pars=[mass,eep,feh]
    Mbol,teff,logg,feh=mist_track.interp_value(pars,['Mbol','Teff','logg','feh'])
    bolometric_correction=bc_grid.interp([teff,logg,feh,av],['G','BP','RP','DECam_g','DECam_r','DECam_i','DECam_z','DECam_Y','DECam_g','DECam_r','DECam_z','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','PS_g','PS_r','PS_i','PS_z'])
    dmodulus=5.*np.log10(distance)-5.
    absmags=Mbol-bolometric_correction
    iso_model_mags=absmags+dmodulus
    sum1=-0.5*np.log(2.*np.pi)*len(keep)
    sum2=np.sum(-0.5*np.log(sigmags[keep]**2))
    sum3=np.sum(-0.5*(mags[keep]-iso_model_mags[keep])**2/(sigmags[keep]**2))
    if sum1+sum2+sum3!=sum1+sum2+sum3:
        return -1.e+30
    else:
        return sum1+sum2+sum3
        
    parameters=['mass','eep','feh','distance','av']
    n_params=len(parameters)
    prefix="chains/poop-"
    #            result=solve(LogLikelihood=myloglike,Prior=myprior,n_dims=n_params,outputfiles_basename=prefix,n_live_points=1000,verbose=True,resume=False)
    sampler=DynamicNestedSampler(myloglike,myprior,n_params,bound='multi')
    #            sampler=NestedSampler(myloglike,myprior,n_params,bound='multi',nlive=400)
    sampler.run_nested(dlogz_init=0.5,nlive_init=500,wt_kwargs={'pfrac':1.0})#dlogz=0.5,maxiter=10000,maxcall=50000)
    res=sampler.results
    pickle.dump(res,open('dynesty.res','wb'))
    mass=result['samples'].T[0]
    eep=result['samples'].T[1]
    distance=result['samples'].T[3]
    teff=[]
    logg=[]
    feh=[]
    for j in range(0,len(result['samples'].T[0])):
        teff0,logg0,feh0=mist_track.interp_value([result['samples'].T[0][j],result['samples'].T[1][j],result['samples'].T[2][j]],['Teff','logg','feh'])
        teff.append(teff0)
        logg.append(logg0)
        feh.append(feh0)
    teff=np.array(teff)
    logg=np.array(logg)
    feh=np.array(feh)
    print(np.median(mass),np.std(mass))
    print(np.median(teff),np.std(teff))
    print(np.median(logg),np.std(logg))
    print(np.median(feh),np.std(feh))
    print(np.median(distance),np.std(distance))
    print(np.median(eep),np.std(eep))
    #            pars[np.median(mass),np.median(]
    #            shite=mist_track.interp_mag(pars,['G','BP','RP'])
    plt.ylim([23,16])
    plt.xlim([-1,3])
    plt.show()
    plt.close()

            
