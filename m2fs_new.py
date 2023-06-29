import numpy as np
import astropy
import sqlutil
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.nddata import NDData
from astropy.nddata import StdDevUncertainty
from astropy.nddata import CCDData
from astropy.coordinates import SkyCoord
import dustmaps.sfd
import astropy.units as u
from astropy.modeling import models
import os
import csv
from os import path
import scipy
import m2fs_process as m2fs
import mycode
import crossmatcher
matplotlib.use('TkAgg')
from matplotlib.patches import Ellipse
import dill as pickle
from isochrones.mist import MIST_Isochrone
from isochrones.mist import MIST_EvolutionTrack
from isochrones.mist.bc import MISTBolometricCorrectionGrid
from pymultinest.solve import solve
from pymultinest import Analyzer
#matplotlib.use('pdf')

get_surveyphot=False#pull photometry for each spectroscopic target from various sky surveys
get_obstable=False#write observation log in latex table
make_skysub=False#generate 1d sky-subtracted spectra?
get_fits=False#fit is done, do post-processing of multinest output
get_catalog_raw=False#catalog post-processed data
get_errors=False#fit error models using repeat observations
get_catalog_final=False#update catalog for calibrated errors
check_badfield=False
get_plots=True

data_directory='/hildafs/projects/phy200028p/mgwalker/m2fs_data/'
chains_directory='/hildafs/projects/phy200028p/mgwalker/m2fs_chains/'
fit_directory='/hildafs/projects/phy200028p/mgwalker/m2fs_chains/'
lambdamin=5127.
lambdamax=5190.

fits_list0='new_m2fshiresian_files'
fits_list='/hildafs/projects/phy200028p/mgwalker/m2fs/'+fits_list0
fits_table_raw_filename=fits_list.split('_files')[0]+'_raw.fits'
fits_table_final_filename='/hildafs/projects/phy200028p/mgwalker/m2fs/m2fs_fits_table.fits'
data_out=fits_list.split('_files')[0]+'.dat'
data_filename='/hildafs/projects/phy200028p/mgwalker/MultiNest_v2.17/'+fits_list0.split('_files')[0]+'_data'
chains_filename='/hildafs/projects/phy200028p/mgwalker/MultiNest_v2.17/'+fits_list0.split('_files')[0]+'_chains'

mist=MIST_Isochrone()
mist_track=MIST_EvolutionTrack()
bc_grid=MISTBolometricCorrectionGrid(['G','BP','RP','DECam_g','DECam_r','DECam_i','DECam_z','DECam_Y','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','PS_g','PS_r','PS_i','PS_z'])

des_est50,des_est16,des_est84=pickle.load(open('/hildafs/projects/phy200028p/mgwalker/gaia_dwarfs/photometry/des_ensemble.pickle','rb'))
decals_est50,decals_est16,decals_est84=pickle.load(open('/hildafs/projects/phy200028p/mgwalker/gaia_dwarfs/photometry/decals_ensemble.pickle','rb'))
gaia_est50,gaia_est16,gaia_est84=pickle.load(open('/hildafs/projects/phy200028p/mgwalker/gaia_dwarfs/photometry/gaia_ensemble.pickle','rb'))
ps1_est50,ps1_est16,ps1_est84=pickle.load(open('/hildafs/projects/phy200028p/mgwalker/gaia_dwarfs/photometry/ps1_ensemble.pickle','rb'))
sdss_est50,sdss_est16,sdss_est84=pickle.load(open('/hildafs/projects/phy200028p/mgwalker/gaia_dwarfs/photometry/sdss_ensemble.pickle','rb'))

thing=['vlos','teffprior_vlos','teff','teffprior_teff','logg','teffprior_logg','feh','teffprior_feh','mgfe','teffprior_mgfe','smooth','teffprior_smooth']
thing_symbol=[r'$v_{\rm los}$',r'$v_{\rm los}$',r'$T_{\rm eff}$',r'$T_{\rm eff}$',r'$\log_{10}[g]$',r'$\log_{10}[g]$',r'[Fe/H]',r'[Fe/H]',r'[Mg/Fe]',r'[Mg/Fe]',r'$\sigma_{\rm LSF}$',r'$\sigma_{\rm LSF}$']

class teffpriorstuff:
    def __init__(self,gaia_teffprior=None,des_teffprior=None,decals_teffprior=None,sdss_teffprior=None,ps1_teffprior=None):
        self.gaia_teffprior=gaia_teffprior
        self.des_teffprior=des_teffprior
        self.decals_teffprior=decals_teffprior
        self.sdss_teffprior=sdss_teffprior
        self.ps1_teffprior=ps1_teffprior

with open(fits_list) as f:
    data=f.readlines()
fitsfile=[]
obj=[]
for line in data:
    p=line.split()
    fitsfile.append(p[0])
    obj.append(p[1])
fitsfile=np.array(fitsfile)
obj=np.array(obj)
for i in range(0,len(fitsfile)):
    this=np.where(fitsfile==fitsfile[i])[0]
    if len(this)>1:
        print('ERROR: multiple listings of same observation in input file!!!')
        print(this)
        np.pause()

if get_surveyphot:

    rrl_id=sqlutil.get('''select source_id from gaia_dr3.vari_rrlyrae''',db='wsdb',host='wsdb.hpc1.cs.cmu.edu',user='matt_walker',password='5UN*ur8Y')[0]#pull gaia source IDs for all stars in RRL catalog
    for i in range(0,len(fitsfile)):
        print(i,' of ',len(fitsfile))
        hdul=fits.open(fitsfile[i])
        fitsobject=m2fs.m2fs_getfromfits(hdul)
        filtername=hdul[0].header['filtername']
        m2fsrun=hdul[0].header['m2fsrun']
        field_name=hdul[0].header['field_name']
        temperature=hdul[0].header['CCD_temp']
        hdul.close()
        surveyphot_filename=fitsfile[i].split('_stackskysub.fits')[0]+'_surveyphot.pkl'
    
        ebv=[]
        for j in range(0,len(fitsobject.radeg)):
            if fitsobject.radeg[j]>-998.:
                ebv.append(dustmaps.sfd.SFDQuery().query_equ(fitsobject.radeg[j],fitsobject.decdeg[j]))
            else:
                ebv.append(-999.)
        ebv=np.array(ebv)

        gaia_objid,gaia_ra,gaia_dec,gaia_gmag,gaia_bpmag,gaia_rpmag,gaia_gflux,gaia_bpflux,gaia_rpflux,gaia_siggflux,gaia_sigbpflux,gaia_sigrpflux,gaia_parallax,gaia_sigparallax,gaia_pmra,gaia_sigpmra,gaia_pmdec,gaia_sigpmdec=crossmatcher.doit('gaia_edr3.gaia_source',fitsobject.radeg,fitsobject.decdeg,'source_id,ra,dec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,phot_g_mean_flux,phot_bp_mean_flux,phot_rp_mean_flux,phot_g_mean_flux_error,phot_bp_mean_flux_error,phot_rp_mean_flux_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error',rad=1.,db='wsdb',host='wsdb.hpc1.cs.cmu.edu',user='matt_walker',password='5UN*ur8Y')
        gaia_rrl=np.zeros(len(gaia_objid),dtype='int')
        for j in range(0,len(gaia_objid)):
            this=np.where(rrl_id==gaia_objid[j])[0]
            if len(this)>0:
                gaia_rrl[j]=1
        des_gmag,des_siggmag,des_rmag,des_sigrmag,des_imag,des_sigimag,des_zmag,des_sigzmag,des_ymag,des_sigymag,des_ebv,des_gspread,des_rspread,des_ispread,des_zspread,des_yspread=crossmatcher.doit('des_dr2.main',fitsobject.radeg,fitsobject.decdeg,'wavg_mag_psf_g,wavg_magerr_psf_g,wavg_mag_psf_r,wavg_magerr_psf_r,wavg_mag_psf_i,wavg_magerr_psf_i,wavg_mag_psf_z,wavg_magerr_psf_z,wavg_mag_psf_y,wavg_magerr_psf_y,ebv_sfd98,wavg_spread_model_g,wavg_spread_model_r,wavg_spread_model_i,wavg_spread_model_z,wavg_spread_model_y',extra='wavg_magerr_psf_g>=0.',rad=1.,db='wsdb',host='wsdb.hpc1.cs.cmu.edu',user='matt_walker',password='5UN*ur8Y')
        decals_gflux,decals_ivargflux,decals_rflux,decals_ivarrflux,decals_zflux,decals_ivarzflux,decals_type=crossmatcher.doit('decals_dr9.main',fitsobject.radeg,fitsobject.decdeg,'flux_g,flux_ivar_g,flux_r,flux_ivar_r,flux_z,flux_ivar_z,type',extra='flux_g>=0',rad=1.,db='wsdb',host='wsdb.hpc1.cs.cmu.edu',user='matt_walker',password='5UN*ur8Y')
        decals_gmag=22.5-2.5*np.log10(decals_gflux)
        decals_rmag=22.5-2.5*np.log10(decals_rflux)
        decals_zmag=22.5-2.5*np.log10(decals_zflux)
        decals_siggmag=2.5/decals_gflux/np.log(10.)/np.sqrt(decals_ivargflux)
        decals_sigrmag=2.5/decals_rflux/np.log(10.)/np.sqrt(decals_ivarrflux)
        decals_sigzmag=2.5/decals_zflux/np.log(10.)/np.sqrt(decals_ivarzflux)    
        ps1_objid,ps1_gmag,ps1_rmag,ps1_imag,ps1_zmag,ps1_siggmag,ps1_sigrmag,ps1_sigimag,ps1_sigzmag,ps1_ebv,ps1_rkronmag=crossmatcher.doit('panstarrs_dr1.stackobjectthin',fitsobject.radeg,fitsobject.decdeg,'objid,gpsfmag,rpsfmag,ipsfmag,zpsfmag,gpsfmagerr,rpsfmagerr,ipsfmagerr,zpsfmagerr,ebv,rkronmag',extra='(rinfoflag3& panstarrs_dr1.detectionflags3(\'STACK_PRIMARY\'))>0 and gpsfmag=gpsfmag',rad=1.,db='wsdb',host='wsdb.hpc1.cs.cmu.edu',user='matt_walker',password='5UN*ur8Y')
        sdss_objid,sdss_umag,sdss_gmag,sdss_rmag,sdss_imag,sdss_zmag,sdss_sigumag,sdss_siggmag,sdss_sigrmag,sdss_sigimag,sdss_sigzmag,sdss_extinction_u,sdss_extinction_g,sdss_extinction_r,sdss_extinction_i,sdss_extinction_z,sdss_type,sdss_mode=crossmatcher.doit('sdssdr9.phototag',fitsobject.radeg,fitsobject.decdeg,'objid,psfmag_u,psfmag_g,psfmag_r,psfmag_i,psfmag_z,psfmagerr_u,psfmagerr_g,psfmagerr_r,psfmagerr_i,psfmagerr_z,extinction_u,extinction_g,extinction_r,extinction_i,extinction_z,type,mode',extra='psfmag_g>-9998.',rad=1.,db='wsdb',host='wsdb.hpc1.cs.cmu.edu',user='matt_walker',password='5UN*ur8Y')

        sdss_uext=ebv*4.239#de-reddened mag will be the psfmag minus this value
        sdss_gext=ebv*3.303
        sdss_rext=ebv*2.285
        sdss_iext=ebv*1.698
        sdss_zext=ebv*1.263
        sdss=m2fs.survey_phot_sdss(sdss_objid,sdss_umag,sdss_gmag,sdss_rmag,sdss_imag,sdss_zmag,sdss_sigumag,sdss_siggmag,sdss_sigrmag,sdss_sigimag,sdss_sigzmag,sdss_type,sdss_mode,ebv,sdss_uext,sdss_gext,sdss_rext,sdss_iext,sdss_zext)
    
        des_gext=ebv*3.237#3.214
        des_rext=ebv*2.176#2.165
        des_iext=ebv*1.595#1.592
        des_zext=ebv*1.217#1.211
        des_yext=ebv*1.058#1.064#de-reddened mag will be the psfmag minus this value
        des=m2fs.survey_phot_des(des_gmag,des_rmag,des_imag,des_zmag,des_ymag,des_siggmag,des_sigrmag,des_sigimag,des_sigzmag,des_sigymag,des_gspread,des_rspread,des_ispread,des_zspread,des_yspread,ebv,des_gext,des_rext,des_iext,des_zext,des_yext)
        
        decals_gext=ebv*3.237#3.214
        decals_rext=ebv*2.176#2.165
        decals_zext=ebv*1.217#1.211
        decals=m2fs.survey_phot_decals(decals_gmag,decals_rmag,decals_zmag,decals_siggmag,decals_sigrmag,decals_sigzmag,decals_type,ebv,decals_gext,decals_rext,decals_zext)
    
        ps1_gext=ebv*3.172
        ps1_rext=ebv*2.271
        ps1_iext=ebv*1.682
        ps1_zext=ebv*1.322
        ps1=m2fs.survey_phot_ps1(ps1_objid,ps1_gmag,ps1_rmag,ps1_imag,ps1_zmag,ps1_siggmag,ps1_sigrmag,ps1_sigimag,ps1_sigzmag,ps1_rkronmag,ebv,ps1_gext,ps1_rext,ps1_iext,ps1_zext)

        gaia_siggmag=(2.5/gaia_gflux/2.30258)**2*gaia_siggflux**2
        gaia_sigbpmag=(2.5/gaia_bpflux/2.30258)**2*gaia_sigbpflux**2
        gaia_sigrpmag=(2.5/gaia_rpflux/2.30258)**2*gaia_sigrpflux**2
        a0=3.1*ebv*0.86#use Eq. 1 from Babusiaux et al 2018 (1804.09378) to get extinction coefficients, apply extinction to update BP-RP, then interate again to get final estimate of extinction coefficient.  The factor of 0.86 comes from Sergey, as explained to me in email from Andrew on July 11, 2022: "I've recently realized that Babusieux's formula was written to be applied to "True" EBV. But the SFD EBV is overestimated according to Schlafly11. So if you use SFD EBV's for the getDust, then the ebv needs to be multiplied by 0.86. (when you work with other surveys it's usually not a problem, because the coefficients are already given wrt SFD, but that I learned doesn't apply for gaia's code."
        kg_c=np.array([0.9761,-0.1704,0.0086,0.0011,-0.0438,0.0013,0.0099])
        kbp_c=np.array([1.1517,-0.0871,-0.0333,0.0173,-0.0230,0.0006,0.0043])
        krp_c=np.array([0.6104,-0.0170,-0.0026,-0.0017,-0.0078,0.00005,0.0006])
        col0=gaia_bpmag-gaia_rpmag
        kg=kg_c[0]+kg_c[1]*col0+kg_c[2]*col0**2+kg_c[3]*col0**3+kg_c[4]*a0+kg_c[5]*a0**2+kg_c[6]*col0*a0
        kbp=kbp_c[0]+kbp_c[1]*col0+kbp_c[2]*col0**2+kbp_c[3]*col0**3+kbp_c[4]*a0+kbp_c[5]*a0**2+kbp_c[6]*col0*a0
        krp=krp_c[0]+krp_c[1]*col0+krp_c[2]*col0**2+krp_c[3]*col0**3+krp_c[4]*a0+krp_c[5]*a0**2+krp_c[6]*col0*a0
        ag=kg*a0
        abp=kbp*a0
        arp=krp*a0
        col1=gaia_bpmag-abp-(gaia_rpmag-arp)
        kg=kg_c[0]+kg_c[1]*col1+kg_c[2]*col1**2+kg_c[3]*col1**3+kg_c[4]*a0+kg_c[5]*a0**2+kg_c[6]*col1*a0
        kbp=kbp_c[0]+kbp_c[1]*col1+kbp_c[2]*col1**2+kbp_c[3]*col1**3+kbp_c[4]*a0+kbp_c[5]*a0**2+kbp_c[6]*col1*a0
        krp=krp_c[0]+krp_c[1]*col1+krp_c[2]*col1**2+krp_c[3]*col1**3+krp_c[4]*a0+krp_c[5]*a0**2+krp_c[6]*col1*a0
        gaia_ag=kg*a0
        gaia_abp=kbp*a0
        gaia_arp=krp*a0
        gaia_gmag_dered=gaia_gmag-gaia_ag
        gaia_bpmag_dered=gaia_bpmag-gaia_abp
        gaia_rpmag_dered=gaia_rpmag-gaia_arp
        
        gaia=m2fs.survey_phot_gaia(gaia_objid,gaia_ra,gaia_dec,gaia_gmag,gaia_bpmag,gaia_rpmag,gaia_siggmag,gaia_sigbpmag,gaia_sigrpmag,ebv,gaia_ag,gaia_abp,gaia_arp,gaia_pmra,gaia_sigpmra,gaia_pmdec,gaia_sigpmdec,gaia_parallax,gaia_sigparallax,gaia_rrl)
        survey_phot=m2fs.survey_phot(gaia,des,decals,sdss,ps1)
        print(i,surveyphot_filename)
        pickle.dump(survey_phot,open(surveyphot_filename,'wb'))

if get_obstable:

    field_center=[]
    field_id=[]
    field_utdate=[]
    field_uttime=[]
    field_exptime=[]
    field_object=[]
    field_ntarget=[]
    field_nsky=[]
    field_subframes=[]

    for i in range(0,len(fitsfile)):
        print(i,' of ',len(fitsfile))
        hdul=fits.open(fitsfile[i])
        fitsobject=m2fs.m2fs_getfromfits(hdul)
        filtername=hdul[0].header['filtername']
        m2fsrun=hdul[0].header['m2fsrun']
        field_name=hdul[0].header['field_name']
        temperature=hdul[0].header['CCD_temp']
        hdul.close()
        surveyphot_filename=fitsfile[i].split('_stackskysub.fits')[0]+'_surveyphot.pkl'
    
        skies=np.where(fitsobject.obj=='SKY')[0]
        targets=np.where(fitsobject.icode>0)[0]

        field_center.append(SkyCoord(hdul[0].header['RA']+hdul[0].header['DEC'],unit=(u.hourangle,u.deg)))
        field_id.append(hdul[0].header['field_name'])
        field_utdate.append(hdul[0].header['UT-DATE'])
        field_uttime.append(hdul[0].header['UT-TIME'])
        field_exptime.append(hdul[0].header['exp_times'])
        field_object.append(hdul[0].header['object'])
        field_ntarget.append(len(targets))
        field_nsky.append(len(skies))
        field_subframes.append(hdul[0].header['sub_frames'])

    field_center=np.array(field_center)
    field_id=np.array(field_id)
    field_utdate=np.array(field_utdate)
    field_uttime=np.array(field_uttime)
    field_exptime=np.array(field_exptime)
    field_object=np.array(field_object)
    field_ntarget=np.array(field_ntarget)
    field_nsky=np.array(field_nsky)
    field_subframes=np.array(field_subframes)

    order=np.lexsort((field_uttime,field_utdate))
    field_center=field_center[order]
    field_id=field_id[order]
    field_utdate=field_utdate[order]
    field_uttime=field_uttime[order]
    field_exptime=field_exptime[order]
    field_object=field_object[order]
    field_ntarget=field_ntarget[order]
    field_nsky=field_nsky[order]
    field_subframes=field_subframes[order]
    
    field_object_formal=np.empty(len(field_object),dtype='object')
    dsph=np.where(field_object=='ant2')[0]
    field_object_formal[dsph]='Antlia II'
    dsph=np.where(field_object=='car')[0]
    field_object_formal[dsph]='Carina'
    dsph=np.where(field_object=='cra')[0]
    field_object_formal[dsph]='Crater'
    dsph=np.where(field_object=='cra2')[0]
    field_object_formal[dsph]='Crater II'
    dsph=np.where(field_object=='for')[0]
    field_object_formal[dsph]='Fornax'
    dsph=np.where(field_object=='gru1')[0]
    field_object_formal[dsph]='Grus I'
    dsph=np.where(field_object=='gru2')[0]
    field_object_formal[dsph]='Grus II'
    dsph=np.where(field_object=='hor2')[0]
    field_object_formal[dsph]='Horologium II'
    dsph=np.where(field_object=='hyd1')[0]
    field_object_formal[dsph]='Hyd 1'
    dsph=np.where(field_object=='ind1')[0]
    field_object_formal[dsph]='Indus I'
    dsph=np.where(field_object=='ind2')[0]
    field_object_formal[dsph]='Indus II'
    dsph=np.where(field_object=='kgo2')[0]
    field_object_formal[dsph]='KGO 2'
    dsph=np.where(field_object=='kgo4')[0]
    field_object_formal[dsph]='KGO 4'
    dsph=np.where(field_object=='kgo7')[0]
    field_object_formal[dsph]='KGO 7'
    dsph=np.where(field_object=='kgo8')[0]
    field_object_formal[dsph]='KGO 8'
    dsph=np.where(field_object=='kgo10')[0]
    field_object_formal[dsph]='KGO 10'
    dsph=np.where(field_object=='kgo13')[0]
    field_object_formal[dsph]='KGO 13'
    dsph=np.where(field_object=='kgo22')[0]
    field_object_formal[dsph]='KGO 22'
    dsph=np.where(field_object=='kop2')[0]
    field_object_formal[dsph]='Koposov 2'    
    dsph=np.where(field_object=='pal5')[0]
    field_object_formal[dsph]='Palomar 5'
    dsph=np.where(field_object=='pho2')[0]
    field_object_formal[dsph]='Phoenix II'
    dsph=np.where(field_object=='ret2')[0]
    field_object_formal[dsph]='Reticulum II'
    dsph=np.where(field_object=='sgr2')[0]
    field_object_formal[dsph]='Sagittarius 2'
    dsph=np.where(field_object=='scl')[0]
    field_object_formal[dsph]='Sculptor'
    dsph=np.where(field_object=='sex')[0]
    field_object_formal[dsph]='Sextans'
    dsph=np.where(field_object=='tuc2')[0]
    field_object_formal[dsph]='Tucana II'
    dsph=np.where(field_object=='tuc3')[0]
    field_object_formal[dsph]='Tucana III'
    dsph=np.where(field_object=='tuc4')[0]
    field_object_formal[dsph]='Tucana IV'
    dsph=np.where(field_object=='tuc5')[0]
    field_object_formal[dsph]='Tucana V'

    out_table=open('m2fs_obs_table.tex','w')
    out_table_short=open('m2fs_obs_table_short.tex','w')
    used=np.zeros(len(field_center),dtype='long')
    fudge_uttime=SkyCoord(field_uttime,np.full(len(field_uttime),'+00:00:00'),unit=(u.hourangle,u.deg)).ra.hour

    out_table.write(r'\begin{deluxetable*}{llllcccl}'+' \n')
    out_table.write(r'\tablewidth{0pt}'+' \n')
    out_table.write(r'\tablecaption{Log of M2FS Observations of Galactic Halo Objects}'+' \label{tab:m2fs_obs_table} \n')
    out_table.write(r'\tablehead{\multicolumn{2}{l}{\colhead{Field Center}}&\colhead{UT date\tablenotemark{a}}&\colhead{UT start\tablenotemark{b}}&\colhead{Exp. Time}&\colhead{$N_{\rm exp}$}&\colhead{$N_{\rm target}$}&\colhead{Object}\\'+' \n')
    out_table.write(r'\colhead{$\alpha_{2000}$ [deg.]}&\colhead{$\delta_{2000}$ [deg.]}&&&\colhead{[sec.]}&}'+' \n')
    out_table.write(r'\startdata'+' \n')

    out_table_short.write(r'\begin{deluxetable*}{llllcccl}'+' \n')
    out_table_short.write(r'\tablewidth{0pt}'+' \n')
    out_table_short.write(r'\tablecaption{Log of M2FS Observations of Galactic Halo Objects (abbreviated---see electronic version for full table)}'+' \label{tab:m2fs_obs_table} \n')
    out_table_short.write(r'\tablehead{\multicolumn{2}{l}{\colhead{Field Center}}&\colhead{UT date\tablenotemark{a}}&\colhead{UT start\tablenotemark{b}}&\colhead{Exp. Time}&\colhead{$N_{\rm exp}$}&\colhead{$N_{\rm target}$}&\colhead{Object}\\'+' \n')
    out_table_short.write(r'\colhead{$\alpha_{2000}$ [deg.]}&\colhead{$\delta_{2000}$ [deg.]}&&&\colhead{[sec.]}&}'+' \n')
    out_table_short.write(r'\startdata'+' \n')

    for i in range(0,len(order)):
        if used[i]==0:
            this=np.where((field_id==field_id[i])&(field_object==field_object[i])&(field_utdate==field_utdate[i])&(np.abs(fudge_uttime-fudge_uttime[i])<0.5))[0]
            used[this]=1
            if field_center[i].ra.deg<10.:
                rastring='$00'+f'{field_center[i].ra.deg:.6f}'
            elif field_center[i].ra.deg<100.:
                rastring='$0'+f'{field_center[i].ra.deg:.6f}'
            else:
                rastring='$'+f'{field_center[i].ra.deg:.6f}'
            if field_center[i].dec.deg<0:
                if np.abs(field_center[i].dec.deg)<10.:
                    decstring='$-00$'+f'{np.abs(field_center[i].dec.deg):.6f}'
                elif np.abs(field_center[i].dec.deg)<100.:
                    decstring='$-0$'+f'{np.abs(field_center[i].dec.deg):.6f}'
                else:
                    decstring='$-$'+f'{np.abs(field_center[i].dec.deg):.6f}'
            else:
                if np.abs(field_center[i].dec.deg)<10.:
                    decstring='$+00$'+f'{np.abs(field_center[i].dec.deg):.6f}'
                elif np.abs(field_center[i].dec.deg)<100.:
                    decstring='$+0$'+f'{np.abs(field_center[i].dec.deg):.6f}'
                else:
                    decstring='$+$'+f'{np.abs(field_center[i].dec.deg):.6f}'
                
            string=rastring+'$ & $'+decstring+'$ & '+field_utdate[i]+' & '+field_uttime[i]+' &$'+str(int(np.sum(np.array(field_exptime[i].split(','),dtype=float))))+'$ & $'+str(len(np.array(field_exptime[i].split(','),dtype=float)))+'$ & $'+str(int(np.sum(field_ntarget[this])))+'$ & '+str(field_object_formal[i])+'\\\ \n'
            out_table.write(string)
            if i<25:
                out_table_short.write(string)
            print(this,field_subframes[this])
    out_table.write(r'\enddata'+' \n')
    out_table.write(r'\tablenotetext{a}{YYYY-MM-DD format}'+' \n')
    out_table.write(r'\tablenotetext{b}{Universal time at start of first exposure; HH:MM:SS format}'+' \n')
    out_table.write(r'\end{deluxetable*}'+' \n')

    out_table_short.write(r'\enddata'+' \n')
    out_table_short.write(r'\tablenotetext{a}{YYYY-MM-DD format}'+' \n')
    out_table_short.write(r'\tablenotetext{b}{Universal time at start of first exposure; HH:MM:SS format}'+' \n')
    out_table_short.write(r'\end{deluxetable*}'+' \n')

    out_table.close()
    out_table_short.close()
    pickle.dump(field_center,open('m2fs_field_center.pkl','wb'))

if make_skysub:

    skysub_list='/hildafs/projects/phy200028p/mgwalker/scripts/'+fits_list0.split('_files')[0]
    g00=open(skysub_list,'w')
    g11=open(data_filename,'w')
    g12=open(chains_filename,'w')

    for i in range(0,len(fitsfile)):
        print(i,' of ',len(fitsfile))
        hdul=fits.open(fitsfile[i])
        fitsobject=m2fs.m2fs_getfromfits(hdul)
        filtername=hdul[0].header['filtername']
        m2fsrun=hdul[0].header['m2fsrun']
        field_name=hdul[0].header['field_name']
        temperature=hdul[0].header['CCD_temp']
        hdul.close()
        surveyphot_filename=fitsfile[i].split('_stackskysub.fits')[0]+'_surveyphot.pkl'
    
        root=[]
        root2=[]
        for j in range(0,len(fitsobject.obj)):
            root.append('m2fs_ian_'+m2fsrun+'_'+obj[i]+'_'+field_name+'_'+filtername+'_ra'+str.format('{0:.6f}',round(fitsobject.radeg[j],6)).zfill(6)+'_dec'+str.format('{0:.6f}',round(fitsobject.decdeg[j],6)).zfill(6)+'_hjd'+str.format('{0:.3f}',round(fitsobject.hjd[j],3))+'_ap'+fitsobject.channel[j]+str(int(fitsobject.aperture[j])).zfill(3))
            root2.append('m2fs_ian_'+m2fsrun+'_'+obj[i]+'_'+field_name+'_'+filtername+'_ra'+str.format('{0:.6f}',round(fitsobject.radeg[j],6)).zfill(6)+'_dec'+str.format('{0:.6f}',round(fitsobject.decdeg[j],6)).zfill(6)+'_hjd'+str.format('{0:.3f}',round(fitsobject.hjd[j],3))+'_ap'+fitsobject.channel[j]+str(int(fitsobject.aperture[j])).zfill(3))
        root=np.array(root)
        root2=np.array(root2)

        skies=np.where(fitsobject.obj=='SKY')[0]
        targets=np.where(fitsobject.icode>0)[0]

        if len(targets)==0:
            print('WARNING: NO TARGETS')
        for j in targets:
            out=data_directory+root[j]+'_skysub.dat'
            if len(np.where(fitsobject.mask[j]==False)[0])>0:#write skysub.dat file only if >100 good pixels in spectrum
                print('writing .skysub.dat file for frame ',i,root[j],j)
                g1=open(out,'w')
                g00.write(out+' '+chains_directory+root2[j]+' \n')
                g11.write(out+' \n')
                g12.write(chains_directory+root2[j]+'post_equal_weights.dat'+' \n')
                for k in range(0,len(fitsobject.wav[j])):
#                    if fitsobject.mask[j][k]==False:
                    string=str(round(fitsobject.wav[j][k],10))+' '+str(round(fitsobject.spec[j][k],3))+' '+str(round(fitsobject.var[j][k],5))+' \n'
                    g1.write(string)
                g1.close()
            else:
                print('NOT writing .skysub.dat file for frame ',i,root[j],j)
    g00.close()
    g11.close()
    g12.close()

if get_fits:

    for i in range(0,len(fitsfile)):
        print(i,' of ',len(fitsfile))
        hdul=fits.open(fitsfile[i])
        fitsobject=m2fs.m2fs_getfromfits(hdul)
        filtername=hdul[0].header['filtername']
        m2fsrun=hdul[0].header['m2fsrun']
        field_name=hdul[0].header['field_name']
        temperature=hdul[0].header['CCD_temp']
        hdul.close()
        surveyphot_filename=fitsfile[i].split('_stackskysub.fits')[0]+'_surveyphot.pkl'
    
        root=[]
        root2=[]
        for j in range(0,len(fitsobject.obj)):
            root.append('m2fs_ian_'+m2fsrun+'_'+obj[i]+'_'+field_name+'_'+filtername+'_ra'+str.format('{0:.6f}',round(fitsobject.radeg[j],6)).zfill(6)+'_dec'+str.format('{0:.6f}',round(fitsobject.decdeg[j],6)).zfill(6)+'_hjd'+str.format('{0:.3f}',round(fitsobject.hjd[j],3))+'_ap'+fitsobject.channel[j]+str(int(fitsobject.aperture[j])).zfill(3))
            root2.append('m2fs_ian_'+m2fsrun+'_'+obj[i]+'_'+field_name+'_'+filtername+'_ra'+str.format('{0:.6f}',round(fitsobject.radeg[j],6)).zfill(6)+'_dec'+str.format('{0:.6f}',round(fitsobject.decdeg[j],6)).zfill(6)+'_hjd'+str.format('{0:.3f}',round(fitsobject.hjd[j],3))+'_ap'+fitsobject.channel[j]+str(int(fitsobject.aperture[j])).zfill(3))
        root=np.array(root)
        root2=np.array(root2)

        skies=np.where(fitsobject.obj=='SKY')[0]
        targets=np.where(fitsobject.icode>0)[0]

        survey_phot=pickle.load(open(surveyphot_filename,'rb'))
        gaia_ra=survey_phot.gaia.ra
        gaia_dec=survey_phot.gaia.dec
        gaia_gmag_dered=survey_phot.gaia.gmag-survey_phot.gaia.gext
        gaia_bpmag_dered=survey_phot.gaia.bpmag-survey_phot.gaia.bpext
        gaia_rpmag_dered=survey_phot.gaia.rpmag-survey_phot.gaia.rpext
        des_gmag_dered=survey_phot.des.gmag-survey_phot.des.gext
        des_rmag_dered=survey_phot.des.rmag-survey_phot.des.rext
        des_zmag_dered=survey_phot.des.zmag-survey_phot.des.zext
        decals_gmag_dered=survey_phot.decals.gmag-survey_phot.decals.gext
        decals_rmag_dered=survey_phot.decals.rmag-survey_phot.decals.gext
        decals_zmag_dered=survey_phot.decals.zmag-survey_phot.decals.gext
        sdss_gmag_dered=survey_phot.sdss.gmag-survey_phot.sdss.gext
        sdss_rmag_dered=survey_phot.sdss.rmag-survey_phot.sdss.gext
        sdss_zmag_dered=survey_phot.sdss.zmag-survey_phot.sdss.gext
        ps1_gmag_dered=survey_phot.ps1.gmag-survey_phot.ps1.gext
        ps1_rmag_dered=survey_phot.ps1.rmag-survey_phot.ps1.gext
        ps1_zmag_dered=survey_phot.ps1.zmag-survey_phot.ps1.gext
        
        gaia_teffprior=((gaia_est50,gaia_est16,gaia_est84),gaia_gmag_dered,gaia_bpmag_dered,gaia_rpmag_dered)
        des_teffprior=((des_est50,des_est16,des_est84),des_gmag_dered,des_rmag_dered,des_zmag_dered)
        decals_teffprior=((decals_est50,decals_est16,decals_est84),decals_gmag_dered,decals_rmag_dered,decals_zmag_dered)
        sdss_teffprior=((sdss_est50,sdss_est16,sdss_est84),sdss_gmag_dered,sdss_rmag_dered,sdss_zmag_dered)
        ps1_teffprior=((ps1_est50,ps1_est16,ps1_est84),ps1_gmag_dered,ps1_rmag_dered,ps1_zmag_dered)

        if len(fitsobject.obj)>0:
            teffpriorstuff0=teffpriorstuff(gaia_teffprior,des_teffprior,decals_teffprior,sdss_teffprior,ps1_teffprior)
            multinest=m2fs.m2fs_multinest_ian(fit_directory,root2,targets,fitsobject,len(fitsobject.wav[0]),teffpriorstuff0)#object containing sample of posterior, moments thereof, bestfit wavelength and bestfit spectrum
            hdr=fits.Header(fitsobject.header)
            primary_hdu=fits.PrimaryHDU(header=hdr)
            new_hdul=fits.HDUList([primary_hdu])

            update_coords_to_gaia=np.where((gaia_ra==gaia_ra)&(gaia_dec==gaia_dec))[0]#if found gaia match, update coords to gaia EDR3
            fitsobject.radeg[update_coords_to_gaia]=gaia_ra[update_coords_to_gaia]
            fitsobject.decdeg[update_coords_to_gaia]=gaia_dec[update_coords_to_gaia]
            col1=fits.Column(name='ra_deg',format='D',array=fitsobject.radeg)
            col2=fits.Column(name='dec_deg',format='D',array=fitsobject.decdeg)
            col3=fits.Column(name='mjd',format='D',array=fitsobject.mjd)
            col4=fits.Column(name='hjd',format='D',array=fitsobject.hjd)
            col5=fits.Column(name='snratio',format='D',array=fitsobject.snratio)
            col6=fits.Column(name='vhelio_correction',format='D',array=fitsobject.vheliocorr)
            col7=fits.Column(name='posterior_moments',format='68D',dim='(17,4)',array=multinest.moments)
            col8=fits.Column(name='posterior_1000',format='17000D',dim='(17,1000)',array=multinest.posterior_1000)
            col9=fits.Column(name='teffprior_posterior_moments',format='68D',dim='(17,4)',array=multinest.teffprior_moments)
            col10=fits.Column(name='teffprior_posterior_1000',format='17000D',dim='(17,1000)',array=multinest.teffprior_posterior_1000)
            col11=fits.Column(name='teffprior_survey',format='A10',array=multinest.teffprior_survey)
            col12=fits.Column(name='objtype',format='A6',array=fitsobject.obj)
            col13=fits.Column(name='run_id',format='A100',array=fitsobject.run_id)
            col14=fits.Column(name='field_name',format='A100',array=fitsobject.field_name)
            col15=fits.Column(name='filtername',format='A100',array=fitsobject.filtername)
            col16=fits.Column(name='wav_npoints',format='PI()',array=fitsobject.wav_npoints)
            col17=fits.Column(name='wav_rms',format='PD()',array=fitsobject.wav_rms)
            col18=fits.Column(name='wav_resolution',format='PD()',array=fitsobject.wav_resolution)
            col19=fits.Column(name='wav_min',format='PD()',array=fitsobject.wav_min)
            col20=fits.Column(name='wav_max',format='PD()',array=fitsobject.wav_max)
            col21=fits.Column(name='row',format='D',array=fitsobject.row)
            col22=fits.Column(name='temperature',format='PD()',array=fitsobject.temperature)
            col23=fits.Column(name='gaia_source_id',format='A30',array=survey_phot.gaia.objid)
            col24=fits.Column(name='gaia_gmag',format='D',array=survey_phot.gaia.gmag)
            col25=fits.Column(name='gaia_siggmag',format='D',array=survey_phot.gaia.siggmag)
            col26=fits.Column(name='gaia_gmag_dered',format='D',array=gaia_gmag_dered)
            col27=fits.Column(name='gaia_bpmag',format='D',array=survey_phot.gaia.bpmag)
            col28=fits.Column(name='gaia_sigbpmag',format='D',array=survey_phot.gaia.sigbpmag)
            col29=fits.Column(name='gaia_bpmag_dered',format='D',array=gaia_bpmag_dered)
            col30=fits.Column(name='gaia_rpmag',format='D',array=survey_phot.gaia.rpmag)
            col31=fits.Column(name='gaia_sigrpmag',format='D',array=survey_phot.gaia.sigrpmag)
            col32=fits.Column(name='gaia_rpmag_dered',format='D',array=gaia_rpmag_dered)
            col33=fits.Column(name='gaia_pmra',format='D',array=survey_phot.gaia.pmra)
            col34=fits.Column(name='gaia_sigpmra',format='D',array=survey_phot.gaia.sigpmra)
            col35=fits.Column(name='gaia_pmdec',format='D',array=survey_phot.gaia.pmdec)
            col36=fits.Column(name='gaia_sigpmdec',format='D',array=survey_phot.gaia.sigpmdec)
            col37=fits.Column(name='gaia_parallax',format='D',array=survey_phot.gaia.parallax)
            col38=fits.Column(name='gaia_sigparallax',format='D',array=survey_phot.gaia.sigparallax)
            col39=fits.Column(name='gaia_rrl',format='D',array=survey_phot.gaia.rrl)
            cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34,col35,col36,col37,col38,col39])
            table_hdu=fits.FITS_rec.from_columns(cols)#BinTableHDU.from_columns(cols)

            primary_hdu=fits.PrimaryHDU([],header=hdr)
            new_hdul=fits.HDUList([primary_hdu])
            new_hdul.append(fits.ImageHDU(fitsobject.wav,name='wavlength'))
            new_hdul.append(fits.ImageHDU(fitsobject.spec,name='sky_subtracted'))
            new_hdul.append(fits.ImageHDU(fitsobject.var,name='variance'))
            new_hdul.append(fits.ImageHDU(fitsobject.mask,name='mask'))
            new_hdul.append(fits.ImageHDU(multinest.bestfit_fit,name='bestfit'))
            new_hdul.append(fits.BinTableHDU(table_hdu,name='data_table'))
            new_hdul.append(fits.ImageHDU(fitsobject.sky_spec,name='sky'))

            outfile=fitsfile[i].split('.fits')[0]+'_postfit_ian.fits'
            print('writing fit results to '+outfile)
            new_hdul.writeto(outfile,overwrite=True)

            for k in range(0,len(new_hdul[4].data)):
                rastring=str('{0:.7f}'.format(new_hdul[6].data['ra_deg'][k]))
                decstring=str('{0:.7f}'.format(new_hdul[6].data['dec_deg'][k]))
                hjdstring=str('{0:.3f}'.format(new_hdul[6].data['hjd'][k]))
                if new_hdul[6].data['dec_deg'][k]>=0.:
                    radechjdstring=rastring+'_+'+decstring+'_'+hjdstring
                else:
                    radechjdstring=rastring+'_'+decstring+'_'+hjdstring

if get_catalog_raw:
    
    postfit_object=[]
    postfit_obsid=[]
    postfit_index=[]
    postfit_objtype=[]
    postfit_radeg=[]
    postfit_decdeg=[]
    postfit_mjd=[]
    postfit_hjd=[]
    postfit_snratio=[]
    postfit_gaiaid=[]
    postfit_gaiagmag=[]
    postfit_siggaiagmag=[]
    postfit_gaiagmag_dered=[]
    postfit_gaiabpmag=[]
    postfit_siggaiabpmag=[]
    postfit_gaiabpmag_dered=[]
    postfit_gaiarpmag=[]
    postfit_siggaiarpmag=[]
    postfit_gaiarpmag_dered=[]
    postfit_gaiapmra=[]
    postfit_siggaiapmra=[]
    postfit_gaiapmdec=[]
    postfit_siggaiapmdec=[]
    postfit_gaiaparallax=[]
    postfit_siggaiaparallax=[]
    postfit_gaiarrl=[]
    postfit_vhelio_correction=[]
    postfit_mean=[]
    postfit_std=[]
    postfit_skew=[]
    postfit_kurtosis=[]        
    postfit_teffprior_mean=[]
    postfit_teffprior_std=[]
    postfit_teffprior_skew=[]
    postfit_teffprior_kurtosis=[]
    postfit_teffprior_survey=[]
    postfit_mediansky=[]
    postfit_stdsky=[]
    postfit_run_id=[]
    postfit_field_name=[]
    postfit_filtername=[]
    postfit_chi2=[]
    postfit_bluechi2=[]
    postfit_redchi2=[]
    postfit_chi22=[]
    postfit_n=[]
    postfit_bluen=[]
    postfit_redn=[]
    postfit_wav_npoints=[]
    postfit_wav_rms=[]
    postfit_wav_resolution=[]
    postfit_wav_min=[]
    postfit_wav_max=[]
    postfit_row=[]
    postfit_temperature=[]
    postfit_filename=[]
    postfit_allmasked=[]
    
    for i in range(0,len(fitsfile)):
        print(i,' of ',len(fitsfile))
        hdul=fits.open(fitsfile[i])
        fitsobject=m2fs.m2fs_getfromfits(hdul)
        filtername=hdul[0].header['filtername']
        m2fsrun=hdul[0].header['m2fsrun']
        field_name=hdul[0].header['field_name']
        temperature=hdul[0].header['CCD_temp']
        hdul.close()
        surveyphot_filename=fitsfile[i].split('_stackskysub.fits')[0]+'_surveyphot.pkl'
    
        root=[]
        root2=[]
        for j in range(0,len(fitsobject.obj)):
            root.append('m2fs_ian_'+m2fsrun+'_'+obj[i]+'_'+field_name+'_'+filtername+'_ra'+str.format('{0:.6f}',round(fitsobject.radeg[j],6)).zfill(6)+'_dec'+str.format('{0:.6f}',round(fitsobject.decdeg[j],6)).zfill(6)+'_hjd'+str.format('{0:.3f}',round(fitsobject.hjd[j],3))+'_ap'+fitsobject.channel[j]+str(int(fitsobject.aperture[j])).zfill(3))
            root2.append('m2fs_ian_'+m2fsrun+'_'+obj[i]+'_'+field_name+'_'+filtername+'_ra'+str.format('{0:.6f}',round(fitsobject.radeg[j],6)).zfill(6)+'_dec'+str.format('{0:.6f}',round(fitsobject.decdeg[j],6)).zfill(6)+'_hjd'+str.format('{0:.3f}',round(fitsobject.hjd[j],3))+'_ap'+fitsobject.channel[j]+str(int(fitsobject.aperture[j])).zfill(3))
        root=np.array(root)
        root2=np.array(root2)

        all_masked=np.zeros(len(fitsobject.mask),dtype='int')
        for j in range(0,len(fitsobject.mask)):
            if np.min(fitsobject.mask[j])>0:
                all_masked[j]=1
                
        skies=np.where((fitsobject.obj=='SKY')&(all_masked==0))[0]
        targets=np.where((fitsobject.icode>0)&(all_masked==0))[0]
    
        if len(fitsobject.obj)>0:
            infile=fitsfile[i].split('.fits')[0]+'_postfit_ian.fits'
            postfit=fits.open(infile)
            mediansky=[]
            for j in range(0,len(postfit[6].data)):
                rastring=str('{0:.7f}'.format(postfit[6].data['ra_deg'][j]))
                decstring=str('{0:.7f}'.format(postfit[6].data['dec_deg'][j]))
                hjdstring=str('{0:.3f}'.format(postfit[6].data['hjd'][j]))
                if postfit[6].data['dec_deg'][j]>=0.:
                    radechjdstring=rastring+'_+'+decstring+'_'+hjdstring
                else:
                    radechjdstring=rastring+'_'+decstring+'_'+hjdstring
                postfit_temperature.append(postfit[6].data['temperature'][j])
                postfit_obsid.append(radechjdstring)
                postfit_index.append(j)
                postfit_objtype.append(postfit[6].data['objtype'][j])
                postfit_radeg.append(postfit[6].data['ra_deg'][j])
                postfit_decdeg.append(postfit[6].data['dec_deg'][j])
                postfit_mjd.append(postfit[6].data['mjd'][j])
                postfit_hjd.append(postfit[6].data['hjd'][j])
                postfit_snratio.append(postfit[6].data['snratio'][j])
                postfit_gaiaid.append(postfit[6].data['gaia_source_id'][j])
                postfit_gaiagmag.append(postfit[6].data['gaia_gmag'][j])
                postfit_siggaiagmag.append(postfit[6].data['gaia_siggmag'][j])
                postfit_gaiagmag_dered.append(postfit[6].data['gaia_gmag_dered'][j])
                postfit_gaiabpmag.append(postfit[6].data['gaia_bpmag'][j])
                postfit_siggaiabpmag.append(postfit[6].data['gaia_sigbpmag'][j])
                postfit_gaiabpmag_dered.append(postfit[6].data['gaia_bpmag_dered'][j])
                postfit_gaiarpmag.append(postfit[6].data['gaia_rpmag'][j])
                postfit_siggaiarpmag.append(postfit[6].data['gaia_sigrpmag'][j])
                postfit_gaiarpmag_dered.append(postfit[6].data['gaia_rpmag_dered'][j])
                postfit_gaiapmra.append(postfit[6].data['gaia_pmra'][j])
                postfit_siggaiapmra.append(postfit[6].data['gaia_sigpmra'][j])
                postfit_gaiapmdec.append(postfit[6].data['gaia_pmdec'][j])
                postfit_siggaiapmdec.append(postfit[6].data['gaia_sigpmdec'][j])
                postfit_gaiaparallax.append(postfit[6].data['gaia_parallax'][j])
                postfit_siggaiaparallax.append(postfit[6].data['gaia_sigparallax'][j])
                postfit_gaiarrl.append(postfit[6].data['gaia_rrl'][j])
                postfit_vhelio_correction.append(postfit[6].data['vhelio_correction'][j])
                postfit_mean.append(postfit[6].data['posterior_moments'][j][0])
                postfit_std.append(postfit[6].data['posterior_moments'][j][1])
                postfit_skew.append(postfit[6].data['posterior_moments'][j][2])
                postfit_kurtosis.append(postfit[6].data['posterior_moments'][j][3])
                postfit_teffprior_mean.append(postfit[6].data['teffprior_posterior_moments'][j][0])
                postfit_teffprior_std.append(postfit[6].data['teffprior_posterior_moments'][j][1])
                postfit_teffprior_skew.append(postfit[6].data['teffprior_posterior_moments'][j][2])
                postfit_teffprior_kurtosis.append(postfit[6].data['teffprior_posterior_moments'][j][3])
                postfit_teffprior_survey.append(postfit[6].data['teffprior_survey'][j])
                keep=np.where((postfit[4].data[j]==0)&(postfit[3].data[j]<1.e+5))[0]
                bluekeep=np.where((postfit[4].data[j]==0)&(postfit[3].data[j]<1.e+5)&(postfit[1].data[j]>lambdamin)&(postfit[1].data[j]<5165.))[0]
                redkeep=np.where((postfit[4].data[j]==0)&(postfit[3].data[j]<1.e+5)&(postfit[1].data[j]<lambdamax)&(postfit[1].data[j]>=5165.))[0]
                mediansky.append(np.median(postfit[7].data[j][keep]))
                postfit_mediansky.append(np.median(postfit[7].data[j][keep]))
                postfit_object.append(obj[i])
                postfit_run_id.append(postfit[6].data['run_id'][j])
                postfit_field_name.append(postfit[6].data['field_name'][j])
                postfit_filtername.append(postfit[6].data['filtername'][j])
                postfit_wav_npoints.append(postfit[6].data['wav_npoints'][j])
                postfit_wav_rms.append(postfit[6].data['wav_rms'][j])
                postfit_wav_resolution.append(postfit[6].data['wav_resolution'][j])
                postfit_wav_min.append(postfit[6].data['wav_min'][j])
                postfit_wav_max.append(postfit[6].data['wav_max'][j])
                postfit_row.append(postfit[6].data['row'][j])
                postfit_temperature.append(postfit[6].data['temperature'][j])
                postfit_filename.append(infile)
                postfit_allmasked.append(all_masked[j])
                
                spec=postfit[2].data[j]
                mask=postfit[4].data[j]
                bestfit=postfit[5].data[j]
                var=postfit[3].data[j]

                chi2=np.sum((spec[keep]-bestfit[keep])**2/var[keep])
                bluechi2=np.median((spec[bluekeep]-bestfit[bluekeep])**2/var[bluekeep])
                redchi2=np.median((spec[redkeep]-bestfit[redkeep])**2/var[redkeep])
                
                postfit_chi2.append(chi2)
                postfit_n.append(len(keep))
                postfit_bluechi2.append(bluechi2)
                postfit_bluen.append(len(bluekeep))
                postfit_redchi2.append(redchi2)
                postfit_redn.append(len(redkeep))

                sigma0=postfit[6].data['posterior_moments'][j][0][14]
                sigma1=postfit[6].data['posterior_moments'][j][0][15]
                var2=10.**sigma1*var+(10.**sigma0)**2
                chi22=np.sum((spec[keep]-bestfit[keep])**2/var2[keep])
                postfit_chi22.append(chi22)
            mediansky=np.array(mediansky)
            stdsky=np.std(mediansky[skies])
            for j in range(0,len(postfit[6].data)):
                postfit_stdsky.append(stdsky)

    postfit_allmasked=np.array(postfit_allmasked)
    postfit_object=np.array(postfit_object)[np.where(postfit_allmasked==0)[0]]
    postfit_obsid=np.array(postfit_obsid)[np.where(postfit_allmasked==0)[0]]
    postfit_index=np.array(postfit_index)[np.where(postfit_allmasked==0)[0]]
    postfit_objtype=np.array(postfit_objtype)[np.where(postfit_allmasked==0)[0]]
    postfit_radeg=np.array(postfit_radeg)[np.where(postfit_allmasked==0)[0]]
    postfit_decdeg=np.array(postfit_decdeg)[np.where(postfit_allmasked==0)[0]]
    postfit_mjd=np.array(postfit_mjd)[np.where(postfit_allmasked==0)[0]]
    postfit_hjd=np.array(postfit_hjd)[np.where(postfit_allmasked==0)[0]]
    postfit_snratio=np.array(postfit_snratio)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiaid=np.array(postfit_gaiaid)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiagmag=np.array(postfit_gaiagmag)[np.where(postfit_allmasked==0)[0]]
    postfit_siggaiagmag=np.array(postfit_siggaiagmag)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiagmag_dered=np.array(postfit_gaiagmag_dered)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiabpmag=np.array(postfit_gaiabpmag)[np.where(postfit_allmasked==0)[0]]
    postfit_siggaiabpmag=np.array(postfit_siggaiabpmag)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiabpmag_dered=np.array(postfit_gaiabpmag_dered)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiarpmag=np.array(postfit_gaiarpmag)[np.where(postfit_allmasked==0)[0]]
    postfit_siggaiarpmag=np.array(postfit_siggaiarpmag)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiarpmag_dered=np.array(postfit_gaiarpmag_dered)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiapmra=np.array(postfit_gaiapmra)[np.where(postfit_allmasked==0)[0]]
    postfit_siggaiapmra=np.array(postfit_siggaiapmra)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiapmdec=np.array(postfit_gaiapmdec)[np.where(postfit_allmasked==0)[0]]
    postfit_siggaiapmdec=np.array(postfit_siggaiapmdec)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiaparallax=np.array(postfit_gaiaparallax)[np.where(postfit_allmasked==0)[0]]
    postfit_siggaiaparallax=np.array(postfit_siggaiaparallax)[np.where(postfit_allmasked==0)[0]]
    postfit_gaiarrl=np.array(postfit_gaiarrl)[np.where(postfit_allmasked==0)[0]]
    postfit_vhelio_correction=np.array(postfit_vhelio_correction)[np.where(postfit_allmasked==0)[0]]
    postfit_mean=np.array(postfit_mean)[np.where(postfit_allmasked==0)[0]]
    postfit_std=np.array(postfit_std)[np.where(postfit_allmasked==0)[0]]
    postfit_skew=np.array(postfit_skew)[np.where(postfit_allmasked==0)[0]]
    postfit_kurtosis=np.array(postfit_kurtosis)[np.where(postfit_allmasked==0)[0]]
    postfit_teffprior_mean=np.array(postfit_teffprior_mean)[np.where(postfit_allmasked==0)[0]]
    postfit_teffprior_std=np.array(postfit_teffprior_std)[np.where(postfit_allmasked==0)[0]]
    postfit_teffprior_skew=np.array(postfit_teffprior_skew)[np.where(postfit_allmasked==0)[0]]
    postfit_teffprior_kurtosis=np.array(postfit_teffprior_kurtosis)[np.where(postfit_allmasked==0)[0]]
    postfit_teffprior_survey=np.array(postfit_teffprior_survey)[np.where(postfit_allmasked==0)[0]]
    postfit_mediansky=np.array(postfit_mediansky)[np.where(postfit_allmasked==0)[0]]
    postfit_stdsky=np.array(postfit_stdsky)[np.where(postfit_allmasked==0)[0]]
    postfit_run_id=np.array(postfit_run_id)[np.where(postfit_allmasked==0)[0]]
    postfit_field_name=np.array(postfit_field_name)[np.where(postfit_allmasked==0)[0]]
    postfit_filtername=np.array(postfit_filtername)[np.where(postfit_allmasked==0)[0]]
    postfit_chi2=np.array(postfit_chi2)[np.where(postfit_allmasked==0)[0]]
    postfit_bluechi2=np.array(postfit_bluechi2)[np.where(postfit_allmasked==0)[0]]
    postfit_redchi2=np.array(postfit_redchi2)[np.where(postfit_allmasked==0)[0]]
    postfit_chi22=np.array(postfit_chi22)[np.where(postfit_allmasked==0)[0]]
    postfit_n=np.array(postfit_n)[np.where(postfit_allmasked==0)[0]]
    postfit_bluen=np.array(postfit_bluen)[np.where(postfit_allmasked==0)[0]]
    postfit_redn=np.array(postfit_redn)[np.where(postfit_allmasked==0)[0]]
    postfit_wav_npoints=np.array(postfit_wav_npoints)[np.where(postfit_allmasked==0)[0]]
    postfit_wav_rms=np.array(postfit_wav_rms)[np.where(postfit_allmasked==0)[0]]
    postfit_wav_resolution=np.array(postfit_wav_resolution)[np.where(postfit_allmasked==0)[0]]
    postfit_wav_min=np.array(postfit_wav_min)[np.where(postfit_allmasked==0)[0]]
    postfit_wav_max=np.array(postfit_wav_max)[np.where(postfit_allmasked==0)[0]]
    postfit_row=np.array(postfit_row)[np.where(postfit_allmasked==0)[0]]
    postfit_temperature=np.array(postfit_temperature)[np.where(postfit_allmasked==0)[0]]
    postfit_filename=np.array(postfit_filename)[np.where(postfit_allmasked==0)[0]]

    v=np.array([postfit_mean[q][0] for q in range(0,len(postfit_mean))])
    v_sig=np.array([postfit_std[q][0] for q in range(0,len(postfit_std))])
    v_skew=np.array([postfit_skew[q][0] for q in range(0,len(postfit_skew))])
    v_kurt=np.array([postfit_kurtosis[q][0] for q in range(0,len(postfit_kurtosis))])
    teff=np.array([postfit_mean[q][1] for q in range(0,len(postfit_mean))])
    teff_sig=np.array([postfit_std[q][1] for q in range(0,len(postfit_std))])
    teff_skew=np.array([postfit_skew[q][1] for q in range(0,len(postfit_skew))])
    teff_kurt=np.array([postfit_kurtosis[q][1] for q in range(0,len(postfit_kurtosis))])
    logg=np.array([postfit_mean[q][2] for q in range(0,len(postfit_mean))])
    logg_sig=np.array([postfit_std[q][2] for q in range(0,len(postfit_std))])
    logg_skew=np.array([postfit_skew[q][2] for q in range(0,len(postfit_skew))])
    logg_kurt=np.array([postfit_kurtosis[q][2] for q in range(0,len(postfit_kurtosis))])
    z=np.array([postfit_mean[q][3] for q in range(0,len(postfit_mean))])
    z_sig=np.array([postfit_std[q][3] for q in range(0,len(postfit_std))])
    z_skew=np.array([postfit_skew[q][3] for q in range(0,len(postfit_skew))])
    z_kurt=np.array([postfit_kurtosis[q][3] for q in range(0,len(postfit_kurtosis))])
    mgfe=np.array([postfit_mean[q][4] for q in range(0,len(postfit_mean))])
    mgfe_sig=np.array([postfit_std[q][4] for q in range(0,len(postfit_std))])
    mgfe_skew=np.array([postfit_skew[q][4] for q in range(0,len(postfit_skew))])
    mgfe_kurt=np.array([postfit_kurtosis[q][4] for q in range(0,len(postfit_kurtosis))])
    smooth=np.array([postfit_mean[q][11] for q in range(0,len(postfit_mean))])
    smooth_sig=np.array([postfit_std[q][11] for q in range(0,len(postfit_std))])
    smooth_skew=np.array([postfit_skew[q][11] for q in range(0,len(postfit_skew))])
    smooth_kurt=np.array([postfit_kurtosis[q][11] for q in range(0,len(postfit_kurtosis))])
    phantom1=np.array([postfit_mean[q][14] for q in range(0,len(postfit_mean))])
    phantom1_sig=np.array([postfit_std[q][14] for q in range(0,len(postfit_std))])
    phantom1_skew=np.array([postfit_skew[q][14] for q in range(0,len(postfit_skew))])
    phantom1_kurt=np.array([postfit_kurtosis[q][14] for q in range(0,len(postfit_kurtosis))])
    phantom2=np.array([postfit_mean[q][15] for q in range(0,len(postfit_mean))])
    phantom2_sig=np.array([postfit_std[q][15] for q in range(0,len(postfit_std))])
    phantom2_skew=np.array([postfit_skew[q][15] for q in range(0,len(postfit_skew))])
    phantom2_kurt=np.array([postfit_kurtosis[q][15] for q in range(0,len(postfit_kurtosis))])
    
    teffprior_v=np.array([postfit_teffprior_mean[q][0] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_v_sig=np.array([postfit_teffprior_std[q][0] for q in range(0,len(postfit_teffprior_std))])
    teffprior_v_skew=np.array([postfit_teffprior_skew[q][0] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_v_kurt=np.array([postfit_teffprior_kurtosis[q][0] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_teff=np.array([postfit_teffprior_mean[q][1] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_teff_sig=np.array([postfit_teffprior_std[q][1] for q in range(0,len(postfit_teffprior_std))])
    teffprior_teff_skew=np.array([postfit_teffprior_skew[q][1] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_teff_kurt=np.array([postfit_teffprior_kurtosis[q][1] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_logg=np.array([postfit_teffprior_mean[q][2] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_logg_sig=np.array([postfit_teffprior_std[q][2] for q in range(0,len(postfit_teffprior_std))])
    teffprior_logg_skew=np.array([postfit_teffprior_skew[q][2] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_logg_kurt=np.array([postfit_teffprior_kurtosis[q][2] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_z=np.array([postfit_teffprior_mean[q][3] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_z_sig=np.array([postfit_teffprior_std[q][3] for q in range(0,len(postfit_teffprior_std))])
    teffprior_z_skew=np.array([postfit_teffprior_skew[q][3] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_z_kurt=np.array([postfit_teffprior_kurtosis[q][3] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_mgfe=np.array([postfit_teffprior_mean[q][4] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_mgfe_sig=np.array([postfit_teffprior_std[q][4] for q in range(0,len(postfit_teffprior_std))])
    teffprior_mgfe_skew=np.array([postfit_teffprior_skew[q][4] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_mgfe_kurt=np.array([postfit_teffprior_kurtosis[q][4] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_smooth=np.array([postfit_teffprior_mean[q][11] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_smooth_sig=np.array([postfit_teffprior_std[q][11] for q in range(0,len(postfit_teffprior_std))])
    teffprior_smooth_skew=np.array([postfit_teffprior_skew[q][11] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_smooth_kurt=np.array([postfit_teffprior_kurtosis[q][11] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_phantom1=np.array([postfit_teffprior_mean[q][14] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_phantom1_sig=np.array([postfit_teffprior_std[q][14] for q in range(0,len(postfit_teffprior_std))])
    teffprior_phantom1_skew=np.array([postfit_teffprior_skew[q][14] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_phantom1_kurt=np.array([postfit_teffprior_kurtosis[q][14] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_phantom2=np.array([postfit_teffprior_mean[q][15] for q in range(0,len(postfit_teffprior_mean))])
    teffprior_phantom2_sig=np.array([postfit_teffprior_std[q][15] for q in range(0,len(postfit_teffprior_std))])
    teffprior_phantom2_skew=np.array([postfit_teffprior_skew[q][15] for q in range(0,len(postfit_teffprior_skew))])
    teffprior_phantom2_kurt=np.array([postfit_teffprior_kurtosis[q][15] for q in range(0,len(postfit_teffprior_kurtosis))])
    teffprior_survey=np.array([postfit_teffprior_survey[q] for q in range(0,len(postfit_teffprior_survey))])
    
    v=v+postfit_vhelio_correction#apply heliocentrcic correction!!!
    teffprior_v=teffprior_v+postfit_vhelio_correction#apply heliocentric correction!!!
    
    g1=open(data_out,'w')
    g1.write('#Object # RA [deg] # Dec [deg] # HJD [days] # median S/N/pix # vhelio [km/s] # err_vhelio [km/s] # skew_vhelio # kurt_vhelio # Teff [K] # err_Teff [K] # skew_Teff # kurt_Teff # logg # err_logg # skew_logg # kurt_logg # Fe/H # err_Fe/H # skew_Fe/H # kurt_Fe/H # [Mg/Fe] # err_[Mg/Fe] # skew_[Mg/Fe] # kurt_[Mg/Fe] # resolution \n ')
    for i in range(0,len(postfit_radeg)):
        if postfit_objtype[i]=='TARGET':
            string=postfit_object[i]+' '+str.format('{0:.7f}',round(postfit_radeg[i],7)).zfill(7)+' '+str.format('{0:.7f}',round(postfit_decdeg[i],7)).zfill(7)+' '+str.format('{0:.3f}',round(postfit_hjd[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(postfit_snratio[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(v[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(v_sig[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(v_skew[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(v_kurt[i],3)).zfill(3)+' '+str.format('{0:.1f}',round(teff[i],1)).zfill(1)+' '+str.format('{0:.1f}',round(teff_sig[i],1)).zfill(1)+' '+str.format('{0:.1f}',round(teff_skew[i],1)).zfill(1)+' '+str.format('{0:.1f}',round(teff_kurt[i],1)).zfill(1)+' '+str.format('{0:.3f}',round(logg[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(logg_sig[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(logg_skew[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(logg_kurt[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(z[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(z_sig[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(z_skew[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(z_kurt[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(mgfe[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(mgfe_sig[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(mgfe_skew[i],3)).zfill(3)+' '+str.format('{0:.3f}',round(mgfe_kurt[i],3)).zfill(3)+' '+postfit_filtername[i]+' \n'
            print(string)
            if ((postfit_snratio[i]==postfit_snratio[i])&(not((v_kurt[i]==0.)&(v_kurt[i]==0.)&(z_kurt[i]==0.)&(z_skew[i]==0.)&(logg_kurt[i]==0.)&(logg_skew[i]==0.)))):
                g1.write(string)
    g1.close()

    targets=np.where((postfit_objtype=='TARGET')&(v_sig>1.e-10)&(teff_sig>1.e-10)&(logg_sig>1.e-10)&(z_sig>1.e-10)&(postfit_filtername=='HiRes'))[0]
    col1=fits.Column(name='target_system',format='A100',array=postfit_object[targets])
    col2=fits.Column(name='obj_id',format='A100',array=postfit_obsid[targets])
    col3=fits.Column(name='ra_deg',format='D',array=postfit_radeg[targets])
    col4=fits.Column(name='dec_deg',format='D',array=postfit_decdeg[targets])
    col5=fits.Column(name='mjd',format='D',array=postfit_mjd[targets])
    col6=fits.Column(name='hjd',format='D',array=postfit_hjd[targets])
    col7=fits.Column(name='sn_ratio',format='D',array=postfit_snratio[targets])
    col8=fits.Column(name='vlos',format='D',array=v[targets])
    col9=fits.Column(name='vlos_error_raw',format='D',array=v_sig[targets])
    col10=fits.Column(name='vlos_error',format='D',array=v_sig[targets])
    col11=fits.Column(name='vlos_skew',format='D',array=v_skew[targets])
    col12=fits.Column(name='vlos_kurtosis',format='D',array=v_kurt[targets])
    col13=fits.Column(name='teff',format='D',array=teff[targets])
    col14=fits.Column(name='teff_error_raw',format='D',array=teff_sig[targets])
    col15=fits.Column(name='teff_error',format='D',array=teff_sig[targets])
    col16=fits.Column(name='teff_skew',format='D',array=teff_skew[targets])
    col17=fits.Column(name='teff_kurtosis',format='D',array=teff_kurt[targets])
    col18=fits.Column(name='logg',format='D',array=logg[targets])
    col19=fits.Column(name='logg_error_raw',format='D',array=logg_sig[targets])
    col20=fits.Column(name='logg_error',format='D',array=logg_sig[targets])
    col21=fits.Column(name='logg_skew',format='D',array=logg_skew[targets])
    col22=fits.Column(name='logg_kurtosis',format='D',array=logg_kurt[targets])
    col23=fits.Column(name='feh',format='D',array=z[targets])
    col24=fits.Column(name='feh_error_raw',format='D',array=z_sig[targets])
    col25=fits.Column(name='feh_error',format='D',array=z_sig[targets])
    col26=fits.Column(name='feh_skew',format='D',array=z_skew[targets])
    col27=fits.Column(name='feh_kurtosis',format='D',array=z_kurt[targets])
    col28=fits.Column(name='mgfe',format='D',array=mgfe[targets])
    col29=fits.Column(name='mgfe_error_raw',format='D',array=mgfe_sig[targets])
    col30=fits.Column(name='mgfe_error',format='D',array=mgfe_sig[targets])
    col31=fits.Column(name='mgfe_skew',format='D',array=mgfe_skew[targets])
    col32=fits.Column(name='mgfe_kurtosis',format='D',array=mgfe_kurt[targets])
    col33=fits.Column(name='smooth',format='D',array=smooth[targets])
    col34=fits.Column(name='smooth_error_raw',format='D',array=smooth_sig[targets])
    col35=fits.Column(name='smooth_error',format='D',array=smooth_sig[targets])
    col36=fits.Column(name='smooth_skew',format='D',array=smooth_skew[targets])
    col37=fits.Column(name='smooth_kurtosis',format='D',array=smooth_kurt[targets])
    col38=fits.Column(name='teffprior_vlos',format='D',array=teffprior_v[targets])
    col39=fits.Column(name='teffprior_vlos_error_raw',format='D',array=teffprior_v_sig[targets])
    col40=fits.Column(name='teffprior_vlos_error',format='D',array=teffprior_v_sig[targets])
    col41=fits.Column(name='teffprior_vlos_skew',format='D',array=teffprior_v_skew[targets])
    col42=fits.Column(name='teffprior_vlos_kurtosis',format='D',array=teffprior_v_kurt[targets])
    col43=fits.Column(name='teffprior_teff',format='D',array=teffprior_teff[targets])
    col44=fits.Column(name='teffprior_teff_error_raw',format='D',array=teffprior_teff_sig[targets])
    col45=fits.Column(name='teffprior_teff_error',format='D',array=teffprior_teff_sig[targets])
    col46=fits.Column(name='teffprior_teff_skew',format='D',array=teffprior_teff_skew[targets])
    col47=fits.Column(name='teffprior_teff_kurtosis',format='D',array=teffprior_teff_kurt[targets])
    col48=fits.Column(name='teffprior_logg',format='D',array=teffprior_logg[targets])
    col49=fits.Column(name='teffprior_logg_error_raw',format='D',array=teffprior_logg_sig[targets])
    col50=fits.Column(name='teffprior_logg_error',format='D',array=teffprior_logg_sig[targets])
    col51=fits.Column(name='teffprior_logg_skew',format='D',array=teffprior_logg_skew[targets])
    col52=fits.Column(name='teffprior_logg_kurtosis',format='D',array=teffprior_logg_kurt[targets])
    col53=fits.Column(name='teffprior_feh',format='D',array=teffprior_z[targets])
    col54=fits.Column(name='teffprior_feh_error_raw',format='D',array=teffprior_z_sig[targets])
    col55=fits.Column(name='teffprior_feh_error',format='D',array=teffprior_z_sig[targets])
    col56=fits.Column(name='teffprior_feh_skew',format='D',array=teffprior_z_skew[targets])
    col57=fits.Column(name='teffprior_feh_kurtosis',format='D',array=teffprior_z_kurt[targets])
    col58=fits.Column(name='teffprior_mgfe',format='D',array=teffprior_mgfe[targets])
    col59=fits.Column(name='teffprior_mgfe_error_raw',format='D',array=teffprior_mgfe_sig[targets])
    col60=fits.Column(name='teffprior_mgfe_error',format='D',array=teffprior_mgfe_sig[targets])
    col61=fits.Column(name='teffprior_mgfe_skew',format='D',array=teffprior_mgfe_skew[targets])
    col62=fits.Column(name='teffprior_mgfe_kurtosis',format='D',array=teffprior_mgfe_kurt[targets])
    col63=fits.Column(name='teffprior_smooth',format='D',array=teffprior_smooth[targets])
    col64=fits.Column(name='teffprior_smooth_error_raw',format='D',array=teffprior_smooth_sig[targets])
    col65=fits.Column(name='teffprior_smooth_error',format='D',array=teffprior_smooth_sig[targets])
    col66=fits.Column(name='teffprior_smooth_skew',format='D',array=teffprior_smooth_skew[targets])
    col67=fits.Column(name='teffprior_smooth_kurtosis',format='D',array=teffprior_smooth_kurt[targets])
    col68=fits.Column(name='teffprior_survey',format='A10',array=teffprior_survey[targets])
    col69=fits.Column(name='median_sky',format='D',array=postfit_mediansky[targets])
    col70=fits.Column(name='standard_deviation_median_sky',format='D',array=postfit_stdsky[targets])
    col71=fits.Column(name='filtername',format='A100',array=postfit_filtername[targets])
    col72=fits.Column(name='chi2',format='D',array=postfit_chi2[targets])
    col73=fits.Column(name='chi2_rescaled',format='D',array=postfit_chi22[targets])
    col74=fits.Column(name='n',format='D',array=postfit_n[targets])
    col75=fits.Column(name='run_id',format='A100',array=postfit_run_id[targets])
    col76=fits.Column(name='field_name',format='A100',array=postfit_field_name[targets])
    col77=fits.Column(name='wav_npoints',format='PI()',array=postfit_wav_npoints[targets])
    col78=fits.Column(name='wav_rms',format='PD()',array=postfit_wav_rms[targets])
    col79=fits.Column(name='wav_resolution',format='PD()',array=postfit_wav_resolution[targets])
    col80=fits.Column(name='wav_min',format='PD()',array=postfit_wav_min[targets])
    col81=fits.Column(name='wav_max',format='PD()',array=postfit_wav_max[targets])
    col82=fits.Column(name='row',format='D',array=postfit_row[targets])
    col83=fits.Column(name='temperature',format='PD()',array=postfit_temperature[targets])
    col84=fits.Column(name='vhelio_correction',format='D',array=postfit_vhelio_correction[targets])
    col85=fits.Column(name='fits_filename',format='A200',array=postfit_filename[targets])
    col86=fits.Column(name='fits_index',format='I',array=postfit_index[targets])
    col87=fits.Column(name='gaia_source_id',format='A30',array=postfit_gaiaid[targets])
    col88=fits.Column(name='gaia_gmag',format='D',array=postfit_gaiagmag[targets])
    col89=fits.Column(name='gaia_siggmag',format='D',array=postfit_siggaiagmag[targets])
    col90=fits.Column(name='gaia_gmag_dered',format='D',array=postfit_gaiagmag_dered[targets])
    col91=fits.Column(name='gaia_bpmag',format='D',array=postfit_gaiabpmag[targets])
    col92=fits.Column(name='gaia_sigbpmag',format='D',array=postfit_siggaiabpmag[targets])
    col93=fits.Column(name='gaia_bpmag_dered',format='D',array=postfit_gaiabpmag_dered[targets])
    col94=fits.Column(name='gaia_rpmag',format='D',array=postfit_gaiarpmag[targets])
    col95=fits.Column(name='gaia_sigrpmag',format='D',array=postfit_siggaiarpmag[targets])
    col96=fits.Column(name='gaia_rpmag_dered',format='D',array=postfit_gaiarpmag_dered[targets])
    col97=fits.Column(name='gaia_pmra',format='D',array=postfit_gaiapmra[targets])
    col98=fits.Column(name='gaia_sigpmra',format='D',array=postfit_siggaiapmra[targets])
    col99=fits.Column(name='gaia_pmdec',format='D',array=postfit_gaiapmdec[targets])
    col100=fits.Column(name='gaia_sigpmdec',format='D',array=postfit_siggaiapmdec[targets])
    col101=fits.Column(name='gaia_parallax',format='D',array=postfit_gaiaparallax[targets])
    col102=fits.Column(name='gaia_sigparallax',format='D',array=postfit_siggaiaparallax[targets])
    col103=fits.Column(name='gaia_rrl',format='I',array=postfit_gaiarrl[targets])
    col104=fits.Column(name='bluechi2',format='D',array=postfit_bluechi2[targets])
    col105=fits.Column(name='redchi2',format='D',array=postfit_redchi2[targets])
    col106=fits.Column(name='bluen',format='D',array=postfit_bluen[targets])
    col107=fits.Column(name='redn',format='D',array=postfit_redn[targets])
    cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34,col35,col36,col37,col38,col39,col40,col41,col42,col43,col44,col45,col46,col47,col48,col49,col50,col51,col52,col53,col54,col55,col56,col57,col58,col59,col60,col61,col62,col63,col64,col65,col66,col67,col68,col69,col70,col71,col72,col73,col74,col75,col76,col77,col78,col79,col80,col81,col82,col83,col84,col85,col86,col87,col88,col89,col90,col91,col92,col93,col94,col95,col96,col97,col98,col99,col100,col101,col102,col103,col104,col105,col106,col107])
    t=fits.BinTableHDU.from_columns(cols)
    #table_hdu=fits.FITS_rec.from_columns(cols)
    t.writeto(fits_table_raw_filename,overwrite=True)

if get_errors:
    
    g0=open('m2fs_error_adjust.tex','w')
    pairs=fits.open('m2fs_pairs.fits')
    catalog=fits.open(fits_table_raw_filename)
    order=np.argsort(catalog[1].data['hjd'])
    catalog_sorted=catalog[1].data[order]
    gaiaid=np.array(catalog_sorted['gaia_source_id'],dtype='int')
    keep=np.where((catalog_sorted['vlos']==catalog_sorted['vlos'])&(catalog_sorted['vlos_error_raw']<5.)&(catalog_sorted['filtername']=='HiRes')&(catalog_sorted['sn_ratio']>0.)&(gaiaid>0.)&(catalog_sorted['gaia_rrl']==0))[0]
    keep=np.where((catalog_sorted['vlos']==catalog_sorted['vlos'])&(catalog_sorted['vlos_error_raw']<5.)&(catalog_sorted['filtername']=='HiRes')&(catalog_sorted['sn_ratio']>0.)&(gaiaid>0.))[0]
    catalog_keep=catalog_sorted[keep]

    used=np.zeros(len(catalog_keep),dtype='int')
    count=np.zeros(len(catalog_keep),dtype='int')
    obs=np.zeros(len(catalog_keep),dtype='int')
    nobs=np.zeros(len(catalog_keep),dtype='int')

    count0=0
    for i in range(0,len(catalog_keep)):
        if used[i]==0:
            count0+=1
            dist=np.sqrt((1./np.cos(catalog_keep['dec_deg']*np.pi/180.)*(catalog_keep['ra_deg']-catalog_keep['ra_deg'][i]))**2+(catalog_keep['dec_deg']-catalog_keep['dec_deg'][i])**2)*3600.
            this=np.where(dist<1.)[0]
            used[this]=1
            nobs[this]=len(this)
            count[this]=count0
            for j in range(0,len(this)):
                obs[this[j]]=j+1

    repeat_ind1=[]
    repeat_ind2=[]
    for i in range(0,len(catalog_keep)):
        if ((obs[i]==1)&(nobs[i]>1)):
            this=np.where(count==count[i])[0]
            for j in range(0,len(this)):
                for k in range(j+1,len(this)):
                    repeat_ind1.append(this[j])
                    repeat_ind2.append(this[k])
    repeat_ind1=np.array(repeat_ind1)
    repeat_ind2=np.array(repeat_ind2)
    string='\\newcommand{\mtwofsnrepeatv}{'+str(len(repeat_ind1))+'} \n'
    g0.write(string)
    
    thing_newcommand1=['mtwofssigvfloor','mtwofsteffsigvfloor','mtwofssigtefffloor','mtwofsteffsigtefffloor','mtwofssigloggfloor','mtwofsteffsigloggfloor','mtwofssigfehfloor','mtwofsteffsigfehfloor','mtwofssigmgfefloor','mtwofsteffsigmgfefloor','mtwofssigresolutionfloor','mtwofsteffsigresolutionfloor']
    thing_newcommand2=['mtwofssigvscale','mtwofsteffsigvscale','mtwofssigteffscale','mtwofsteffsigteffscale','mtwofssigloggscale','mtwofsteffsigloggscale','mtwofssigfehscale','mtwofsteffsigfehscale','mtwofssigmgfescale','mtwofsteffsigmgfescale','mtwofssigresolutionscale','mtwofsteffsigresolutionscale']
    thing_newcommand3=['mtwofssigvmixfrac','mtwofsteffsigvmixfrac','mtwofssigteffmixfrac','mtwofsteffsigteffmixfrac','mtwofssigloggmixfrac','mtwofsteffsigloggmixfrac','mtwofssigfehmixfrac','mtwofsteffsigfehmixfrac','mtwofssigmgfemixfrac','mtwofsteffsigmgfemixfrac','mtwofssigresolutionmixfrac','mtwofsteffsigresolutionmixfrac']
    thing_newcommand4=['mtwofssigvsigout','mtwofsteffsigvsigout','mtwofssigteffsigout','mtwofsteffsigteffsigout','mtwofssigloggsigout','mtwofsteffsigloggsigout','mtwofssigfehsigout','mtwofsteffsigfehsigout','mtwofssigmgfesigout','mtwofsteffsigmgfesigout','mtwofssigresolutionsigout','mtwofsteffsigresolutionsigout']
    thing_newcommand5=['mtwofsvnpairs','mtwofsteffvnpairs','mtwofsteffnpairs','mtwofsteffteffnpairs','mtwofsloggnpairs','mtwofsteffloggnpairs','mtwofsfehnpairs','mtwofstefffehnpairs','mtwofsmgfenpairs','mtwofsteffmgfenpairs','mtwofsresolutionnpairs','mtwofsteffresolutionnpairs']
    thing_dmax=np.array([100.,100.,2000.,2000.,2.5,2.5,2.5,2.5,1.,1.,0.03,0.03])
    floor=np.array([[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.],[-2.,2.]])
    scale=np.array([[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.]])
    frac=np.array([[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.],[0.51,1.]])
    sigmaout=np.array([[-1.,4.],[-1.,4.],[-1.,4.],[-1.,4.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.],[-1.,1.]])
    for k in range(0,len(thing)):
        x=catalog_keep[thing[k]]
        sigx=catalog_keep[thing[k]+'_error_raw']
        dx=x[repeat_ind1]-x[repeat_ind2]
        keep=np.where(np.abs(dx)<=thing_dmax[k])[0]
        result,bestfit=m2fs.get_error_adjust(x,sigx,floor[k],scale[k],frac[k],sigmaout[k],repeat_ind1,repeat_ind2,thing[k],thing_symbol[k],thing_dmax[k])
        pickle.dump((x,sigx,repeat_ind1,repeat_ind2,thing[k],thing_symbol[k],thing_dmax[k],result,bestfit),open('m2fs_ian_'+thing[k]+'_adjust.pkl','wb'))
        string='\\newcommand{\\'+thing_newcommand1[k]+'}{$'+str('{0:.2f}'.format(np.mean(10.**result['samples'].T[0])))+' \pm '+str('{0:.2f}'.format(np.std(10.**result['samples'].T[0])))+'$} \n'
        g0.write(string)
        string='\\newcommand{\\'+thing_newcommand2[k]+'}{$'+str('{0:.2f}'.format(np.mean(10.**result['samples'].T[1])))+' \pm '+str('{0:.2f}'.format(np.std(10.**result['samples'].T[1])))+'$} \n'
        g0.write(string)
        string='\\newcommand{\\'+thing_newcommand3[k]+'}{$'+str('{0:.2f}'.format(1.-np.mean(result['samples'].T[2])))+' \pm '+str('{0:.2f}'.format(np.std(result['samples'].T[2])))+'$} \n'
        g0.write(string)
        string='\\newcommand{\\'+thing_newcommand4[k]+'}{$'+str('{0:.2f}'.format(np.mean(10.**result['samples'].T[3])))+' \pm '+str('{0:.2f}'.format(np.std(10.**result['samples'].T[3])))+'$} \n'
        g0.write(string)
        string='\\newcommand{\\'+thing_newcommand5[k]+'}{$'+str(len(keep))+'$} \n'
        g0.write(string)
    g0.close()
    
if get_catalog_final:

    catalog_raw=fits.open(fits_table_raw_filename)
    catalog=catalog_raw[1].data
    for k in range(0,len(thing)):
        crap,crap,crap,crap,crap,crap,crap,samples,bestfit=pickle.load(open('m2fs_ian_'+thing[k]+'_adjust.pkl','rb'))
        catalog[thing[k]+'_error']=10.**bestfit['parameters'][0]+10.**bestfit['parameters'][1]*catalog_raw[1].data[thing[k]+'_error_raw']

    used=np.zeros(len(catalog),dtype='int')
    count=np.zeros(len(catalog),dtype='int')
    obs=np.zeros(len(catalog),dtype='int')
    goodobs=np.zeros(len(catalog),dtype='int')
    goodvobs=np.zeros(len(catalog),dtype='int')
    nobs=np.zeros(len(catalog),dtype='int')
    goodnobs=np.zeros(len(catalog),dtype='int')
    count0=0

    v_mean=np.zeros(len(catalog))
    teff_mean=np.zeros(len(catalog))
    logg_mean=np.zeros(len(catalog))
    z_mean=np.zeros(len(catalog))
    mgfe_mean=np.zeros(len(catalog))
    sigv_mean=np.zeros(len(catalog))
    sigteff_mean=np.zeros(len(catalog))
    siglogg_mean=np.zeros(len(catalog))
    sigfeh_mean=np.zeros(len(catalog))
    sigmgfe_mean=np.zeros(len(catalog))
    stdv_mean=np.zeros(len(catalog))
    stdteff_mean=np.zeros(len(catalog))
    stdlogg_mean=np.zeros(len(catalog))
    stdz_mean=np.zeros(len(catalog))
    stdmgfe_mean=np.zeros(len(catalog))

    teffprior_v_mean=np.zeros(len(catalog))
    teffprior_teff_mean=np.zeros(len(catalog))
    teffprior_logg_mean=np.zeros(len(catalog))
    teffprior_z_mean=np.zeros(len(catalog))
    teffprior_mgfe_mean=np.zeros(len(catalog))
    teffprior_sigv_mean=np.zeros(len(catalog))
    teffprior_sigteff_mean=np.zeros(len(catalog))
    teffprior_siglogg_mean=np.zeros(len(catalog))
    teffprior_sigfeh_mean=np.zeros(len(catalog))
    teffprior_sigmgfe_mean=np.zeros(len(catalog))
    teffprior_stdv_mean=np.zeros(len(catalog))
    teffprior_stdteff_mean=np.zeros(len(catalog))
    teffprior_stdlogg_mean=np.zeros(len(catalog))
    teffprior_stdz_mean=np.zeros(len(catalog))
    teffprior_stdmgfe_mean=np.zeros(len(catalog))

    for i in range(0,len(catalog)):
        if used[i]==0:
            count0+=1
            dist=np.sqrt((1./np.cos(catalog['dec_deg']*np.pi/180.)*(catalog['ra_deg']-catalog['ra_deg'][i]))**2+(catalog['dec_deg']-catalog['dec_deg'][i])**2)*3600.
            this=np.where(dist<1.)[0]
            this_goodv=np.where((catalog['sn_ratio']>0.)&(catalog['vlos']==catalog['vlos'])&(dist<1.)&(catalog['vlos_error_raw']<5.))[0]
            this_goodteff=np.where((catalog['sn_ratio']>0.)&(catalog['vlos']==catalog['vlos'])&(dist<1.)&(catalog['vlos_error_raw']<5.)&(catalog['teff_error_raw']<100000.))[0]
            this_goodlogg=np.where((catalog['sn_ratio']>0.)&(catalog['vlos']==catalog['vlos'])&(dist<1.)&(catalog['vlos_error_raw']<5.)&(catalog['logg_error_raw']<1000.5))[0]
            this_goodz=np.where((catalog['sn_ratio']>0.)&(catalog['vlos']==catalog['vlos'])&(dist<1.)&(catalog['vlos_error_raw']<5.)&(catalog['logg_error_raw']<1000.5)&(catalog['feh_error_raw']<1000.5))[0]
            this_goodmgfe=np.where((catalog['sn_ratio']>0.)&(catalog['vlos']==catalog['vlos'])&(dist<1.)&(catalog['vlos_error_raw']<5.)&(catalog['logg_error_raw']<1000.5)&(catalog['feh_error_raw']<1000.5)&(catalog['mgfe_error_raw']<1000.3))[0]
            used[this]=1
            nobs[this]=len(this)
            goodnobs[this]=len(this_goodv)
            count[this]=count0
            for j in range(0,len(this)):
                obs[this[j]]=j+1        
            for j in range(0,len(this_goodv)):
                goodobs[this_goodv[j]]=j+1

            if len(this_goodv)>0:
                mean,sigmean,std=mycode.weightedmean(catalog['vlos'][this_goodv],catalog['vlos_error'][this_goodv])
                v_mean[this]=mean
                sigv_mean[this]=sigmean
                stdv_mean[this]=std
                mean,sigmean,std=mycode.weightedmean(catalog['teffprior_vlos'][this_goodv],catalog['teffprior_vlos_error'][this_goodv])
                teffprior_v_mean[this]=mean
                teffprior_sigv_mean[this]=sigmean
                teffprior_stdv_mean[this]=std
            else:
                v_mean[this]=np.nan
                sigv_mean[this]=np.nan
                stdv_mean[this]=np.nan
                teffprior_v_mean[this]=np.nan
                teffprior_sigv_mean[this]=np.nan
                teffprior_stdv_mean[this]=np.nan
                
            if len(this_goodteff)>0:
                mean,sigmean,std=mycode.weightedmean(catalog['teff'][this_goodteff],catalog['teff_error'][this_goodteff])
                teff_mean[this]=mean
                sigteff_mean[this]=sigmean
                stdteff_mean[this]=std
                mean,sigmean,std=mycode.weightedmean(catalog['teffprior_teff'][this_goodteff],catalog['teffprior_teff_error'][this_goodteff])
                teffprior_teff_mean[this]=mean
                teffprior_sigteff_mean[this]=sigmean
                teffprior_stdteff_mean[this]=std
            else:
                teff_mean[this]=np.nan
                sigteff_mean[this]=np.nan
                stdteff_mean[this]=np.nan
                teffprior_teff_mean[this]=np.nan
                teffprior_sigteff_mean[this]=np.nan
                teffprior_stdteff_mean[this]=np.nan

            if len(this_goodlogg)>0:
                mean,sigmean,std=mycode.weightedmean(catalog['logg'][this_goodlogg],catalog['logg_error'][this_goodlogg])
                logg_mean[this]=mean
                siglogg_mean[this]=sigmean
                stdlogg_mean[this]=std
                mean,sigmean,std=mycode.weightedmean(catalog['teffprior_logg'][this_goodlogg],catalog['teffprior_logg_error'][this_goodlogg])
                teffprior_logg_mean[this]=mean
                teffprior_siglogg_mean[this]=sigmean
                teffprior_stdlogg_mean[this]=std
            else:
                logg_mean[this]=np.nan
                siglogg_mean[this]=np.nan
                stdlogg_mean[this]=np.nan
                teffprior_logg_mean[this]=np.nan
                teffprior_siglogg_mean[this]=np.nan
                teffprior_stdlogg_mean[this]=np.nan

            if len(this_goodz)>0:
                mean,sigmean,std=mycode.weightedmean(catalog['feh'][this_goodz],catalog['feh_error'][this_goodz])
                z_mean[this]=mean
                sigfeh_mean[this]=sigmean
                stdz_mean[this]=std
                mean,sigmean,std=mycode.weightedmean(catalog['teffprior_feh'][this_goodz],catalog['teffprior_feh_error'][this_goodz])
                teffprior_z_mean[this]=mean
                teffprior_sigfeh_mean[this]=sigmean
                teffprior_stdz_mean[this]=std
            else:
                z_mean[this]=np.nan
                sigfeh_mean[this]=np.nan
                stdz_mean[this]=np.nan
                teffprior_z_mean[this]=np.nan
                teffprior_sigfeh_mean[this]=np.nan
                teffprior_stdz_mean[this]=np.nan

            if len(this_goodmgfe)>0:
                mean,sigmean,std=mycode.weightedmean(catalog['mgfe'][this_goodmgfe],catalog['mgfe_error'][this_goodmgfe])
                mgfe_mean[this]=mean
                sigmgfe_mean[this]=sigmean
                stdmgfe_mean[this]=std
                mean,sigmean,std=mycode.weightedmean(catalog['teffprior_mgfe'][this_goodmgfe],catalog['teffprior_mgfe_error'][this_goodmgfe])
                teffprior_mgfe_mean[this]=mean
                teffprior_sigmgfe_mean[this]=sigmean
                teffprior_stdmgfe_mean[this]=std
            else:
                mgfe_mean[this]=np.nan
                sigmgfe_mean[this]=np.nan
                stdmgfe_mean[this]=np.nan
                teffprior_mgfe_mean[this]=np.nan
                teffprior_sigmgfe_mean[this]=np.nan
                teffprior_stdmgfe_mean[this]=np.nan

    col1=fits.Column(name='vlos_mean',format='D',array=v_mean)
    col2=fits.Column(name='vlos_mean_error',format='D',array=sigv_mean)
    col3=fits.Column(name='vlos_mean_scatter',format='D',array=stdv_mean)
    col4=fits.Column(name='teff_mean',format='D',array=teff_mean)
    col5=fits.Column(name='teff_mean_error',format='D',array=sigteff_mean)
    col6=fits.Column(name='teff_mean_scatter',format='D',array=stdteff_mean)
    col7=fits.Column(name='logg_mean',format='D',array=logg_mean)
    col8=fits.Column(name='logg_mean_error',format='D',array=siglogg_mean)
    col9=fits.Column(name='logg_mean_scatter',format='D',array=stdlogg_mean)
    col10=fits.Column(name='feh_mean',format='D',array=z_mean)
    col11=fits.Column(name='feh_mean_error',format='D',array=sigfeh_mean)
    col12=fits.Column(name='feh_mean_scatter',format='D',array=stdz_mean)
    col13=fits.Column(name='mgfe_mean',format='D',array=mgfe_mean)
    col14=fits.Column(name='mgfe_mean_error',format='D',array=sigmgfe_mean)
    col15=fits.Column(name='mgfe_mean_scatter',format='D',array=stdmgfe_mean)

    col16=fits.Column(name='teffprior_vlos_mean',format='D',array=teffprior_v_mean)
    col17=fits.Column(name='teffprior_vlos_mean_error',format='D',array=teffprior_sigv_mean)
    col18=fits.Column(name='teffprior_vlos_mean_scatter',format='D',array=teffprior_stdv_mean)
    col19=fits.Column(name='teffprior_teff_mean',format='D',array=teffprior_teff_mean)
    col20=fits.Column(name='teffprior_teff_mean_error',format='D',array=teffprior_sigteff_mean)
    col21=fits.Column(name='teffprior_teff_mean_scatter',format='D',array=teffprior_stdteff_mean)
    col22=fits.Column(name='teffprior_logg_mean',format='D',array=teffprior_logg_mean)
    col23=fits.Column(name='teffprior_logg_mean_error',format='D',array=teffprior_siglogg_mean)
    col24=fits.Column(name='teffprior_logg_mean_scatter',format='D',array=teffprior_stdlogg_mean)
    col25=fits.Column(name='teffprior_feh_mean',format='D',array=teffprior_z_mean)
    col26=fits.Column(name='teffprior_feh_mean_error',format='D',array=teffprior_sigfeh_mean)
    col27=fits.Column(name='teffprior_feh_mean_scatter',format='D',array=teffprior_stdz_mean)
    col28=fits.Column(name='teffprior_mgfe_mean',format='D',array=teffprior_mgfe_mean)
    col29=fits.Column(name='teffprior_mgfe_mean_error',format='D',array=teffprior_sigmgfe_mean)
    col30=fits.Column(name='teffprior_mgfe_mean_scatter',format='D',array=teffprior_stdmgfe_mean)

    col31=fits.Column(name='obs',format='I',array=obs)
    col32=fits.Column(name='n_obs',format='I',array=nobs)
    col33=fits.Column(name='good_obs',format='I',array=goodobs)
    col34=fits.Column(name='good_n_obs',format='I',array=goodnobs)

    raw_cols=catalog.columns
    new_cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34])
    cols=raw_cols+new_cols    
    table_hdu=fits.BinTableHDU.from_columns(cols)#BinTableHDU.from_columns(cols)
        
    table_hdu.writeto(fits_table_final_filename,overwrite=True)

    catalog=fits.open('all_m2fshiresian_final.fits')[1].data

    g1=open('all_m2fshiresian_final.dat','w')
    g1.write('#Object # RA [deg] # Dec [deg] # HJD [days] # median S/N/pix # vhelio [km/s] # err_vhelio [km/s] # Teff [K] # err_Teff [K]  # logg # err_logg # Fe/H # err_Fe/H # [Mg/Fe] # err_[Mg/Fe] # vhelio_mean [km/s] # err_vhelio_mean [km/s] # Teff_mean [K] # err_Teff_mean [K]  # logg_mean # err_logg_mean # Fe/H_mean # err_Fe/H_mean # [Mg/Fe]_mean # err_[Mg/Fe]_mean # good_observation_index   \n ')
    for i in range(0,len(catalog)):
            if catalog['good_n_obs'][i]>0:
                    string=catalog['target_system'][i]+' '+str.format('{0:.7f}',round(catalog['ra_deg'][i],7)).zfill(7)+' '+str.format('{0:.7f}',round(catalog['dec_deg'][i],7)).zfill(7)+' '+str.format('{0:.3f}',round(catalog['hjd'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['sn_ratio'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['vlos'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['vlos_error'][i],3)).zfill(3)+' '+str.format('{0:.1f}',round(catalog['teff'][i],1)).zfill(1)+' '+str.format('{0:.1f}',round(catalog['teff_error'][i],1)).zfill(1)+' '+str.format('{0:.3f}',round(catalog['logg'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['logg_error'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['feh'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['feh_error'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['mgfe'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['mgfe_error'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['vlos_mean'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['vlos_mean_error'][i],3)).zfill(3)+' '+str.format('{0:.1f}',round(catalog['teff_mean'][i],1)).zfill(1)+' '+str.format('{0:.1f}',round(catalog['teff_mean_error'][i],1)).zfill(1)+' '+str.format('{0:.3f}',round(catalog['logg_mean'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['logg_mean_error'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['feh_mean'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['feh_mean_error'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['mgfe_mean'][i],3)).zfill(3)+' '+str.format('{0:.3f}',round(catalog['mgfe_mean_error'][i],3)).zfill(3)+' '+str(catalog['good_obs'][i])+' '+' \n'
                    g1.write(string)
    g1.close()

    
if check_badfield:

    m2fs=fits.open(fits_table_final_filename)[1].data

    obs=np.zeros(len(m2fs),dtype='int')
    count=0
    rate=[]
    for i in range(0,len(m2fs)):
        if obs[i]==0:
            this=np.where(m2fs['fits_filename']==m2fs['fits_filename'][i])[0]
            if len(this)>0:
                obs[this]=1
                keep=np.where(m2fs['good_obs'][this]==1)[0]
                if len(keep)>0:
                    count+=1
                    for j in range(0,len(keep)):
                        stat=m2fs['vlos_mean_scatter'][this][keep]/m2fs['vlos_mean_error'][this][keep]
#                        stat=m2fs['vlos_mean_scatter'][this][keep]
                    bad=np.where(stat>10.)[0]
                    rate0=len(bad)/len(keep)
                    rate.append(rate0)
                    if rate0>0.15:
                        print(i,count,m2fs['fits_filename'][i],rate0)
    rate=np.array(rate)
                        

                    
                              
    
    
if get_plots:

    data_table=open('m2fs_data_table.tex','w')
    data_table_short=open('m2fs_data_table_short.tex','w')
    field_center=pickle.load(open('m2fs_field_center.pkl','rb'))
    labels=np.array(['name','host','ra','dec','rhalf','sigrhalf','pa','sigpa','ellipticity','sigellipticity','ref_structure','dmodulus','sigdmodulus','distance','sigdistance','ref_distance','rhalfpc','rhalfpcsphere','absvmag','sigabsvmag','ref_abvsmag','vlos','sigvlos','vdisp','sigvdisp','vdisplimit','ref_vlos','pmra_dr3','sigpmra_dr3','pmdec_dr3','sigpmdec_dr3','ref_pm_dr3','pmra_dr2','sigpmra_dr2','pmdec_dr2','sigpmdec_dr2','ref_pm_dr2','feh','sigfeh','fehdisp','sigfehdisp','fehdisplimit','ref_z'],dtype='str')
    cols=[]
    for i in range(0,len(labels)):
        in_csvfile=open('general_dsph_info_all.csv','r',newline='')
        in_obj=csv.reader(in_csvfile)
        xxx=[]
        for row in in_obj:
            if in_obj.line_num>1:
                xxx.append(row[i])
        xxx=np.array(xxx)
        if ((not 'ref' in labels[i])&(not 'name' in labels[i])&(not 'host' in labels[i])):
            bad=np.where(xxx=='')[0]
            xxx[bad]=np.nan
            cols.append(fits.Column(name=labels[i],format='D',array=xxx))
        else:
            bad=np.where(xxx=='')[0]
            xxx[bad]='none'
            cols.append(fits.Column(name=labels[i],format='A100',array=xxx))
    coldefs=fits.ColDefs([cols[q] for q in range(0,len(cols))])
    table_hdu=fits.BinTableHDU.from_columns([cols[q] for q in range(0,len(cols))])
    table_hdu.writeto('dsph_parameters.fits',overwrite=True)
    dsph=fits.open('dsph_parameters.fits')
#    dsph=fits.open('NearbyGalaxies_Jan2021_PUBLIC.fits')

    catalog=fits.open(fits_table_final_filename)[1].data
    coords=SkyCoord(catalog['ra_deg'],catalog['dec_deg'],unit=(u.deg,u.deg))

    age=np.zeros(len(dsph[1].data['name']))+1.e+10
    
    field_object_andrew=np.empty(len(catalog['target_system']),dtype='object')
    field_object_formal=np.empty(len(catalog['target_system']),dtype='object')
    this=np.where(catalog['target_system']=='ant2')[0]
    field_object_andrew[this]='antlia_2'
    field_object_formal[this]='Antlia 2'
    this=np.where(catalog['target_system']=='car')[0]
    field_object_andrew[this]='carina_1'
    field_object_formal[this]='Carina'
    this=np.where(catalog['target_system']=='cra')[0]
    field_object_andrew[this]='crater'
    field_object_formal[this]='Crater'
    this=np.where(catalog['target_system']=='cra2')[0]
    field_object_andrew[this]='crater_2'
    field_object_formal[this]='Crater 2'
    this=np.where(catalog['target_system']=='for')[0]
    field_object_andrew[this]='fornax_1'
    field_object_formal[this]='Fornax'
    this=np.where(catalog['target_system']=='gru1')[0]
    field_object_andrew[this]='grus_1'
    field_object_formal[this]='Grus 1'
    this=np.where(catalog['target_system']=='gru2')[0]
    field_object_andrew[this]='grus_2'
    field_object_formal[this]='Grus 2'
    this=np.where(catalog['target_system']=='hor2')[0]
    field_object_andrew[this]='horologium_2'
    field_object_formal[this]='Horologium 2'
    this=np.where(catalog['target_system']=='hyd1')[0]
    field_object_andrew[this]='hydrus_1'
    field_object_formal[this]='Hydrus 1'
    this=np.where(catalog['target_system']=='ind1')[0]
    field_object_andrew[this]='indus_1'
    field_object_formal[this]='Indus 1'
    this=np.where(catalog['target_system']=='ind2')[0]
    field_object_andrew[this]='indus_2'
    field_object_formal[this]='Indus 2'
    this=np.where(catalog['target_system']=='kgo2')[0]
    field_object_andrew[this]='kgo2'
    field_object_formal[this]='KGO 2'
    this=np.where(catalog['target_system']=='kgo4')[0]
    field_object_andrew[this]='kgo4'
    field_object_formal[this]='KGO 4'
    this=np.where(catalog['target_system']=='kgo7')[0]
    field_object_andrew[this]='kgo7'
    field_object_formal[this]='KGO 7'
    this=np.where(catalog['target_system']=='kgo8')[0]
    field_object_andrew[this]='kgo8'
    field_object_formal[this]='KGO 8'
    this=np.where(catalog['target_system']=='kgo10')[0]
    field_object_andrew[this]='kgo10'
    field_object_formal[this]='KGO 10'
    this=np.where(catalog['target_system']=='kgo13')[0]
    field_object_andrew[this]='kgo13'
    field_object_formal[this]='KGO 13'
    this=np.where(catalog['target_system']=='kgo22')[0]
    field_object_andrew[this]='kgo22'
    field_object_formal[this]='KGO 22'
    this=np.where(catalog['target_system']=='kop2')[0]
    field_object_andrew[this]='koposov_2'
    field_object_formal[this]='Koposov 2'
    this=np.where(catalog['target_system']=='pal5')[0]
    field_object_andrew[this]='palomar_5'
    field_object_formal[this]='Palomar 5'
    this=np.where(catalog['target_system']=='pho2')[0]
    field_object_andrew[this]='phoenix_2'
    field_object_formal[this]='Phoenix 2'
    this=np.where(catalog['target_system']=='ret2')[0]
    field_object_andrew[this]='reticulum_2'
    field_object_formal[this]='Reticulum 2'
    this=np.where(catalog['target_system']=='sgr2')[0]
    field_object_andrew[this]='sagittarius_2'
    field_object_formal[this]='Sagittarius 2'
    this=np.where(catalog['target_system']=='scl')[0]
    field_object_andrew[this]='sculptor_1'
    field_object_formal[this]='Sculptor'
    this=np.where(catalog['target_system']=='sex')[0]
    field_object_andrew[this]='sextans_1'
    field_object_formal[this]='Sextans'
    this=np.where(catalog['target_system']=='tuc2')[0]
    field_object_andrew[this]='tucana_2'
    field_object_formal[this]='Tucana 2'
    this=np.where(catalog['target_system']=='tuc3')[0]
    field_object_andrew[this]='tucana_3'
    field_object_formal[this]='Tucana 3'
    this=np.where(catalog['target_system']=='tuc4')[0]
    field_object_andrew[this]='tucana_4'
    field_object_formal[this]='Tucana 4'
    this=np.where(catalog['target_system']=='tuc5')[0]
    field_object_andrew[this]='tucana_5'
    field_object_formal[this]='Tucana 5'

    order=np.argsort(catalog['hjd'])
    catalog_sorted=catalog[order]
    keep=np.where((catalog_sorted['vlos_error']<5.)&(catalog_sorted['sn_ratio']>0.))[0]
    catalog_keep=catalog_sorted[keep]
    field_object_formal_keep=field_object_formal[order][keep]

    used=np.zeros(len(catalog_keep),dtype='int')
    count=np.zeros(len(catalog_keep),dtype='int')
    obs=np.zeros(len(catalog_keep),dtype='int')
    nobs=np.zeros(len(catalog_keep),dtype='int')

    count0=0
    for i in range(0,len(catalog_keep)):
        if used[i]==0:
            count0+=1
            dist=np.sqrt((1./np.cos(catalog_keep['dec_deg']*np.pi/180.)*(catalog_keep['ra_deg']-catalog_keep['ra_deg'][i]))**2+(catalog_keep['dec_deg']-catalog_keep['dec_deg'][i])**2)*3600.
            this=np.where(dist<1.)[0]
            used[this]=1
            nobs[this]=len(this)
            count[this]=count0
            for j in range(0,len(this)):
                obs[this[j]]=j+1

    nshort=0

    catalog_ra=[]
    catalog_dec=[]
    catalog_gaia=[]
    catalog_hjd=[]
    catalog_vlos=[]
    catalog_vlos_error=[]
    catalog_teff=[]
    catalog_teff_error=[]
    catalog_logg=[]
    catalog_logg_error=[]
    catalog_z=[]
    catalog_z_error=[]
    catalog_mgfe=[]
    catalog_mgfe_error=[]
    catalog_obj=[]
    
    for i in range(0,len(catalog_keep)):
        if obs[i]==1:
            dist=np.sqrt((1./np.cos(catalog_keep['dec_deg']*np.pi/180.)*(catalog_keep['ra_deg']-catalog_keep['ra_deg'][i]))**2+(catalog_keep['dec_deg']-catalog_keep['dec_deg'][i])**2)*3600.
            this=np.where(dist<1.)[0]
            if catalog_keep['ra_deg'][i]<10.:
                rastring='$00'+str('{0:.6f}'.format(catalog_keep['ra_deg'][i]))+'$'
            elif catalog_keep['ra_deg'][i]<100.:
                rastring='$0'+str('{0:.6f}'.format(catalog_keep['ra_deg'][i]))+'$'
            else:
                rastring='$'+str('{0:.6f}'.format(catalog_keep['ra_deg'][i]))+'$'
            if catalog_keep['dec_deg'][i]<0:
                if np.abs(catalog_keep['dec_deg'][i])<10.:
                    decstring='$-00'+str('{0:.6f}'.format(np.abs(catalog_keep['dec_deg'][i])))+'$'
                elif np.abs(catalog_keep['dec_deg'][i])<100.:
                    decstring='$-0'+str('{0:.6f}'.format(np.abs(catalog_keep['dec_deg'][i])))+'$'
                else:
                    decstring='$-'+str('{0:.6f}'.format(np.abs(catalog_keep['dec_deg'][i])))+'$'
            else:
                if np.abs(catalog_keep['dec_deg'][i])<10.:
                    decstring='$+00'+str('{0:.6f}]'.format(np.abs(catalog_keep['dec_deg'][i])))+'$'
                elif np.abs(catalog_keep['dec_deg'][i])<100.:
                    decstring='$+0'+str('{0:.6f}'.format(np.abs(catalog_keep['dec_deg'][i])))+'$'
                else:
                    decstring='$+'+str('{0:.6f}'.format(np.abs(catalog_keep['dec_deg'][i])))+'$'
            for j in range(0,len(this)):
                hjdstring='$'+str('{0:.3f}'.format(catalog_keep['hjd'][this[j]]))+'$'
                vstring='$'+str('{0:.1f}'.format(catalog_keep['vlos'][this[j]]))+' \pm '+str('{0:.1f}'.format(catalog_keep['vlos_error'][this[j]]))+'$'
                teffstring='$'+str('{0:.0f}'.format(catalog_keep['teff'][this[j]]))+' \pm '+str('{0:.0f}'.format(catalog_keep['teff_error'][this[j]]))+'$'
                loggstring='$'+str('{0:.2f}'.format(catalog_keep['logg'][this[j]]))+' \pm '+str('{0:.2f}'.format(catalog_keep['logg_error'][this[j]]))+'$'
                zstring='$'+str('{0:.2f}'.format(catalog_keep['feh'][this[j]]))+' \pm '+str('{0:.2f}'.format(catalog_keep['feh_error'][this[j]]))+'$'
                mgfestring='$'+str('{0:.2f}'.format(catalog_keep['mgfe'][this[j]]))+' \pm '+str('{0:.2f}'.format(catalog_keep['mgfe_error'][this[j]]))+'$'
                gaiastring='$'+catalog_keep['gaia_source_id'][this[j]]+'$'
                if np.float(catalog_keep['gaia_source_id'][this[j]])<0:
                    gaiastring='$\\nodata$'
                objstring=str(field_object_formal_keep[this[j]])
                if j==0:
                    data_table.write(rastring+'& '+decstring+'& '+gaiastring+'& '+hjdstring+'& '+vstring+'& '+teffstring+'& '+loggstring+'& '+zstring+'& '+mgfestring+'& '+objstring+'\\\ \n')
                    if nshort<20:
                        data_table_short.write(rastring+'& '+decstring+'& '+gaiastring+'& '+hjdstring+'& '+vstring+'& '+teffstring+'& '+loggstring+'& '+zstring+'& '+mgfestring+'& '+objstring+'\\\ \n')
                        nshort+=1

                    catalog_ra.append(rastring)
                    catalog_dec.append(decstring)
                    catalog_gaia.append(gaiastring)
                    catalog_hjd.append(hjdstring)
                    catalog_vlos.append(catalog_keep['vlos'][this[j]])
                    catalog_vlos_error.append(catalog_keep['vlos_error'][this[j]])
                    catalog_teff.append(catalog_keep['teff'][this[j]])
                    catalog_teff_error.append(catalog_keep['teff_error'][this[j]])
                    catalog_logg.append(catalog_keep['logg'][this[j]])
                    catalog_logg_error.append(catalog_keep['logg_error'][this[j]])
                    catalog_z.append(catalog_keep['feh'][this[j]])
                    catalog_z_error.append(catalog_keep['feh_error'][this[j]])
                    catalog_mgfe.append(catalog_keep['mgfe'][this[j]])
                    catalog_mgfe_error.append(catalog_keep['mgfe_error'][this[j]])
                    catalog_obj.append(objstring)
                    
                else:
                    data_table.write('&&& '+hjdstring+'& '+vstring+'& '+teffstring+'& '+loggstring+'& '+zstring+'& '+mgfestring+'&& '+'\\\ \n')
                    if nshort<20:
                        data_table_short.write('&&& '+hjdstring+'& '+vstring+'& '+teffstring+'& '+loggstring+'& '+zstring+'& '+mgfestring+'&& '+'\\\ \n')
                        nshort+=1
                        
                    catalog_ra.append(' ')
                    catalog_dec.append(' ')
                    catalog_gaia.append(' ')
                    catalog_hjd.append(hjdstring)
                    catalog_vlos.append(catalog_keep['vlos'][this[j]])
                    catalog_vlos_error.append(catalog_keep['vlos_error'][this[j]])
                    catalog_teff.append(catalog_keep['teff'][this[j]])
                    catalog_teff_error.append(catalog_keep['teff_error'][this[j]])
                    catalog_logg.append(catalog_keep['logg'][this[j]])
                    catalog_logg_error.append(catalog_keep['logg_error'][this[j]])
                    catalog_z.append(catalog_keep['feh'][this[j]])
                    catalog_z_error.append(catalog_keep['feh_error'][this[j]])
                    catalog_mgfe.append(catalog_keep['mgfe'][this[j]])
                    catalog_mgfe_error.append(catalog_keep['mgfe_error'][this[j]])
                    catalog_obj.append(' ')
                        
    data_table.close()
    data_table_short.close()
    pickle.dump((catalog_ra,catalog_dec,catalog_gaia,catalog_hjd,catalog_vlos,catalog_vlos_error,catalog_teff,catalog_teff_error,catalog_logg,catalog_logg_error,catalog_z,catalog_z_error,catalog_mgfe,catalog_mgfe_error,catalog_obj),open('m2fs_catalog_before_zeropointshift.pkl','wb'))
    
    
    ra=[]
    dec=[]
    rad=[]
    v=[]
    teff=[]
    logg=[]
    z=[]
    mgfe=[]
    sigv=[]
    sigteff=[]
    siglogg=[]
    sigfeh=[]
    sigmgfe=[]
    teffprior_v=[]
    teffprior_teff=[]
    teffprior_logg=[]
    teffprior_z=[]
    teffprior_mgfe=[]
    teffprior_sigv=[]
    teffprior_sigteff=[]
    teffprior_siglogg=[]
    teffprior_sigfeh=[]
    teffprior_sigmgfe=[]
    isomem=[]
    vmem=[]
    pmmem=[]
    vpmmem=[]
    isovpmmem=[]
    obj=[]
    v_mean=[]
    sigv_mean=[]
    stdv_mean=[]
    teff_mean=[]
    sigteff_mean=[]
    stdteff_mean=[]
    logg_mean=[]
    siglogg_mean=[]
    stdlogg_mean=[]
    z_mean=[]
    sigfeh_mean=[]
    stdz_mean=[]
    mgfe_mean=[]
    sigmgfe_mean=[]
    stdmgfe_mean=[]
    teffprior_v_mean=[]
    teffprior_sigv_mean=[]
    teffprior_stdv_mean=[]
    teffprior_teff_mean=[]
    teffprior_sigteff_mean=[]
    teffprior_stdteff_mean=[]
    teffprior_logg_mean=[]
    teffprior_siglogg_mean=[]
    teffprior_stdlogg_mean=[]
    teffprior_z_mean=[]
    teffprior_sigfeh_mean=[]
    teffprior_stdz_mean=[]
    teffprior_mgfe_mean=[]
    teffprior_sigmgfe_mean=[]
    teffprior_stdmgfe_mean=[]
    obs=[]
    nobs=[]
    goodobs=[]
    goodnobs=[]
    
    for i in range(0,len(dsph[1].data['name'])):
        print(i,dsph[1].data['name'][i])
    
    for i in [21]:
        age=1.e+10
        maglim=[21,14]
        maglimticks=[21,20,19,18,17,16,15,14]
        print(i,dsph[1].data['name'][i],dsph[1].data['rhalf'][i],dsph[1].data['feh'][i],dsph[1].data['ellipticity'][i],dsph[1].data['pa'][i],dsph[1].data['pmra_dr3'][i],dsph[1].data['pmdec_dr3'][i],)
        if dsph[1].data['feh'][i]!=dsph[1].data['feh'][i]:
            dsph[1].data['feh'][i]=-2.
            dsph[1].data['feh'][i]=999.
        if dsph[1].data['ellipticity'][i]>2.:
            dsph[1].data['ellipticity'][i]=0.
        if dsph[1].data['feh'][i]>2.:
            dsph[1].data['feh'][i]=-2.
        if dsph[1].data['vdisp'][i]>900.:
            dsph[1].data['vdisp'][i]=5.
        center=SkyCoord(dsph[1].data['ra'][i],dsph[1].data['dec'][i],unit=(u.deg,u.deg))
        xi,eta=mycode.etaxiarr(coords.ra.rad,coords.dec.rad,center.ra.rad,center.dec.rad)
#        this=np.where(catalog['target_system']==dsph[1].data['name'][i])[0]
        r=np.sqrt(xi**2+eta**2)
        this=np.where((field_object_andrew==dsph[1].data['name'][i])&(catalog['vlos_error']<=5.))[0]
        this1=np.where((field_object_andrew==dsph[1].data['name'][i])&(catalog['vlos_error']<=5.)&(catalog['good_obs']==1))[0]
        thismem=np.where((field_object_andrew==dsph[1].data['name'][i])&(catalog['vlos_error']<=5.)&(catalog['good_obs']==1)&(catalog['logg_mean']<=2.))[0]
        thismem3=np.where((field_object_andrew==dsph[1].data['name'][i])&(catalog['vlos_error']<=5.)&(catalog['good_obs']==1)&(catalog['logg_mean']<=2.)&(r>dsph[1].data['rhalf'][i]*3))[0]
        thismem4=np.where((field_object_andrew==dsph[1].data['name'][i])&(catalog['vlos_error']<=5.)&(catalog['good_obs']==1)&(catalog['logg_mean']<=2.)&(r>dsph[1].data['rhalf'][i]*4))[0]
        thismem5=np.where((field_object_andrew==dsph[1].data['name'][i])&(catalog['vlos_error']<=5.)&(catalog['good_obs']==1)&(catalog['logg_mean']<=2.)&(r>dsph[1].data['rhalf'][i]*5))[0]
        rmax=np.max(r[this1])/60.
        rsearch=np.min(np.array([5.,3.*rmax]))
        rplot=np.max(np.array([3.*dsph[1].data['rhalf'][i]/60.,np.max(r[this1])/60.*1.1,13./60.]))
        if dsph[1].data['name'][i]=='indus_1':
            age=1.e9
        if dsph[1].data['name'][i]=='fornax_1':
            age=5.e9
        if dsph[1].data['name'][i]=='hydrus_1':
            age=1.e10
#            maglim=[21,14]
        if dsph[1].data['name'][i]=='reticulum_2':
            age=1.3e10
#            maglim=[21,14]
        if dsph[1].data['name'][i]=='tucana_3':
            age=1.3e10
#            maglim=[21,14]
        if dsph[1].data['name'][i]=='grus_1':
            age=1.3e10
#            maglim=[21,14]
#        if dsph[1].data['name'][i]=='tucana_4':
#            maglim=[21,15]
        if dsph[1].data['name'][i]=='tucana_5':
            age=10.e+9
#            maglim=[21,15]
        if dsph[1].data['name'][i]=='antlia_2':
            rplot=1.7
        
        Q=str.format('''select source_id,ra,dec,pmra,pmdec,pmra_error,pmdec_error,pmra_pmdec_corr,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,phot_g_mean_flux,phot_g_mean_flux_error,phot_bp_mean_flux,phot_bp_mean_flux_error,phot_rp_mean_flux,phot_rp_mean_flux_error,parallax,parallax_error,astrometric_excess_noise,duplicated_source from gaia_edr3.gaia_source where phot_g_mean_mag=phot_g_mean_mag and q3c_radial_query(ra,dec,{center.ra.deg},{center.dec.deg},{rsearch}) and astrometric_excess_noise<exp(greatest(1.5,1.5+0.3*(phot_g_mean_mag-18))) and parallax<3*parallax_error limit 10000000''', **locals())
        gaia_id,gaia_radeg,gaia_decdeg,gaia_pmra,gaia_pmdec,gaia_sigpmra,gaia_sigpmdec,gaia_pmra_pmdec_corr,gaia_gmag,gaia_bpmag,gaia_rpmag,gaia_gflux,gaia_siggflux,gaia_bpflux,gaia_sigbpflux,gaia_rpflux,gaia_sigrpflux,gaia_parallax,gaia_sigparallax,gaia_astrometric_excess_noise,gaia_duplicated_source0=sqlutil.get(Q,db='wsdb',host='wsdb.hpc1.cs.cmu.edu',user='matt_walker',password='5UN*ur8Y')
        gaia_coords=SkyCoord(gaia_radeg,gaia_decdeg,unit=(u.deg,u.deg))
        gaia_xi,gaia_eta=mycode.etaxiarr(gaia_coords.ra.rad,gaia_coords.dec.rad,center.ra.rad,center.dec.rad)
        gaia_r=np.sqrt(gaia_xi**2+gaia_eta**2)
        gaia_ebv=dustmaps.sfd.SFDQuery().query_equ(gaia_radeg,gaia_decdeg)
        gaia_siggmag=(2.5/gaia_gflux/2.30258)**2*gaia_siggflux**2
        gaia_sigbpmag=(2.5/gaia_bpflux/2.30258)**2*gaia_sigbpflux**2
        gaia_sigrpmag=(2.5/gaia_rpflux/2.30258)**2*gaia_sigrpflux**2
        a0=3.1*gaia_ebv*0.86#use Eq. 1 from Babusiaux et al 2018 (1804.09378) to get extinction coefficients, apply extinction to update BP-RP, then interate again to get final estimate of extinction coefficient.  The factor of 0.86 comes from Sergey, as explained to me in email from Andrew on July 11, 2022: "I've recently realized that Babusieux's formula was written to be applied to "True" EBV. But the SFD EBV is overestimated according to Schlafly11. So if you use SFD EBV's for the getDust, then the ebv needs to be multiplied by 0.86. (when you work with other surveys it's usually not a problem, because the coefficients are already given wrt SFD, but that I learned doesn't apply for gaia's code."
        kg_c=np.array([0.9761,-0.1704,0.0086,0.0011,-0.0438,0.0013,0.0099])
        kbp_c=np.array([1.1517,-0.0871,-0.0333,0.0173,-0.0230,0.0006,0.0043])
        krp_c=np.array([0.6104,-0.0170,-0.0026,-0.0017,-0.0078,0.00005,0.0006])
        col0=gaia_bpmag-gaia_rpmag
        kg=kg_c[0]+kg_c[1]*col0+kg_c[2]*col0**2+kg_c[3]*col0**3+kg_c[4]*a0+kg_c[5]*a0**2+kg_c[6]*col0*a0
        kbp=kbp_c[0]+kbp_c[1]*col0+kbp_c[2]*col0**2+kbp_c[3]*col0**3+kbp_c[4]*a0+kbp_c[5]*a0**2+kbp_c[6]*col0*a0
        krp=krp_c[0]+krp_c[1]*col0+krp_c[2]*col0**2+krp_c[3]*col0**3+krp_c[4]*a0+krp_c[5]*a0**2+krp_c[6]*col0*a0
        ag=kg*a0
        abp=kbp*a0
        arp=krp*a0
        col1=gaia_bpmag-abp-(gaia_rpmag-arp)
        kg=kg_c[0]+kg_c[1]*col1+kg_c[2]*col1**2+kg_c[3]*col1**3+kg_c[4]*a0+kg_c[5]*a0**2+kg_c[6]*col1*a0
        kbp=kbp_c[0]+kbp_c[1]*col1+kbp_c[2]*col1**2+kbp_c[3]*col1**3+kbp_c[4]*a0+kbp_c[5]*a0**2+kbp_c[6]*col1*a0
        krp=krp_c[0]+krp_c[1]*col1+krp_c[2]*col1**2+krp_c[3]*col1**3+krp_c[4]*a0+krp_c[5]*a0**2+krp_c[6]*col1*a0
        gaia_ag=kg*a0
        gaia_abp=kbp*a0
        gaia_arp=krp*a0
        gaia_gmag_dered=gaia_gmag-gaia_ag
        gaia_bpmag_dered=gaia_bpmag-gaia_abp
        gaia_rpmag_dered=gaia_rpmag-gaia_arp

        g0=[]
        bp0=[]
        rp0=[]
        eep0=np.linspace(202,707,400)
        for j in range(0,len(eep0)):
            pars=[eep0[j],np.log10(age),dsph[1].data['feh'][i]]
            shite=mist.interp_mag(pars+[10.**((dsph[1].data['dmodulus'][i]+5.)/5.),0.0],['G','BP','RP'])
            if shite[3][0]==shite[3][0]:
                g0.append(shite[3][0])
                bp0.append(shite[3][1])
                rp0.append(shite[3][2])
        g0=np.array(g0)
        bp0=np.array(bp0)
        rp0=np.array(rp0)

        iso_color=bp0-rp0
        iso_mag=g0
        gaia_mag=gaia_gmag_dered
        gaia_color=gaia_bpmag_dered-gaia_rpmag_dered
        gaia_err=np.sqrt(gaia_siggmag**2+gaia_sigbpmag**2+gaia_sigrpmag**2)

        phot_target=np.zeros(len(gaia_id),dtype='int')
        for j in range(0,len(gaia_id)):
            if gaia_color[j]==gaia_color[j]:
                mindist=np.min(np.sqrt((gaia_color[j]-iso_color)**2+(gaia_mag[j]-iso_mag)**2))
                if mindist<=np.max(np.array([0.15,gaia_err[j]])):
                    phot_target[j]=1

#        def myprior(cube):
#            prior=[]
#            prior.append([-2.,2.])#log of scaling constant that sets member fraction
#            prior.append([-100.,100.])
#            prior.append([-100.,100.])
#            prior.append([-100.,100.])
#            prior.append([-100.,100.])
#            prior.append([-1.,2.])
#            prior.append([-1.,2.])
#            prior=np.array(prior)
#            x=np.array(cube)
#            for i in range(0,len(x)):
#                x[i]=prior[i][0]+(prior[i][1]-prior[i][0])*cube[i]
#            return x

#        def myloglike(cube):
#            n1=1./np.sqrt(2.*np.pi*gaia_sigpmra[phot_target==1]**2)*np.exp(-0.5*(gaia_pmra[phot_target==1]-cube[1])**2/(gaia_sigpmra[phot_target==1]**2))
#            n2=1./np.sqrt(2.*np.pi*gaia_sigpmdec[phot_target==1]**2)*np.exp(-0.5*(gaia_pmdec[phot_target==1]-cube[2])**2/(gaia_sigpmdec[phot_target==1]**2))
#            n3=1./np.sqrt(2.*np.pi*(gaia_sigpmra[phot_target==1]**2+(10.**cube[5])**2))*np.exp(-0.5*(gaia_pmra[phot_target==1]-cube[3])**2/(gaia_sigpmra[phot_target==1]**2+(10.**cube[5])**2))
#            n4=1./np.sqrt(2.*np.pi*(gaia_sigpmdec[phot_target==1]**2+(10.**cube[6])**2))*np.exp(-0.5*(gaia_pmdec[phot_target==1]-cube[4])**2/(gaia_sigpmdec[phot_target==1]**2+(10.**cube[6])**2))
#            logl=np.sum(np.log(((10.**cube[0])*(1.+(gaia_r[phot_target==1]/dsph_rhalf[i])**2)**(-2.)*n1*n2+n3*n4)/(1.+(10.**cube[0])*(1.+(gaia_r[phot_target==1]/dsph_rhalf[i])**2)**(-2))))
#            return logl

#        parameters=['1','2','3','4','5','6','7']
#        n_params=len(parameters)
#        prefix='chains/'+dsph_name[i]
#        result=solve(LogLikelihood=myloglike,Prior=myprior,n_dims=n_params,outputfiles_basename=prefix,verbose=True,resume=True)
#        a=Analyzer(n_params,outputfiles_basename=prefix)
#        bestfit=a.get_best_fit()

#        p1=1./np.sqrt(2.*np.pi*gaia_sigpmra**2)*np.exp(-0.5*(gaia_pmra-bestfit['parameters'][1])**2/(gaia_sigpmra**2))
#        p2=1./np.sqrt(2.*np.pi*gaia_sigpmdec**2)*np.exp(-0.5*(gaia_pmdec-bestfit['parameters'][2])**2/(gaia_sigpmdec**2))
#        p3=1./np.sqrt(2.*np.pi*(gaia_sigpmra**2+(10.**bestfit['parameters'][5])**2))*np.exp(-0.5*(gaia_pmra-bestfit['parameters'][3])**2/(gaia_sigpmra**2+(10.**bestfit['parameters'][5])**2))
#        p4=1./np.sqrt(2.*np.pi*(gaia_sigpmdec**2+(10.**bestfit['parameters'][6])**2))*np.exp(-0.5*(gaia_pmdec-bestfit['parameters'][4])**2/(gaia_sigpmdec**2+(10.**bestfit['parameters'][6])**2))
#        plummer_integral=np.pi*dsph_rhalf[i]**2*rsearch**2/dsph_rhalf[i]**2/(1.+rsearch**2/dsph_rhalf[i]**2)
#        pm=10.**(bestfit['parameters'][0])*plummer_integral/(1.+10.**(bestfit['parameters'][0])*plummer_integral)
#        pn=1./(10.**bestfit['parameters'][0]*(1.+rsearch**2/dsph_rhalf[i]**2)**(-2.)+1.)
#        gaia_pmem=p1*p2*pm/(p1*p2*pm+p3*p4*pn)#membership probability for gaia-queried stars based only on proper motion, regardless of CMD filter or position
#        dpm=np.sqrt((gaia_pmra-dsph[1].data['pmra_dr3'][i])**2+(gaia_pmdec-dsph[1].data['pmdec_dr3'][i]))**2
#        var_dpm=((gaia_pmra-dsph[1].data['pmra_dr3'][i])**2*(gaia_sigpmra**2+dsph[1].data['sigpmra_dr3'][i]**2)+(gaia_pmdec-dsph[1].data['pmdec_dr3'][i])**2*(gaia_sigpmdec**2+dsph[1].data['sigpmdec_dr3'][i]**2))/((gaia_pmra-dsph[1].data['pmra_dr3'][i])**2+(gaia_pmdec-dsph[1].data['pmdec_dr3'][i])**2)
#        sig_dpm=np.sqrt(var_dpm)
#        dpmspec=np.sqrt((catalog['gaia_pmra']-dsph[1].data['pmra_dr3'][i])**2+(catalog['gaia_pmdec']-dsph[1].data['pmdec_dr3'][i]))**2
#        var_dpmspec=((catalog['gaia_pmra']-dsph[1].data['pmra_dr3'][i])**2*(catalog['gaia_sigpmra']**2+dsph[1].data['sigpmra_dr3'][i]**2)+(catalog['gaia_pmdec']-dsph[1].data['pmdec_dr3'][i])**2*(catalog['gaia_sigpmdec']**2+dsph[1].data['sigpmdec_dr3'][i]**2))/((catalog['gaia_pmra']-dsph[1].data['pmra_dr3'][i])**2+(catalog['gaia_pmdec']-dsph[1].data['pmdec_dr3'][i])**2)
#        sig_dpmspec=np.sqrt(var_dpmspec)
#        dv=catalog['vlos']-dsph[1].data['vlos'][i]
#        var_dv=catalog['vlos_error']**2+dsph[1].data['sigvlos'][i]**2+dsph[1].data['vdisp'][i]**2
#        sig_dv=np.sqrt(var_dv)
#        
#        dpm_norm=np.abs(dpm/sig_dpm)
#        dv_norm=np.abs(dv/sig_dv)
#        dpmspec_norm=np.abs(dpmspec/sig_dpmspec)
#        pm_mem=np.where(dpm_norm<=3.)[0]
#        v_mem=np.where(dv_norm<=3.)[0]
#        pmspec_mem=np.where(dpmspec_norm<=3.)[0]
#        vpmspec_mem=np.where((dv_norm<=3.)&(dpmspec_norm<=3.))[0]

        plt.rc('xtick',labelsize='6')
        plt.rc('ytick',labelsize='6')
        gs=plt.GridSpec(22,22)
        gs.update(wspace=0,hspace=0)
        fig=plt.figure(figsize=(6,6))
        ax1=fig.add_subplot(gs[0:4,0:4])
        ax2=fig.add_subplot(gs[0:4,6:10])
        ax3=fig.add_subplot(gs[0:4,12:16])
        ax4=fig.add_subplot(gs[0:4,18:22])

        fake_x=np.array([-999,-999])
        fake_y=np.array([-999,-999])
        fake_logg=np.array([0,5])

        print(len(this),len(this1),len(thismem),len(thismem3),len(thismem4),len(thismem5),' bbbbbbbbbbbbbbb')
        
        order=np.argsort(catalog['logg_mean'][this1])[::-1]
        ax2.scatter(gaia_bpmag_dered[gaia_r<60.]-gaia_rpmag_dered[gaia_r<60.],gaia_gmag_dered[gaia_r<60.],s=0.5,alpha=0.3,color='0.7',label=r'$R<1^{\circ}$',rasterized=True)
#        ax1.scatter(gaia_bpmag_dered[((gaia_r<60.)&(dpm_norm<3.))]-gaia_rpmag_dered[((gaia_r<60.)&(dpm_norm<3.))],gaia_gmag_dered[((gaia_r<60.)&(dpm_norm<3.))],s=1,alpha=0.5,color='k',label=r'$R<1^{\circ}+$PM',rasterized=True)
#        plt.scatter(gaia_bpmag_dered[((phot_target==1))]-gaia_rpmag_dered[((phot_target==1))],gaia_gmag_dered[((phot_target==1))],s=3,alpha=0.3,color='r')
        c2=ax2.scatter(np.concatenate([catalog['gaia_bpmag_dered'][this1][order]-catalog['gaia_rpmag_dered'][this1][order],fake_x]),np.concatenate([catalog['gaia_gmag_dered'][this1][order],fake_y]),s=1,c=np.concatenate([catalog['logg_mean'][this1][order],fake_logg]),alpha=0.5,cmap='jet',rasterized=True)
#        ax1.scatter(catalog['gaia_bpmag_dered'][thisnon]-catalog['gaia_rpmag_dered'][thisnon],catalog['gaia_gmag_dered'][thisnon],s=1,c=catalog['logg_mean'][thisnon],alpha=0.5,cmap='jet',label='observed',rasterized=True)
#        ax1.scatter(catalog['gaia_bpmag_dered'][thismem]-catalog['gaia_rpmag_dered'][thismem],catalog['gaia_gmag_dered'][thismem],s=1,c=catalog['logg_mean'][thismem],alpha=0.5,cmap='jet',label='observed',rasterized=True)
#        clb=plt.colorbar(c1,location='right',ax=ax1)
#        clb.ax.set_title(label=r'$\log g$',fontsize=6)
        ax2.plot(iso_color,iso_mag,color='k',lw=0.5,linestyle='--')
        ax2.set_xlim([-0.5,2])
        ax2.set_ylim(maglim)
        ax2.set_yticks(maglimticks)
        ax2.set_yticklabels(maglimticks,fontsize=6)
        ax2.set_xticks([0,1,2])
        ax2.set_xticklabels([0,1,2],fontsize=6)
        ax2.set_xlabel('BP-RP',fontsize=6)
        ax2.set_ylabel('G',fontsize=6,labelpad=-0.5)
#        if dsph[1].data['name'][i]=='antlia_2':
#            ax1.legend(loc=2,fontsize=6,borderaxespad=0)

        for j in range(0,len(field_center)):
            field_center_xi,field_center_eta=mycode.etaxiarr(field_center[j].ra.rad,field_center[j].dec.rad,center.ra.rad,center.dec.rad)
            axis_major=0.5
            axis_minor=0.5
            ell=Ellipse(xy=[field_center_xi/60.,field_center_eta/60],height=axis_major,width=axis_minor,angle=0,fc='0.7',ec=None,fill=True,alpha=1)
            ax1.add_artist(ell)
        axis_major=2.*2.*dsph[1].data['rhalf'][i]/np.sqrt(1.-dsph[1].data['ellipticity'][i])/60.#twice Rhalf
        axis_minor=axis_major*(1.-dsph[1].data['ellipticity'][i])
        ell=Ellipse(xy=[0,0],height=axis_major,width=axis_minor,angle=-dsph[1].data['PA'][i],fc=None,ec='k',fill=False,linestyle='--')

#        fig,ax=plt.subplots()
        ax1.scatter(gaia_xi[phot_target==1]/60.,gaia_eta[phot_target==1]/60.,s=0.5,alpha=0.3,color='0.5',label='isochrone',rasterized=True)
#        ax2.scatter(gaia_xi[((phot_target==1)&(dpm_norm<3.))]/60.,gaia_eta[((phot_target==1)&(dpm_norm<3.))]/60.,s=1,alpha=0.5,color='k',label='isocrhone+PM',rasterized=True)
        c1=ax1.scatter(np.concatenate([xi[this1][order]/60.,fake_x]),np.concatenate([eta[this1][order]/60.,fake_y]),c=np.concatenate([catalog['logg_mean'][this1][order],fake_logg]),s=1,alpha=0.5,cmap='jet',rasterized=True)
#        clb=plt.colorbar(c2,location='right',ax=ax2)
#        ax2.scatter(xi[thisnon]/60.,eta[thisnon]/60.,c=catalog['logg_mean'][thisnon],s=1,alpha=0.5,cmap='jet',label='observed',rasterized=True)
#        ax2.scatter(xi[thismem]/60.,eta[thismem]/60.,c=catalog['logg_mean'][thismem],s=1,alpha=0.5,cmap='jet',label='observed',rasterized=True)
        ax1.add_artist(ell)
        ax1.set_xlim([rplot,-rplot])
        ax1.set_ylim([-rplot,rplot])
#        ax1.text(0.05,0.94,field_object_formal[this1[0]],horizontalalignment='left',verticalalignment='top',transform=ax1.transAxes,fontsize=6)
        ax1.set_xlabel(r'$\Delta$R.A. [deg]',fontsize=6)
        ax1.set_ylabel(r'$\Delta$Dec. [deg]',labelpad=-2,fontsize=6)
        ax1.text(0.05,1.15,field_object_formal[this1[0]],horizontalalignment='left',verticalalignment='top',transform=ax1.transAxes,fontsize=6)


        xxx=np.concatenate([catalog['gaia_pmra'][this1][order],fake_x])
        yyy=np.concatenate([catalog['gaia_pmdec'][this1][order],fake_y])
        ggg=np.concatenate([catalog['logg_mean'][this1][order],fake_logg])
#        c4=ax4.scatter(np.concatenate([catalog['vlos_mean'][this1][order],fake_x]),np.concatenate([catalog['feh_mean'][this1][order],fake_y]),s=1,c=np.concatenate([catalog['logg_mean'][this1][order],fake_logg]),alpha=0.5,cmap='jet',rasterized=True)
        c3=ax3.scatter(xxx,yyy,s=1,c=ggg,alpha=0.5,cmap='jet',rasterized=True,vmin=0,vmax=5)
#        ax4.scatter(catalog['gaia_bpmag_dered'][thisnon]-catalog['gaia_rpmag_dered'][thisnon],catalog['gaia_gmag_dered'][thisnon],s=1,c=catalog['logg_mean'][thisnon],alpha=0.5,cmap='jet',label='observed',rasterized=True)
#        ax4.scatter(catalog['gaia_bpmag_dered'][thismem]-catalog['gaia_rpmag_dered'][thismem],catalog['gaia_gmag_dered'][thismem],s=1,c=catalog['logg_mean'][thismem],alpha=0.5,cmap='jet',label='observed',rasterized=True)
#        clb.ax.set_yticklabels([0,1,2,3,4,5])
        ax3.set_xlim([-5,5])
        ax3.set_ylim([-5,5])
        axis_major=2.*dsph[1].data['sigpmra_dr3'][i]
        axis_minor=2.*dsph[1].data['sigpmdec_dr3'][i]
        if dsph[1].data['pmra_dr3'][i]==dsph[1].data['pmra_dr3'][i]:
            ax3.plot([dsph[1].data['pmra_dr3'][i],dsph[1].data['pmra_dr3'][i]],[-999,999],color='k',linestyle='--',lw=0.5)
            ax3.plot([-999,999],[dsph[1].data['pmdec_dr3'][i],dsph[1].data['pmdec_dr3'][i]],color='k',linestyle='--',lw=0.5)
#        ax3.add_artist(ell)
#        ax3.set_xticks([-300,-200,-100,0,100,200,300])
#        ax3.set_xticklabels([-300,-200,-100,0,100,200,300],fontsize=6,rotation=90)
#        ax3.set_ylim([-1,1])
#        ax3.set_yticks([-4,-3,-2,-1,0,1])
#        ax3.set_yticklabels([-4,-3,-2,-1,0,1],fontsize=6)

        ax3.set_xlabel(r'$\mu_{\alpha *}$ [mas/year]',fontsize=6)
        ax3.set_ylabel(r'$\mu_{\delta}$ [mas/year]',fontsize=6,labelpad=-2)
#        if dsph[1].data['name'][i]=='antlia_2':
#            ax4.legend(loc=2,fontsize=6,borderaxespad=0)
        
        xxx=np.concatenate([catalog['vlos_mean'][this1][order],fake_x])
        yyy=np.concatenate([catalog['feh_mean'][this1][order],fake_y])+0.27
        ggg=np.concatenate([catalog['logg_mean'][this1][order],fake_logg])
#        c4=ax4.scatter(np.concatenate([catalog['vlos_mean'][this1][order],fake_x]),np.concatenate([catalog['feh_mean'][this1][order],fake_y]),s=1,c=np.concatenate([catalog['logg_mean'][this1][order],fake_logg]),alpha=0.5,cmap='jet',rasterized=True)
        c4=ax4.scatter(xxx,yyy,s=1,c=ggg,alpha=0.5,cmap='jet',rasterized=True,vmin=0,vmax=5)
#        ax4.scatter(catalog['gaia_bpmag_dered'][thisnon]-catalog['gaia_rpmag_dered'][thisnon],catalog['gaia_gmag_dered'][thisnon],s=1,c=catalog['logg_mean'][thisnon],alpha=0.5,cmap='jet',label='observed',rasterized=True)
#        ax4.scatter(catalog['gaia_bpmag_dered'][thismem]-catalog['gaia_rpmag_dered'][thismem],catalog['gaia_gmag_dered'][thismem],s=1,c=catalog['logg_mean'][thismem],alpha=0.5,cmap='jet',label='observed',rasterized=True)
#        clb.ax.set_yticklabels([0,1,2,3,4,5])
        clb=plt.colorbar(c4,location='right',ax=ax4,ticks=[0,1,2,3,4,5])
        clb.ax.set_title(label=r'$\log g$',fontsize=6)
        xmin=dsph[1].data['vlos'][i]-100.
        xmax=dsph[1].data['vlos'][i]+100.
        if(dsph[1].data['vlos'][i]==dsph[1].data['vlos'][i]):
            plt.plot([dsph[1].data['vlos'][i],dsph[1].data['vlos'][i]],[-999,999],color='k',linestyle='--',lw=0.5)
            ax4.set_xlim([xmin,xmax])
        else:
            ax4.set_xlim([-300,300])
#        ax4.set_xticks([-300,-200,-100,0,100,200,300])
#        ax4.set_xticklabels([-300,-200,-100,0,100,200,300],fontsize=6,rotation=90)
        ax4.set_ylim([-4,1])
        ax4.set_yticks([-4,-3,-2,-1,0,1])
        ax4.set_yticklabels([-4,-3,-2,-1,0,1],fontsize=6)
        ax4.set_xlabel(r'$V_{\rm LOS}$ [km/s]',fontsize=6)
        ax4.set_ylabel('[Fe/H]',fontsize=6)
#        if dsph[1].data['name'][i]=='antlia_2':
#            ax4.legend(loc=2,fontsize=6,borderaxespad=0)


        plt.savefig(dsph[1].data['name'][i]+'_m2fs_cmdmap.pdf',dpi=200)
        plt.show()
        plt.close()

        for j in range(0,len(this)):
            ra.append(catalog['ra_deg'][this[j]])
            dec.append(catalog['dec_deg'][this[j]])
            rad.append(r[this[j]])
            v.append(catalog['vlos'][this[j]])
            teff.append(catalog['teff'][this[j]])
            logg.append(catalog['logg'][this[j]])
            z.append(catalog['feh'][this[j]])
            mgfe.append(catalog['mgfe'][this[j]])
            sigv.append(catalog['vlos_error'][this[j]])
            sigteff.append(catalog['teff_error'][this[j]])
            siglogg.append(catalog['logg_error'][this[j]])
            sigfeh.append(catalog['feh_error'][this[j]])
            sigmgfe.append(catalog['mgfe_error'][this[j]])
            teffprior_v.append(catalog['teffprior_vlos'][this[j]])
            teffprior_teff.append(catalog['teffprior_teff'][this[j]])
            teffprior_logg.append(catalog['teffprior_logg'][this[j]])
            teffprior_z.append(catalog['teffprior_feh'][this[j]])
            teffprior_mgfe.append(catalog['teffprior_mgfe'][this[j]])
            teffprior_sigv.append(catalog['teffprior_vlos_error'][this[j]])
            teffprior_sigteff.append(catalog['teffprior_teff_error'][this[j]])
            teffprior_siglogg.append(catalog['teffprior_logg_error'][this[j]])
            teffprior_sigfeh.append(catalog['teffprior_feh_error'][this[j]])
            teffprior_sigmgfe.append(catalog['teffprior_mgfe_error'][this[j]])            
            obj.append(catalog['target_system'][this[j]])
            v_mean.append(catalog['vlos_mean'][this[j]])
            sigv_mean.append(catalog['vlos_mean_error'][this[j]])
            stdv_mean.append(catalog['vlos_mean_scatter'][this[j]])
            teff_mean.append(catalog['teff_mean'][this[j]])
            sigteff_mean.append(catalog['teff_mean_error'][this[j]])
            stdteff_mean.append(catalog['teff_mean_scatter'][this[j]])
            logg_mean.append(catalog['logg_mean'][this[j]])
            siglogg_mean.append(catalog['logg_mean_error'][this[j]])
            stdlogg_mean.append(catalog['logg_mean_scatter'][this[j]])
            z_mean.append(catalog['feh_mean'][this[j]])
            sigfeh_mean.append(catalog['feh_mean_error'][this[j]])
            stdz_mean.append(catalog['feh_mean_scatter'][this[j]])
            mgfe_mean.append(catalog['mgfe_mean'][this[j]])
            sigmgfe_mean.append(catalog['mgfe_mean_error'][this[j]])
            stdmgfe_mean.append(catalog['mgfe_mean_scatter'][this[j]])
            teffprior_v_mean.append(catalog['teffprior_vlos_mean'][this[j]])
            teffprior_sigv_mean.append(catalog['teffprior_vlos_mean_error'][this[j]])
            teffprior_stdv_mean.append(catalog['teffprior_vlos_mean_scatter'][this[j]])
            teffprior_teff_mean.append(catalog['teffprior_teff_mean'][this[j]])
            teffprior_sigteff_mean.append(catalog['teffprior_teff_mean_error'][this[j]])
            teffprior_stdteff_mean.append(catalog['teffprior_teff_mean_scatter'][this[j]])
            teffprior_logg_mean.append(catalog['teffprior_logg_mean'][this[j]])
            teffprior_siglogg_mean.append(catalog['teffprior_logg_mean_error'][this[j]])
            teffprior_stdlogg_mean.append(catalog['teffprior_logg_mean_scatter'][this[j]])
            teffprior_z_mean.append(catalog['teffprior_feh_mean'][this[j]])
            teffprior_sigfeh_mean.append(catalog['teffprior_feh_mean_error'][this[j]])
            teffprior_stdz_mean.append(catalog['teffprior_feh_mean_scatter'][this[j]])
            teffprior_mgfe_mean.append(catalog['teffprior_mgfe_mean'][this[j]])
            teffprior_sigmgfe_mean.append(catalog['teffprior_mgfe_mean_error'][this[j]])
            teffprior_stdmgfe_mean.append(catalog['teffprior_mgfe_mean_scatter'][this[j]])
            obs.append(catalog['obs'][this[j]])
            nobs.append(catalog['n_obs'][this[j]])
            goodobs.append(catalog['good_obs'][this[j]])
            goodnobs.append(catalog['good_n_obs'][this[j]])
            
#            if dv_norm[this[j]]<=3.:
#                vmem.append(1)
#            else:
#                vmem.append(0)
#            if dpmspec_norm[this[j]]<=3.:
#                pmmem.append(1)
#            else:
#                pmmem.append(0)
#            if ((dv_norm[this[j]]<=3.)&(dpmspec_norm[this[j]]<=3)):
#                vpmmem.append(1)
#            else:
#                vpmmem.append(0)
#            if 4==4:
#                isovpmmem.append(1)
#            else:
#                isovpmmem.append(0)

    ra=np.array(ra)
    dec=np.array(dec)
    rad=np.array(rad)
    v=np.array(v)
    teff=np.array(teff)
    logg=np.array(logg)
    z=np.array(z)
    mgfe=np.array(mgfe)
    sigv=np.array(sigv)
    sigteff=np.array(sigteff)
    siglogg=np.array(siglogg)
    sigfeh=np.array(sigfeh)
    sigmgfe=np.array(sigmgfe)
    teffprior_v=np.array(teffprior_v)
    teffprior_teff=np.array(teffprior_teff)
    teffprior_logg=np.array(teffprior_logg)
    teffprior_z=np.array(teffprior_z)
    teffprior_mgfe=np.array(teffprior_mgfe)
    teffprior_sigv=np.array(teffprior_sigv)
    teffprior_sigteff=np.array(teffprior_sigteff)
    teffprior_siglogg=np.array(teffprior_siglogg)
    teffprior_sigfeh=np.array(teffprior_sigfeh)
    teffprior_sigmgfe=np.array(teffprior_sigmgfe)
    isomem=np.array(isomem)
    vmem=np.array(vmem)
    pmmem=np.array(pmmem)
    vpmmem=np.array(vpmmem)
    isovpmmem=np.array(isovpmmem)
    obj=np.array(obj)
    v_mean=np.array(v_mean)
    sigv_mean=np.array(sigv_mean)
    stdv_mean=np.array(stdv_mean)
    teff_mean=np.array(teff_mean)
    sigteff_mean=np.array(sigteff_mean)
    stdteff_mean=np.array(stdteff_mean)
    logg_mean=np.array(logg_mean)
    siglogg_mean=np.array(siglogg_mean)
    stdlogg_mean=np.array(stdlogg_mean)
    z_mean=np.array(z_mean)
    sigfeh_mean=np.array(sigfeh_mean)
    stdz_mean=np.array(stdz_mean)
    mgfe_mean=np.array(mgfe_mean)
    sigmgfe_mean=np.array(sigmgfe_mean)
    stdmgfe_mean=np.array(stdmgfe_mean)
    teffprior_v_mean=np.array(teffprior_v_mean)
    teffprior_sigv_mean=np.array(teffprior_sigv_mean)
    teffprior_stdv_mean=np.array(teffprior_stdv_mean)
    teffprior_teff_mean=np.array(teffprior_teff_mean)
    teffprior_sigteff_mean=np.array(teffprior_sigteff_mean)
    teffprior_stdteff_mean=np.array(teffprior_stdteff_mean)
    teffprior_logg_mean=np.array(teffprior_logg_mean)
    teffprior_siglogg_mean=np.array(teffprior_siglogg_mean)
    teffprior_stdlogg_mean=np.array(teffprior_stdlogg_mean)
    teffprior_z_mean=np.array(teffprior_z_mean)
    teffprior_sigfeh_mean=np.array(teffprior_sigfeh_mean)
    teffprior_stdz_mean=np.array(teffprior_stdz_mean)
    teffprior_mgfe_mean=np.array(teffprior_mgfe_mean)
    teffprior_sigmgfe_mean=np.array(teffprior_sigmgfe_mean)
    teffprior_stdmgfe_mean=np.array(teffprior_stdmgfe_mean)
    obs=np.array(obs)
    nobs=np.array(nobs)
    goodobs=np.array(goodobs)
    goodnobs=np.array(goodnobs)

    pickle.dump((ra,dec,v,sigv,teff,sigteff,logg,siglogg,z,sigfeh,mgfe,sigmgfe,teffprior_v,teffprior_sigv,teffprior_teff,teffprior_sigteff,teffprior_logg,teffprior_siglogg,teffprior_z,teffprior_sigfeh,teffprior_mgfe,teffprior_sigmgfe,vpmmem,obj,v_mean,sigv_mean,stdv_mean,teff_mean,sigteff_mean,stdteff_mean,logg_mean,siglogg_mean,stdlogg_mean,z_mean,sigfeh_mean,stdz_mean,mgfe_mean,sigmgfe_mean,stdmgfe_mean,teffprior_v_mean,teffprior_sigv_mean,teffprior_stdv_mean,teffprior_teff_mean,teffprior_sigteff_mean,teffprior_stdteff_mean,teffprior_logg_mean,teffprior_siglogg_mean,teffprior_stdlogg_mean,teffprior_z_mean,teffprior_sigfeh_mean,teffprior_stdz_mean,teffprior_mgfe_mean,teffprior_sigmgfe_mean,teffprior_stdmgfe_mean,obs,nobs,goodobs,goodnobs),open('m2fs_field_of_halos.pkl','wb'))
    

    
