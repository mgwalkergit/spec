import astropy
import matplotlib
import matplotlib.pyplot as plt
import mycode
import scipy
import numpy as np
matplotlib.use('TkAgg')

nsample=10000
ncalibrate_dark=1000
ncalibrate_bias=1000

err_read=2.6
bias_level=15000.
signal_dark=0.
logsignal_source=np.linspace(-1,4,100)
signal_source=10.**logsignal_source

dark_calibrate=np.random.poisson(lam=signal_dark,size=ncalibrate_dark)+np.random.normal(size=ncalibrate_dark,loc=0.,scale=err_read)+bias_level
bias_calibrate=np.random.normal(size=ncalibrate_bias,loc=0.,scale=err_read)+bias_level

bias_estimate=np.mean(bias_calibrate)
err_bias_estimate=np.std(bias_calibrate)/np.sqrt(ncalibrate_bias)

dark_estimate=np.mean(dark_calibrate-bias_estimate)
err_dark_estimate=np.std(dark_calibrate-bias_estimate)/np.sqrt(ncalibrate_dark)

chi2_obs0=[]
chi2_obs=[]
chi2_obs_fudge=[]
chi2_true=[]
for i in range(0,len(signal_source)):
    print(i)
    read=np.random.normal(size=nsample,loc=0.,scale=err_read)
    source=np.random.poisson(lam=signal_source[i],size=nsample)
    dark=np.random.poisson(lam=signal_dark,size=nsample)
    count=read+bias_level+source+dark 
    count_estimate=count-bias_estimate-dark_estimate

    var_read=err_read**2
    var_bias_estimate=err_bias_estimate**2
    var_dark_estimate=err_dark_estimate**2
    var_source=signal_source[i]
    var_obs0=np.array([np.max([count_estimate[q]+var_read,var_read]) for q in range(0,len(source))])
    var_obs=np.array([np.max([count_estimate[q]+dark_estimate+var_read+var_bias_estimate+var_dark_estimate,var_read+var_bias_estimate+var_dark_estimate]) for q in range(0,len(source))])
    var_obs_fudge=np.array([np.max([count_estimate[q]+dark_estimate+2.+var_read+var_bias_estimate+var_dark_estimate,0.6*(var_read+var_bias_estimate+var_dark_estimate)]) for q in range(0,len(source))])
    var_true=var_source+var_read+signal_dark

    chi2_obs0.append(np.sum((count-bias_estimate-dark_estimate-signal_source[i])**2/var_obs0)/nsample)
    chi2_obs.append(np.sum((count-bias_estimate-dark_estimate-signal_source[i])**2/var_obs)/nsample)
    chi2_obs_fudge.append(np.sum((count-bias_estimate-dark_estimate-signal_source[i])**2/var_obs_fudge)/nsample)
    chi2_true.append(np.sum((count-bias_level-signal_dark-signal_source[i])**2/var_true)/nsample)

chi2_obs0=np.array(chi2_obs0)
chi2_obs=np.array(chi2_obs)
chi2_obs_fudge=np.array(chi2_obs_fudge)
chi2_true=np.array(chi2_true)

plt.plot(signal_source,chi2_true,color='k',label=r'$\mathrm{Var}=\overline{S}+\overline{D}+\sigma^2$')
plt.plot(signal_source,chi2_obs,color='b',label=r'$\mathrm{Var}=\mathrm{max}[C-B+\sigma^2,\sigma^2]$')
plt.plot(signal_source,chi2_obs_fudge,color='g',label=r'$\mathrm{Var}=\mathrm{max}[C-B+2\sigma^2,0.6\sigma^2]$')
plt.plot([0,10000],[1,1],linestyle='--')
plt.xlim([0.1,10000])
plt.ylim([0.5,1.5])
plt.legend()
plt.xscale('log')
plt.show()
plt.close()
