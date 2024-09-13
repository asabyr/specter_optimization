import numpy as np
import pickle
import constants as const
import copy
import sys
sys.path.append('../../software/sd_foregrounds_optimize/')
import fisher 
from foregrounds_fisher import *
from spectral_distortions import *

def factors_to_values(test_keys, factors):
    #convert array of factors into dictionary with values
    final_dict={}
    for param,value in const.params_all.items():
        ind=np.where(test_keys==param)[0]
        if len(ind)==1:
            final_dict[param]=factors[ind][0]*value
        else:
            final_dict[param]=value
    
    return final_dict

def fisher_for_contours(settings_dict,bands, dets,
                       sys_err=[], arg_dict=[]):
    
    #set up settings
    fid_values=copy.deepcopy(const.params_all)

    #unpack dictionary
    sky=settings_dict['fraction of the sky observed']
    dur_months=settings_dict['mission integration [months]']
    fg=settings_dict['include all foregrounds']
    put_priors=settings_dict['set priors']
    sensfile=settings_dict['sensitivity file']
    distortion=settings_dict['spectral distortion']
    

    if len(arg_dict)>0:
        fid_args=arg_dict.copy()
    else:
        fid_args={}
    
    #run fisher
    fish = fisher.FisherEstimation(fsky=sky,
                               duration=dur_months,
                               bandpass=False, instrument='specter',
                               file_prefix=None,
                               freq_bands=bands,
                               Ndet_arr=dets,
                               hemt_amps=True, 
                               hemt_freq=10,
                               noisefile=sensfile,
                               priors={}, systematic_error=[],
                               arg_dict=fid_args)
    fish.run_fisher_calculation()
    tot_spectrum=fish.tot_spectrum
    covariance=fish.cov
    
    #convert to array
    fid_values_arr=[]
    errors_arr=[]
    
    for key, value in fid_values.items():
    
        fid_values_arr.append(fish.argvals[key])
        errors_arr.append(fish.errors[key])
        
    #return np.array(fid_values_arr), np.array(errors_arr), tot_spectrum
    return np.array(fid_values_arr), covariance, tot_spectrum

def save_sky_sigma(target, data, snr_fid, run=2):
    
    
    
    ind_sigma=np.where((data.snr_sort-snr_fid>target-0.15) & (data.snr_sort-snr_fid<target+0.15))[0]
    
    sigma_dict={}
    sigma_dict['params']=data.params_snr_sort[ind_sigma].copy()
    sigma_dict['SNRs']=data.snr_sort[ind_sigma]
    
    print(sigma_dict['SNRs'])
    
    if target>0:
        plus_minus="plus"
    else:
        plus_minus="minus"
    
    with open('../files/FG_20percent_run'+str(run)+"_"+plus_minus+str(np.abs(target))+'sigma_params.pkl', 'wb') as pkl_file:
        pickle.dump(sigma_dict, pkl_file)
    
    return sigma_dict

def read_sky_sigma(fname, nu_fg):
    
    sigma_change=np.load("../files/"+fname, allow_pickle=True)
    n_sigma_change=len(sigma_change['SNRs'])
    fg_all_sigma=np.ones((n_sigma_change, len(nu_fg)))
    
    for i in range(n_sigma_change):
        
        FG_keys=np.array(['Bd' ,'Td' ,'Bcib', 'Tcib', 'alps','w2s'])
        new_param_dict=factors_to_values(FG_keys, sigma_change['params'][i])
        nu, bias_FG=plot_tot_FG(new_param_dict)
        fg_all_sigma[i,:]=bias_FG
        
    return np.unique(fg_all_sigma, axis=0)

def plot_tot_FG(param_dict):

    nu_fg=np.arange(1.*10**9, 2000.*10**9,1.e9)
    Aco=param_dict['Aco']
    Aff=param_dict['EM']
    Td=param_dict['Td']
    Bd=param_dict['Bd']
    Ad=param_dict['Bd']
    Acib=param_dict['Acib']
    As=param_dict['As']
    Bs=param_dict['alps']
    w2s=param_dict['w2s']
    Asd=param_dict['Asd']
    Acib=param_dict['Acib']
    Bcib=param_dict['Bcib']
    Tcib=param_dict['Tcib']
    
    corad=co_rad(nu_fg, Aco)
    synch=jens_synch_rad(nu_fg, As, Bs, w2s)
    free=jens_freefree_rad(nu_fg, Aff)
    dust=thermal_dust_rad(nu_fg, Ad, Bd, Td)
    ame=spinning_dust(nu_fg, Asd)
    cib=cib_rad(nu_fg, Acib, Bcib, Tcib)
    tot_fg=corad+synch+free+dust+ame+cib
    
    return nu_fg, tot_fg