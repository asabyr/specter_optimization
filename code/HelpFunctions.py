import sys
sys.path.append('../../software/sd_foregrounds_optimize/')
import fisher
import spectral_distortions as sd
import foregrounds_fisher as fgf
from foregrounds_fisher import *
from spectral_distortions import *
import numpy as np
import constants as const

def cost(bands,dets):
    """
    calculates focal plane area for a set of bands and detectors.
    Assumes size of 150GHz is (5.25mm)^2 and others scale as (1/frequency)^2
    Assumes dual-polarization for frequencies < 500GHz 
    (i.e. 2 detectors have an area of just 1 detector)
    
    args:
    bands: frequency bands [Hz]
    dets: detector counts

    output:
    focal plane area [meters^2]
    """
    nu0=150.
    size150=5.25e-3**2. #size of 150GHz detectors
    nu=np.mean(bands, axis=1) #find center frequencies

    dual_pol_inds=np.where(nu<500)
    dets_dual_pol=np.copy(dets)
    dets_dual_pol[dual_pol_inds]=dets_dual_pol[dual_pol_inds]/2.0
    area=np.sum(nu0**2*size150/nu**2*dets_dual_pol)

    return area

def calc_fisher(settings_dict, bands, dets, sys_err=np.array([]),
                arg_dict={}, return_spectrum=False):
    """
    calculates snr/bias/cost for mu or y distortion using sd_foregrounds code
    (the function became bigger than initially intended so a bit long)
    ############
    args:
    settings_dict: dictionary with settings
    ############
    settings_dict['fraction of the sky observed']-- float
    dur_months=settings_dict['mission integration [months]']-- float
    fg=settings_dict['include all foregrounds']-- False or array of FGs to marginalize over
    put_priors=settings_dict['set priors']-- False or dictionary of priors
    prefix=settings_dict['bolocalc file prefix']--str
    sensfile=settings_dict['sensitivity file']--str or False
    distortion=settings_dict['spectral distortion']--'mu' or 'y'
    hemt_amps=settings_dict['hemt amps']-- True/False
    hemt_freq=settings_dict['hemt freq'] -- Frequency division for hemt/bolometers
    binstep=settings_dict['bin step'] -- number of values to use per pand for bandpass integration
    
    ############
    bands: frequency bands [Hz]
    dets: detector counts
    sys_err=np.array([]): systematic error (in this case, calibration error)
    arg_dict={}: sky model dictionary for fisher
    return_spectrum=False: whether to return sky spectrum from Fisher or not
    ############
    return:
    snr, focal plane area [meters^2] (tot_spectrum if asked [Jy/sr]) -- default
    bias [%], bias [\sigma], (tot_spectrum if asked [Jy/sr]) -- if sys_err is set
    
    """
    #most settings, probably better way to organize this
    sky=settings_dict['fraction of the sky observed']
    dur_months=settings_dict['mission integration [months]']
    fg=settings_dict['include all foregrounds']
    put_priors=settings_dict['set priors']
    prefix=settings_dict['bolocalc file prefix']
    sensfile=settings_dict['sensitivity file']
    distortion=settings_dict['spectral distortion']
    hemt_amps=settings_dict['hemt amps']
    hemt_freq=settings_dict['hemt freq']

    if 'bin step' in settings_dict.keys():
        binstep=settings_dict['bin step']
        bandpass=True
    else:
        binstep=0
        bandpass=False

    #if changing fisher sky model
    if len(arg_dict)>0:
        fid_args=arg_dict.copy()
    else:
        fid_args={}

    #use priors or not
    if put_priors!=False:
        fish = fisher.FisherEstimation(fsky=sky, duration=dur_months, bandpass=bandpass, instrument='specter',file_prefix=prefix,freq_bands=bands,Ndet_arr=dets,hemt_amps=hemt_amps, hemt_freq=hemt_freq,noisefile=sensfile, priors=put_priors, systematic_error=sys_err, arg_dict=fid_args, binstep=binstep)
    else:
        fish = fisher.FisherEstimation(fsky=sky, duration=dur_months, bandpass=bandpass, instrument='specter',file_prefix=prefix,freq_bands=bands,Ndet_arr=dets,hemt_amps=hemt_amps, hemt_freq=hemt_freq, noisefile=sensfile, priors={}, systematic_error=sys_err, arg_dict=fid_args, binstep=binstep)

    #option to marginalize over some FGs
    if isinstance(fg,np.ndarray):
        print('setting custom signals')
        fish.set_signals(fncs=fg)
    elif fg==False:
        fish.set_signals(fncs=[sd.DeltaI_mu,sd.DeltaI_DeltaT,sd.DeltaI_reltSZ_2param_yweight])

    #run fisher, get sky spectrum
    fish.run_fisher_calculation()
    tot_spectrum=fish.tot_spectrum

    #if no systematic errors assumed
    if len(sys_err)==0:
        
        area=cost(bands, dets)
        if distortion=="y":
            yfid=1.77e-6
            dist_snr=yfid/fish.errors['y_tot']
            

        elif distortion=="mu":
            mufid=2.e-8
            dist_snr=mufid/fish.errors['mu_amp']
        else:
            sys.exit("choose y or mu")

        #return SNR_i, area and if specified, sky spectrum
        if return_spectrum==True:
            return dist_snr, area, tot_spectrum
        else:
            return dist_snr, area

    #assuming some systematic errors, calculate bias
    elif len(sys_err)==len(dets):
        
        fish.calculate_systematic_bias()

        if distortion=="y":
            yfid=1.77e-6
            bias=fish.B["y_tot"]
            fish_err=fish.errors["y_tot"]
            bias_percent=bias/yfid*100.0
            bias_sigma=bias/fish_err
        elif distortion=="mu":
            mufid=2.e-8
            bias=fish.B["mu_amp"]
            fish_err=fish.errors["mu_amp"]
            bias_percent=bias/mufid*100.0
            bias_sigma=bias/fish_err
        else:
            sys.exit("choose y or mu")

        #return bias in percent diff from fid, bias in sigma & sky spectrum if asked
        if return_spectrum==True:
            return bias_percent, bias_sigma, tot_spectrum
        else:
            return bias_percent, bias_sigma

def convert_to_dict(param_values,param_keys, factor=True):
    """function to make a list of dictionaries,
    that can be an input to fisher calculations. Same as convert_one_dict
    but for many arrays.
    """

    #if an array of arrays is the input
    if len(np.shape(param_values))>1:

        loop_n=len(param_values[:,0])
        dicts=[]
        #loop through each array
        for j in range(loop_n):

            final_dict=convert_one_dict(param_values_one=param_values[j],
                                        param_keys_one=param_keys, factor_one=factor)
            dicts.append(final_dict.copy())

        return dicts
   
    #if only one array is an input
    else:
        final_dict=convert_one_dict(param_values_one=param_values,
                                        param_keys_one=param_keys, factor_one=factor)
        return final_dict
        

def convert_one_dict(param_values_one, param_keys_one, factor_one=True):
    """function to make a dictionary from an array of values,
    that can be an input to fisher calculations.
    
    args:
    param_values [array of floats] -- array of parameter values
    param_keys [array of strings] -- array of parameter names 
    (must match sd_foregrounds_optimize conventions)
    factor [bool] -- True/False, 
    whether the parameter values are specified as some factor of fiducial values.
    
    return:
    final_dict: dictionary that
    includes all sky parameters used in fisher (i.e makes a dictionary in the order that
    fisher.py needs and includes all parameters that are kept at fiducial values)

    """
    
    #amplitude parameters
    amps=np.array(['Ad','As','Acib','EM','Asd','Aco'])

    #convert to actual values if only factors specified
    param_dict={}
    for i in range(len(param_keys_one)):
        if param_keys_one[i] in amps or factor_one==True:
            param_dict[param_keys_one[i]]=param_values_one[i]*const.params_all[param_keys_one[i]]
        else:
            param_dict[param_keys_one[i]]=param_values_one[i]

    #add missing parameters
    for param, param_value in const.params_all.items():
        if param not in param_dict:
            param_dict[param]=param_value

    #put in right order
    final_dict={}
    for param, param_value in const.params_all.items():
        final_dict[param]=param_dict[param]

    return final_dict

def sort_arg_dict(param_dict):

    final_dict={}
    for param, param_value in const.params_all.items():
        final_dict[param]=param_dict[param]

    return final_dict

def return_new_args(factor):

    final_dict={}
    for param, param_value in const.params_all.items():
        final_dict[param]=param_value*factor

    return final_dict
