import numpy as np
import os
file_path=os.path.dirname(os.path.abspath(__file__))
root_path=file_path.replace('/specter_optimization/code','/software/bolocalc-space/')
sd_fg_path=file_path.replace('/specter_optimization/code','/software/sd_foregrounds_optimize')
import sys
sys.path.append(root_path+'analyze-bc')
from gen_bolos import GenBolos
import pandas as pd

bc_fp = root_path+'calcBolos.py'
exp_fp = root_path+'Experiments/specter_v1/'
ndp=np.float64
Tcmb=2.7255 #SI
kb=1.380649e-23 #SI
hplanck=6.62607015e-34 #SI
clight=299792458. #SI

def muKtoJypersr(muK,f):
    """
    convert from muK units to Jy/sr based on https://arxiv.org/pdf/1303.5070.pdf & https://arxiv.org/pdf/2010.16405.pdf
    i.e. by taking a derivative of the Planck function

    args:
    muK: sensitivity [muK]
    f: frequency [Hz]

    output:
    sensitivity [Jy/sr]
    """
    x=hplanck*f/(kb*Tcmb)
    factor=2.*hplanck**2/(clight**2*kb*(Tcmb*1.e6)*Tcmb)*1.e26
    Jypersr_muK=factor*(f)**4*np.exp(x)/(np.exp(x)-1.)**2

    return muK*Jypersr_muK

def muKtoWm2perHzsr(muK,f):
    """
    convert from muK units to W*m^2/(Hz*sr)

    args:
    muK: sensitivity [muK]
    f: frequency [Hz]

    output:
    sensitivity [W*m^2/(Hz*sr)]
    """
    Jypersr=muKtoJypersr(muK,f)
    return Jypersr*1.e-26

def JypersrtomuK(Jysr,f):
    """
    convert from Jy/sr units to muK

    args:
    Jysr: sensitivity [Jy/sr]
    f: frequency [Hz]

    output:
    sensitivity [muK]
    """
    Jysr_onemuK=muKtoJypersr(1.,f)
    return Jysr*1./Jysr_onemuK

def muKtoJypersr_RJ(muK,f):
    """
    convert from muK units to Jy/sr in the RJ limit
    """
    dB_dT=2*kb*f**2/clight**2

    return muK*dB_dT*1.e26/1.e6

def JypersrtomuK_RJ(Jysr,f):

    Jysr_onemuK=muKtoJypersr_RJ(1.,f)
    return Jysr*1./Jysr_onemuK

def getnoise_raw(path, prefix,bands, hemt_amps = True, hemt_freq = 10):
    """
    calculate instantenous sensitivity for a set of frequency bands using bolocalc calculator
    (saves the sensitivity in path+prefix+"_raw_sensitivity.txt")

    args:
    path: directory path to the project folder
    prefix: file prefix for bolocalc calculator
    bands: frequency bands [GHz]
    hemt_amps: True/False
    hemt_freq: below this frequency HEMTs will be used if hemt_amps set to True

    output:
    frequencies [GHz]
    sensitivity [muK-rtSec]

    """

    bolos = GenBolos(bc_fp=bc_fp, exp_fp=exp_fp, band_edges=bands, file_prefix=prefix, hemt_amps=hemt_amps, hemt_freq=hemt_freq)
    sensitivity_dict = bolos.calc_bolos()
    noise_df = pd.DataFrame(sensitivity_dict).T

    freqs=np.array(noise_df['Center Frequency'], dtype=ndp) #GHz
    raw_noise= np.array(noise_df['Detector NET_CMB'], dtype=ndp)#muK-rts


    raw_sens_file=path+prefix+"_raw_sensitivity.txt"
    np.savetxt(raw_sens_file,np.column_stack([freqs,raw_noise]))


    return freqs, raw_noise

def getnoise_nominal(prefix, bands, dets, hemt_amps = True, hemt_freq = 10, precompute=False):
    """
    calculate nominal specter sensitivity for a set of frequency bands using bolocalc calculator
    (assumes 6 months of integration, 100% of the sky)

    args:
    prefix: file prefix for bolocalc calculator
    bands: frequency bands [GHz]
    dets: detector counts
    hemt_amps: True/False
    hemt_freq: below this frequency HEMTs will be used if hemt_amps set to True
    precompute: sensitivity file (if precalculated, otherwise False)

    output:
    frequencies [Hz]
    sensitivity [Jy/sr*deg]

    """
    nominal_exposure=6. #months
    nominal_frac_sky=1. #fraction

    exposure_sec=365.25*24.*3600.*nominal_exposure/12. #seconds
    sky_area_deg=4.*np.pi*(180./np.pi)**2*nominal_frac_sky #deg^2


    if precompute:

        noise=np.loadtxt(precompute, dtype=ndp)
        freqs=noise[:,0]*1.e9 #Hz
        noise_muK_rtsec=noise[:,1] #muK-rtSec

    else:

        bolos = GenBolos(bc_fp=bc_fp, exp_fp=exp_fp, band_edges=bands, file_prefix=prefix, hemt_amps=hemt_amps, hemt_freq=hemt_freq)
        sensitivity_dict = bolos.calc_bolos()
        noise_df = pd.DataFrame(sensitivity_dict).T
        freqs=np.array(noise_df['Center Frequency'], dtype=ndp)*1.e9 #Hz
        noise_muK_rtsec= np.array(noise_df['Detector NET_CMB'], dtype=ndp)#muK-rtSec


    noise_muK_deg=noise_muK_rtsec/np.sqrt(exposure_sec)* np.sqrt(sky_area_deg) #muK*deg
    noise_Jysr_deg = muKtoJypersr(noise_muK_deg,freqs) #Jypersr*deg
    noise=noise_Jysr_deg/np.sqrt(dets)

    return freqs, noise


def specter_sens(prefix, bands, dets, hemt_amps = True, hemt_freq = 10, precompute=False, frac_sky=1., exposure_months=6.):
    """
    calculate specter sensitivity

    args:
    prefix: file prefix for bolocalc calculator
    bands: frequency bands [GHz]
    dets: detector counts
    hemt_amps: True/False
    hemt_freq: below this frequency HEMTs will be used if hemt_amps set to True
    precompute: sensitivity file (if precalculated, otherwise False)
    frac_sky: fraction of the sky
    exposure_months: integration time [months]

    output:
    frequencies [Hz]
    sensitivity [Jy/sr]

    """
    freqs, noise = getnoise_nominal(bands=bands, prefix=prefix, dets=dets, hemt_amps = hemt_amps, hemt_freq = hemt_freq, precompute=precompute)
    skysr = 4. * np.pi * (180. / np.pi) ** 2 *frac_sky
    scaled_sens=noise/ np.sqrt(skysr) * np.sqrt(6./exposure_months)

    return freqs, scaled_sens

def pixie_sens(nu, fsky=1.,duration=6., mult=1.):
    """
    calculate pixie sensitivity 
    (taken from fisher.py from https://github.com/mabitbol/sd_foregrounds/blob/master/fisher.py)

    args:
    nu: central frequencies [Hz]
    fsky: fraction of the sky
    duration: integration time [months]
    mult: overall multiplicative factor

    output:
    frequencies [Hz]
    sensitivity [Jy/sr]

    """
    from scipy import interpolate
    sdata=np.loadtxt(sd_fg_path+'/templates/Sensitivities.dat',dtype=ndp)
    fs = sdata[:, 0] * 1.e9
    sens = sdata[:, 1]
    template = interpolate.interp1d(np.log10(fs), np.log10(sens), bounds_error=False, fill_value="extrapolate")
    skysr = 4. * np.pi * (180. / np.pi) ** 2 * fsky

    return 10. ** template(np.log10(nu)) / np.sqrt(skysr) * np.sqrt(15. / duration) * mult * 1.e26


def getnoise_map_domain_muK_arcmin(prefix, bands, dets, hemt_amps = True, hemt_freq = 10, precompute=False, 
                        exposure_months=12.):
    """
    calculate map domain noise in muK-arcmin.

    args:
    prefix: file prefix for bolocalc calculator
    bands: frequency bands [GHz]
    dets: detector counts
    hemt_amps: True/False
    hemt_freq: below this frequency HEMTs will be used if hemt_amps set to True
    precompute: sensitivity file (if precalculated, otherwise False)
    exposure_months: integration time [months]

    output:
    frequencies [Hz]
    sensitivity [muK*arcmin]

    """
    
    
    exposure_sec=365.25*24.*3600.*exposure_months/12. #seconds
    sky_area_arcmin=4.*np.pi*(180./np.pi)**2*60.**2 #arcmin^2


    if precompute:

        noise=np.loadtxt(precompute, dtype=ndp)
        freqs=noise[:,0]*1.e9 #Hz
        noise_muK_rtsec=noise[:,1] #muK-rtSec

    else:

        bolos = GenBolos(bc_fp=bc_fp, exp_fp=exp_fp, band_edges=bands, file_prefix=prefix, hemt_amps=hemt_amps, hemt_freq=hemt_freq)
        sensitivity_dict = bolos.calc_bolos()
        noise_df = pd.DataFrame(sensitivity_dict).T
        freqs=np.array(noise_df['Center Frequency'], dtype=ndp)*1.e9 #Hz
        noise_muK_rtsec= np.array(noise_df['Detector NET_CMB'], dtype=ndp)#muK-rtSec


    noise_muK_arcmin=noise_muK_rtsec/np.sqrt(exposure_sec)* np.sqrt(sky_area_arcmin) #muK*arcmin
    noise=noise_muK_arcmin/np.sqrt(dets)

    return freqs, noise