import sys
import configparser
import numpy as np
from DetOptimizationFunctions import Ndet_multichroic, make_det_combinations
from HelpFunctions import calc_fisher
from NoiseFunctions import *

#paths
project_dir='/burg/home/as6131/specter_optimization/'
files_dir=project_dir+"files/"
calib_files_dir=project_dir+"calib_files/"

#arguments
ini_file=str(sys.argv[1])
run=str(sys.argv[2])

#read ini file
config=configparser.ConfigParser()
config.read(project_dir+"ini_files/"+ini_file)

#read ini file
#general fisher settings
snr_settings = {}
snr_settings['bolocalc file prefix']=config['General']['bolocalc file prefix']
snr_settings['spectral distortion']=config['General']['spectral distortion']
snr_settings['set priors']=config.getboolean('General','set priors')
snr_settings['include all foregrounds']=config.getboolean('General','include all foregrounds')
snr_settings['mission integration [months]']=config.getfloat('General','mission integration [months]')
snr_settings['fraction of the sky observed']=config.getfloat('General','fraction of the sky observed')
snr_settings['sensitivity file']=False

#hemt additions
snr_settings['hemt amps']=config.getboolean('General','hemt amps')

if snr_settings['hemt amps']==True:
    snr_settings['hemt freq']=config.getfloat('General','hemt freq')
else:
    snr_settings['hemt freq']=0

#bands & calibration specifications
bands_file_name=config.get('Calibration Optimization Run '+run,'optimized bands')
bands_file=files_dir+bands_file_name
initial_bands_file=np.loadtxt(bands_file, delimiter=',')
initial_bands=initial_bands_file[:,:2].reshape(len(initial_bands_file[:,0]),2)

detectors=np.asarray(config.get('Calibration Optimization Run '+run,'detectors').split(','), dtype=np.float64)
print(detectors)
err_type=config.get('Calibration Optimization Run '+run,'error type')
err_mins=np.asarray(config.get('Calibration Optimization Run '+run,'minimum error').split(','),dtype=np.float64)
err_maxs=np.asarray(config.get('Calibration Optimization Run '+run,'maximum error').split(','),dtype=np.float64)
Ns=np.asarray(config.get('Calibration Optimization Run '+run,'linspace n').split(','),dtype=int)
times=np.asarray(config.get('Calibration Optimization Run '+run,'grid dimension').split(','),dtype=int)
cores=int(config['Calibration Optimization Run '+run]['number of cores'])
nodes=int(config['Calibration Optimization Run '+run]['number of nodes'])
split=int(config['Calibration Optimization Run '+run]['split jobs'])
fprefix=bands_file_name.replace('_optimized_bands.txt','')
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"

#make all possible calib configurations
print("Making & saving calibration combinations")

paramlist=make_det_combinations(err_mins, err_maxs, Ns)
multichroic_list=Ndet_multichroic(paramlist,times.astype(np.int64))
param_file=calib_files_dir+fprefix+"_calib_combinations_run"+run+".txt"

if err_type=='fraction':
    freqs, sens=specter_sens(prefix=snr_settings['bolocalc file prefix'], bands=initial_bands, dets=detectors, hemt_amps = snr_settings['hemt amps'], hemt_freq = snr_settings['hemt freq'], precompute=snr_settings['sensitivity file'], frac_sky=snr_settings['fraction of the sky observed'], exposure_months=snr_settings['mission integration [months]'])
    np.savetxt(files_dir+fprefix+"_run_"+run+"_nominal_sensitivity.txt", np.column_stack([freqs,sens]))
    multichroic_list=multichroic_list*sens
    np.savetxt(param_file,multichroic_list)
elif err_type=='Jysr':
    np.savetxt(param_file,multichroic_list)
elif err_type=='muK_RJ':
    freqs_GHz, noise=np.loadtxt(files_dir+fprefix+"_raw_sensitivity.txt", unpack=True)
    freqs=freqs_GHz*1.e9
    sys_term_Jysr=muKtoJypersr_RJ(multichroic_list, freqs)
    np.savetxt(param_file,sys_term_Jysr)
else:
    sys.exit("invalid error type")

#quick check that maximum detector configuration matches what you've tested previously
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"
multichroic_min=Ndet_multichroic(err_mins,times.astype(np.int64))
print(multichroic_min)
if err_type=='fraction':
    multichroic_min=multichroic_min*sens
elif err_type=='muK_RJ':
    multichroic_min_RJ=np.copy(multichroic_min)
    multichroic_min=muKtoJypersr_RJ(multichroic_min_RJ, freqs)

bias_min=calc_fisher(settings_dict=snr_settings, bands=initial_bands, dets=detectors, sys_err=multichroic_min)
print(multichroic_min)
print("minimum bias:")
print(bias_min)
