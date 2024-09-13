import sys
import configparser
import numpy as np
from DetOptimizationFunctions import Ndet_multichroic, make_det_combinations
from HelpFunctions import calc_fisher,convert_to_dict
from NoiseFunctions import *

#paths
project_dir='/burg/home/as6131/specter_optimization/'
files_dir=project_dir+"files/"
grid_files_dir=project_dir+"grid_files/"

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
bands_file_name=config.get('Sky Exploration Run '+run,'optimized bands')
bands_file=files_dir+bands_file_name
initial_bands_file=np.loadtxt(bands_file, delimiter=',')
initial_bands=initial_bands_file[:,:2].reshape(len(initial_bands_file[:,0]),2)

detectors=np.asarray(config.get('Sky Exploration Run '+run,'detectors').split(','), dtype=np.float64)
param_mins=np.asarray(config.get('Sky Exploration Run '+run,'minimum prior').split(','),dtype=np.float64)
param_maxs=np.asarray(config.get('Sky Exploration Run '+run,'maximum prior').split(','),dtype=np.float64)
Ns=np.asarray(config.get('Sky Exploration Run '+run,'linspace n').split(','),dtype=int)
times=np.asarray(config.get('Sky Exploration Run '+run,'grid dimension').split(','),dtype=int)

fid_params=np.array(['mu_amp', 'y_tot', 'kT_yweight', 'DeltaT_amp', 'Ad', 'Bd', 'Td', 'Acib', 'Bcib', 'Tcib', 'EM', 'As', 'alps', 'w2s', 'Asd', 'Aco'])
param_keys=np.asarray(config.get('Sky Exploration Run '+run,'param keys').split(','))

if set(param_keys).issubset(set(fid_params))==False:
    sys.exit("unknown parameter keys")

factor_change=config.getboolean('Sky Exploration Run '+run,'factor')
cores=int(config['Sky Exploration Run '+run]['number of cores'])
split=int(config['Sky Exploration Run '+run]['split jobs'])
fprefix=bands_file_name.replace('_optimized_bands.txt','')
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"

#make all possible calib configurations
print("Making & saving sky exploration combinations")
paramlist=make_det_combinations(param_mins, param_maxs, Ns)
multichroic_list=Ndet_multichroic(paramlist,times.astype(np.int64))
param_file=grid_files_dir+fprefix+"_sky_combinations_run"+run+".txt"
np.savetxt(param_file,multichroic_list)

#quick check that maximum detector configuration matches what you've tested previously
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"

multichroic_max=Ndet_multichroic(param_maxs,times.astype(np.int64))
max_arg_dict=convert_to_dict(multichroic_max, param_keys,factor=factor_change)

snr, cost=calc_fisher(settings_dict=snr_settings, bands=initial_bands, dets=detectors, arg_dict=max_arg_dict)
print("maximum parameter values lead to snr:")
print(snr)
snr, cost=calc_fisher(settings_dict=snr_settings, bands=initial_bands, dets=detectors)
print("fid snr:")
print(snr)
