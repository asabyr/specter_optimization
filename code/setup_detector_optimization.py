import sys
import configparser
import numpy as np
from DetOptimizationFunctions import Ndet_multichroic, make_det_combinations
from HelpFunctions import calc_fisher
from NoiseFunctions import getnoise_raw

#paths
project_dir='/burg/home/as6131/specter_optimization/'
files_dir=project_dir+"files/"

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

#bands & detectors specifications
bands_file=config['Detector Optimization Run '+run]['optimized bands']
det_mins=np.asarray(config.get('Detector Optimization Run '+run,'minimum detectors').split(','),dtype=np.float64)
det_maxs=np.asarray(config.get('Detector Optimization Run '+run,'maximum detectors').split(','),dtype=np.float64)
Ns=np.asarray(config.get('Detector Optimization Run '+run,'linspace n').split(','),dtype=np.int64)
times=np.asarray(config.get('Detector Optimization Run '+run,'grid dimension').split(','),dtype=np.int64)
fprefix=bands_file.replace('_optimized_bands.txt','')
#compute raw noise for fisher calculations
print("Precomputing raw noise")
initial_bands_file=np.loadtxt(project_dir+"files/"+bands_file, delimiter=',')
initial_bands=initial_bands_file[:,:2].reshape(len(initial_bands_file[:,0]),2)

specter_freqs, specter_noise=getnoise_raw(path=files_dir, bands=initial_bands, prefix=fprefix, hemt_amps = snr_settings['hemt amps'], hemt_freq = snr_settings['hemt freq'])
print("initial_bands")
print(initial_bands)
#make all possible detector configurations
print("Making & saving detector combinations")
paramlist=make_det_combinations(det_mins, det_maxs, Ns)
multichroic_list=Ndet_multichroic(paramlist,times.astype(np.int64))
param_file=files_dir+fprefix+"_ndet_combinations_run"+run+".txt"
np.savetxt(param_file,multichroic_list)

#quick check that maximum detector configuration matches what you've tested previously
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"
print("sensitivity file")
print(files_dir+fprefix+"_raw_sensitivity.txt")
multichroic_max=Ndet_multichroic(det_maxs,times.astype(np.int64))
SNRmax=calc_fisher(settings_dict=snr_settings, bands=initial_bands, dets=multichroic_max)
print(det_maxs)
print(snr_settings)
print(multichroic_max)
print("maximum SNR and cost:")
print(SNRmax)
