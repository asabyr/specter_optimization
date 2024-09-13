import sys
import configparser
import numpy as np
from HelpFunctions import calc_fisher
import multiprocessing
from functools import partial


#paths
project_dir='/burg/home/as6131/specter_optimization/'
files_dir=project_dir+"files/"

#command line arguments
ini_file=str(sys.argv[1])
run=str(sys.argv[2])

#if split between nodes, can have additional argument
if len(sys.argv)>4:
    print("Too many arguments")
    sys.exit()
elif len(sys.argv)==4:
    n_chunk=int(sys.argv[3])


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

#hemt additions
snr_settings['hemt amps']=config.getboolean('General','hemt amps')

if snr_settings['hemt amps']==True:
    snr_settings['hemt freq']=config.getfloat('General','hemt freq')
else:
    snr_settings['hemt freq']=0

#detector specifications
bands_file_name=config.get('Detector Optimization Run '+run,'optimized bands')
bands_file=files_dir+bands_file_name
det_mins=np.asarray(config.get('Detector Optimization Run '+run,'minimum detectors').split(','))
det_maxs=np.asarray(config.get('Detector Optimization Run '+run,'maximum detectors').split(','))
Ns=np.asarray(config.get('Detector Optimization Run '+run,'linspace n').split(','),dtype=int)
times=np.asarray(config.get('Detector Optimization Run '+run,'grid dimension').split(','),dtype=int)
cores=int(config['Detector Optimization Run '+run]['number of cores'])
split=int(config['Detector Optimization Run '+run]['split jobs'])
fprefix=bands_file_name.replace('_optimized_bands.txt','')
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"

#make a function to run fisher calculations in parallel
def parallel_runs(detectors):
    pool = multiprocessing.Pool(processes=cores)
    initial_bands_file=np.loadtxt(bands_file, delimiter=',')
    initial_bands=initial_bands_file[:,:2].reshape(len(initial_bands_file[:,0]),2)
    SNR=partial(calc_fisher,snr_settings, initial_bands)

    return pool.map(SNR, detectors)

if __name__ == '__main__':

    #split across nodes
    if split > 1:
        print("splitting into multiple files")
        multichroic_list=np.loadtxt(files_dir+fprefix+"_ndet_combinations_run"+run+".txt")
        comb_chunks=np.array_split(multichroic_list, split, axis=0)
        snr=parallel_runs(comb_chunks[n_chunk])
        SNRs_file=open(files_dir+fprefix+"SNRs_run"+str(run)+"_"+str(n_chunk)+".txt","w")
        np.savetxt(SNRs_file,snr)
    
    #single node
    else:
        multichroic_list=np.loadtxt(files_dir+fprefix+"_ndet_combinations_run"+run+".txt")
        snr=parallel_runs(multichroic_list)
        SNRs_file=files_dir+fprefix+"SNRs_run"+str(run)+".txt"
        np.savetxt(SNRs_file,snr)
