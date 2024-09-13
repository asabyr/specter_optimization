import sys
import configparser
import numpy as np
from HelpFunctions import calc_fisher,convert_to_dict
import multiprocessing
from functools import partial

#paths
project_dir='/burg/home/as6131/specter_optimization/'
files_dir=project_dir+"files/"
grid_files_dir=project_dir+"grid_files/"

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

#sky exploration specifications
bands_file_name=config.get('Sky Exploration Run '+run,'optimized bands')
bands_file=files_dir+bands_file_name
initial_bands_file=np.loadtxt(bands_file, delimiter=',')
initial_bands=initial_bands_file[:,:2].reshape(len(initial_bands_file[:,0]),2)

detectors=np.asarray(config.get('Sky Exploration Run '+run,'detectors').split(','), dtype=np.float64)
param_keys=np.asarray(config.get('Sky Exploration Run '+run,'param keys').split(','))
factor_change=config.getboolean('Sky Exploration Run '+run,'factor')
cores=int(config['Sky Exploration Run '+run]['number of cores'])
split=int(config['Sky Exploration Run '+run]['split jobs'])
fprefix=bands_file_name.replace('_optimized_bands.txt','')
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"

fid_snr, fid_cost, fid_spectrum=calc_fisher(settings_dict=snr_settings, bands=initial_bands, dets=detectors, return_spectrum=True)
print(f'fiducial snr:{fid_snr}')
#sys.exit(0)
def mod_func(tup):
    func, kwargs=tup
    return func(**kwargs)

#make a function to run in fisher calculations in parallel
def parallel_runs_snr_bias(priors):
    pool = multiprocessing.Pool(processes=cores)
    prior_dicts=convert_to_dict(priors,param_keys, factor=factor_change)
    func_tuple=[(calc_fisher, {"settings_dict":snr_settings, "bands":initial_bands, "dets":detectors, 'sys_err':{}, 'arg_dict':val, 'return_spectrum':True}) for val in prior_dicts]
    return pool.map(mod_func, func_tuple)

if __name__ == '__main__':

    #split across nodes
    if split > 1:
        print("splitting into multiple files")
        multichroic_list=np.loadtxt(grid_files_dir+fprefix+"_sky_combinations_run"+run+".txt")
        comb_chunks=np.array_split(multichroic_list, split, axis=0)
        snr_bias_values=parallel_runs_snr_bias(comb_chunks[n_chunk])
        
        snrs=[fish_values[0] for fish_values in snr_bias_values]
        #print('snrs:')
        #print(snrs)
        
        snr_file=open(grid_files_dir+fprefix+"sky_SNR_bias_run"+str(run)+"_"+str(n_chunk)+".txt","w")
        np.savetxt(snr_file,snrs)
        #print('snrs-fid_snr:')
        #print(snrs-fid_snr)

        #biased_Inu=np.array([fish_values[2] for fish_values in snr_bias_values])
        #bias_Inu=biased_Inu-fid_spectrum
        #bias_spectrum_file=open(grid_files_dir+fprefix+"sky_Inu_bias_run"+str(run)+"_"+str(n_chunk)+".txt","w")
        #np.savetxt(bias_spectrum_file,bias_Inu)

    #single node
    else:
        multichroic_list=np.loadtxt(grid_files_dir+fprefix+"_sky_combinations_run"+run+".txt")
        snr_bias_values=parallel_runs_snr_bias(multichroic_list)
        snrs=[fish_values[0] for fish_values in snr_bias_values]
        snr_file=grid_files_dir+fprefix+"sky_SNR_bias_run"+str(run)+".txt"
        np.savetxt(snr_file,snrs)
        
        #bias_Inu=[fish_values[2] for fish_values in snr_bias_values]
        #bias_spectrum_file=open(grid_files_dir+fprefix+"sky_Inu_bias_run"+str(run)+".txt","w")
        #np.savetxt(bias_spectrum_file,bias_Inu)
