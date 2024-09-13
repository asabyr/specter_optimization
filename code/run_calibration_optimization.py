import sys
import configparser
import numpy as np
from HelpFunctions import calc_fisher
import multiprocessing
from functools import partial

#paths
project_dir='/burg/home/as6131/specter_optimization/'
files_dir=project_dir+"files/"
calib_files_dir=project_dir+"calib_files/"

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

#calibration specifications
bands_file_name=config.get('Calibration Optimization Run '+run,'optimized bands')
bands_file=files_dir+bands_file_name
initial_bands_file=np.loadtxt(bands_file, delimiter=',')
initial_bands=initial_bands_file[:,:2].reshape(len(initial_bands_file[:,0]),2)
detectors=np.asarray(config.get('Calibration Optimization Run '+run,'detectors').split(','), dtype=int)
err_type=config.get('Calibration Optimization Run '+run,'error type')
cores=int(config['Calibration Optimization Run '+run]['number of cores'])
split=int(config['Calibration Optimization Run '+run]['split jobs'])
fprefix=bands_file_name.replace('_optimized_bands.txt','')
snr_settings['sensitivity file']=files_dir+fprefix+"_raw_sensitivity.txt"

def mod_func(tup):
    func, kwargs=tup
    return func(**kwargs)

#make a function to run in fisher calculations in parallel
def parallel_runs(calib_errors):
    pool = multiprocessing.Pool(processes=cores)
    func_tuple=[(calc_fisher, {"settings_dict":snr_settings, "bands":initial_bands, "dets":detectors, 'sys_err':val}) for val in calib_errors]
    return pool.map(mod_func, func_tuple)

if __name__ == '__main__':

    #split across nodes
    if split > 1:
        print("splitting into multiple files")
        multichroic_list=np.loadtxt(calib_files_dir+fprefix+"_calib_combinations_run"+run+".txt")
        comb_chunks=np.array_split(multichroic_list, split, axis=0)
        bias=parallel_runs(comb_chunks[n_chunk])
        bias_file=open(calib_files_dir+fprefix+"Bs_run"+str(run)+"_"+str(n_chunk)+".txt","w")
        np.savetxt(bias_file,bias)
    
    #single node
    else:
        multichroic_list=np.loadtxt(calib_files_dir+fprefix+"_calib_combinations_run"+run+".txt")
        bias=parallel_runs(multichroic_list)
        bias_file=calib_files_dir+fprefix+"Bs_run"+str(run)+".txt"
        np.savetxt(bias_file,bias)
