import sys
import configparser
import numpy as np
from BandOptimizationFunctions import ratio_initial_bands,equal_initial_bands,set_initial_bands, comb_optimize
from HelpFunctions import calc_fisher

project_dir='/burg/home/as6131/specter_optimization/'
bolocalc_sens_dir='/burg/home/as6131/software/bolocalc-space/Experiments/specter_v1/'
files_dir=project_dir+"files/"

ini_file=str(sys.argv[1])

config=configparser.ConfigParser()
config.read(project_dir+"ini_files/"+ini_file)

snr_settings = {}

snr_settings['bolocalc file prefix']=config['General']['bolocalc file prefix']
snr_settings['spectral distortion']=config['General']['spectral distortion']
snr_settings['set priors']=config.getboolean('General','set priors')
snr_settings['include all foregrounds']=config.getboolean('General','include all foregrounds')
snr_settings['mission integration [months]']=config.getfloat('General','mission integration [months]')
snr_settings['fraction of the sky observed']=config.getfloat('General','fraction of the sky observed')
snr_settings['sensitivity file']=False
snr_settings['hemt amps']=config.getboolean('General','hemt amps')

if snr_settings['hemt amps']==True:
    snr_settings['hemt freq']=config.getfloat('General','hemt freq')
else:
    snr_settings['hemt freq']=0

min_freq_bolo=config.getfloat('Band Optimization','lowest frequency edge (bolo)')
max_freq_bolo=config.getfloat('Band Optimization','highest frequency edge (bolo)')
dfreq_bolo=config.getfloat('Band Optimization','lowest frequency band width (bolo)')
freq_step_bolo=config.getfloat('Band Optimization','frequency step increase (bolo)')

if snr_settings['hemt amps']==True:
    min_freq_hemt=config.getfloat('Band Optimization','lowest frequency edge (hemt)')
    max_freq_hemt=config.getfloat('Band Optimization','highest frequency edge (hemt)')
    dfreq_hemt=config.getfloat('Band Optimization','lowest frequency band width (hemt)')
    freq_step_hemt=config.getfloat('Band Optimization','frequency step increase (hemt)')

    if max_freq_hemt==snr_settings['hemt freq']:
        print("hemt frequency settings are consistent")
    else:
        print("hemt frequency settings are NOT consistent")
        sys.exit()

old_start=config.getboolean('Band Optimization','resume optimization')
snr_thresh=config.getfloat('Band Optimization','snr threshold [fractional]')

if config.has_option('Band Optimization', 'specify initial bands & dets')==True:
    print("using specified initial bands and detectors")
    specify_initial=config['Band Optimization']['specify initial bands & dets']

    start_data=np.loadtxt(files_dir+specify_initial, delimiter=',')
    start_bands=start_data[:,:2].reshape(len(start_data[:,0]),2)
    start_dets=start_data[:,2]

else:
    if snr_settings['hemt amps']==True and min_freq_bolo!=0.0:
        if config.has_option('Band Optimization', 'dfreq')==True:
            dfreq=config.getfloat('Band Optimization', 'dfreq')
            start_bands_bolo=equal_initial_bands(min_freq_bolo,max_freq_bolo,dfreq)
            start_bands_hemts=equal_initial_bands(min_freq_hemt,max_freq_hemt,dfreq)
        elif config.has_option('Band Optimization','ratio')==True:
            ratio=config.getfloat('Band Optimization', 'ratio')
            start_bands_bolo=ratio_initial_bands(min_freq_bolo,max_freq_bolo,ratio)
            start_bands_hemts=ratio_initial_bands(min_freq_hemt,max_freq_hemt,ratio)
        else:
            start_bands_bolo=set_initial_bands(low_freq_edge=min_freq_bolo,high_freq_edge=max_freq_bolo,initial_dfreq=dfreq_bolo,interval=freq_step_bolo)
            start_bands_hemts=set_initial_bands(low_freq_edge=min_freq_hemt,high_freq_edge=max_freq_hemt,initial_dfreq=dfreq_hemt,interval=freq_step_hemt)
        start_bands=np.concatenate((start_bands_hemts, start_bands_bolo))

    elif snr_settings['hemt amps']==True and min_freq_bolo==0.0:
        if config.has_option('Band Optimization', 'dfreq')==True:
            dfreq=config.getfloat('Band Optimization', 'dfreq')
            start_bands=equal_initial_bands(min_freq_hemt,max_freq_hemt,dfreq)
        elif config.has_option('Band Optimization','ratio')==True:
            ratio=config.getfloat('Band Optimization', 'ratio')
            start_bands=ratio_initial_bands(min_freq_hemt,max_freq_hemt,ratio)
        else:
            start_bands=set_initial_bands(low_freq_edge=min_freq_hemt,high_freq_edge=max_freq_hemt,initial_dfreq=dfreq_hemt,interval=freq_step_hemt)
    elif snr_settings['hemt amps']==False and min_freq_bolo>0.0:

        if config.has_option('Band Optimization', 'dfreq')==True:
            dfreq=config.getfloat('Band Optimization', 'dfreq')
            start_bands=equal_initial_bands(min_freq_bolo,max_freq_bolo,dfreq)
        elif config.has_option('Band Optimization','ratio')==True:
            ratio=config.getfloat('Band Optimization', 'ratio')
            start_bands=ratio_initial_bands(min_freq_bolo,max_freq_bolo,ratio)
        else:
            start_bands=set_initial_bands(low_freq_edge=min_freq_bolo,high_freq_edge=max_freq_bolo,initial_dfreq=dfreq_bolo,interval=freq_step_bolo)
    else:
        print("check initial band specification")
        sys.exit(0)
    
    start_dets=np.ones(len(start_bands))

print("started combining")
opt_bands, opt_dets=comb_optimize(bolocalc_sensfile_path=bolocalc_sens_dir,output_path=files_dir,initial_bands=start_bands, initial_dets=start_dets, snr_thresh_frac=snr_thresh, gen_snr_settings=snr_settings, restart=old_start)
