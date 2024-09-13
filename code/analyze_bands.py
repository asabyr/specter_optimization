import numpy as np
from HelpFunctions import calc_fisher
from BandOptimizationFunctions import *
from NoiseFunctions import specter_sens, getnoise_raw
import sys

class AnalyzeBands:

    def __init__(self, filename, #optimized bands file
                 snr_settings, #settings for calc_fisher function
                 #initial band details like in the ini file
                 min_freq_1=1, max_freq_1=10, dfreq_1=1, freq_step_1=100,
                min_freq_2=10, max_freq_2=2000, dfreq_2=3, freq_step_2=100, 
                #whether to save instantenous noise
                raw_noise_calc=True,
                #some options for initializing bands differently
                ratio=0, equal=0):

        self.filename=filename
        self.snr_settings=snr_settings
        self.snr_settings['sensitivity file']=None #make sure noise is always computed using bolocalc 
        
        self.min_freq_1=min_freq_1
        self.min_freq_2=min_freq_2

        self.max_freq_1=max_freq_1
        self.max_freq_2=max_freq_2

        self.dfreq_1=dfreq_1
        self.dfreq_2=dfreq_2

        self.freq_step_1=freq_step_1
        self.freq_step_2=freq_step_2

        self.raw_noise_calc=True

        self.ratio=ratio
        self.equal=equal

    def read_optimized_bands(self):

        optimized_data=np.loadtxt(self.filename, delimiter=',')
        self.bands=optimized_data[:,:2].reshape(len(optimized_data[:,0]),2)

        self.dets=optimized_data[:,2]
        file_prefix=self.filename.replace('../files/','').replace('_optimized_bands.txt','')

        snr_settings=self.snr_settings.copy()
        snr_settings['bolocalc file prefix']=file_prefix

        self.snr_cost=calc_fisher(snr_settings, self.bands, self.dets)
        self.freqs, self.noise=specter_sens(file_prefix, self.bands, self.dets, hemt_amps=snr_settings['hemt amps'], hemt_freq=snr_settings['hemt freq'])

        if self.raw_noise_calc==True:
            freqs_raw, noise_raw=getnoise_raw(path="../files/", bands=self.bands, prefix=file_prefix, hemt_amps = snr_settings['hemt amps'], hemt_freq=snr_settings['hemt freq'])

        return self.output()

    def read_thresh_bands(self, thresh):

        file_prefix="snr_"+str(int(thresh*100))+"_percent_"+self.filename.replace('../files/','').replace('_optimized_bands.txt','')
        comb_file=self.filename.replace('_optimized_bands','_combined_bands')
        snr_settings=self.snr_settings.copy()
        snr_settings['bolocalc file prefix']=file_prefix

        if snr_settings['hemt amps']==True:
            
            if self.ratio>0:
                test_comb_1=ratio_initial_bands(self.min_freq_1, self.max_freq_1, self.ratio)
            elif self.equal>0:
                test_comb_1=equal_initial_bands(self.min_freq_1, self.max_freq_1, self.equal)
            else:
                test_comb_1=set_initial_bands(self.min_freq_1, self.max_freq_1, self.dfreq_1, self.freq_step_1)
            
            if self.min_freq_2>0:
                if self.ratio>0:
                    test_comb_2=ratio_initial_bands(self.min_freq_2, self.max_freq_2, self.ratio)
                elif self.equal>0:
                    test_comb_2=equal_initial_bands(self.min_freq_2, self.max_freq_2, self.equal)
                else:
                    test_comb_2=set_initial_bands(self.min_freq_2, self.max_freq_2, self.dfreq_2, self.freq_step_2)
            
                test_comb_2=set_initial_bands(self.min_freq_2, self.max_freq_2, self.dfreq_2, self.freq_step_2)
                test_comb=np.concatenate((test_comb_1, test_comb_2))
                print("considering both hemts and bolometers")
            else:
                test_comb=np.copy(test_comb_1)
                print("considering only hemts")
        else:
            if self.ratio>0:
                test_comb_2=ratio_initial_bands(self.min_freq_2, self.max_freq_2, self.ratio)
            elif self.equal>0:
                test_comb_2=equal_initial_bands(self.min_freq_2, self.max_freq_2, self.equal)
            else:
                test_comb_2=set_initial_bands(self.min_freq_2, self.max_freq_2, self.dfreq_2, self.freq_step_2)
            
            test_comb_2=set_initial_bands(self.min_freq_2, self.max_freq_2, self.dfreq_2, self.freq_step_2)
            test_comb=np.copy(test_comb_2)
            print("considering only bolometers")


        test_ndet=np.ones(len(test_comb))

        snr_initial=calc_fisher(settings_dict=snr_settings, bands=test_comb, dets=test_ndet)

        indices, snrs, areas=np.loadtxt(comb_file, delimiter=',', unpack=True, dtype=np.float64)


        snr_thresh=thresh*snr_initial[0]
        
        ind_thresh_first=np.where(snrs<snr_thresh)[0]

        if len(ind_thresh_first)!=0:
            ind_thresh=ind_thresh_first[0]
        else:
            print("cannot return bands for the specified threshold")
            return 0
        

        for ind in indices[:ind_thresh]:

            new_band, new_ndet=combine_band_right(test_comb,test_ndet,int(ind))
            test_comb=np.copy(new_band)
            test_ndet=np.copy(new_ndet)


        bands_thresh=np.copy(test_comb)
        dets_thresh=np.copy(test_ndet)

        filename_thresh='../files/'+file_prefix+"_optimized_bands.txt"
        np.savetxt(filename_thresh, np.column_stack([bands_thresh, dets_thresh]), fmt='%1.1f,%1.1f,%1.1f')

        #output sens + snr results
        thresh_data=np.loadtxt(filename_thresh, delimiter=',')
        self.bands=thresh_data[:,:2].reshape(len(thresh_data[:,0]),2)
        self.dets=thresh_data[:,2]

        self.snr_cost=calc_fisher(snr_settings, self.bands, self.dets)
        self.freqs, self.noise=specter_sens(file_prefix, self.bands, self.dets, hemt_amps=snr_settings['hemt amps'], hemt_freq=snr_settings['hemt freq'])

        if self.raw_noise_calc==True:
            freqs_raw, noise_raw=getnoise_raw(path="../files/", bands=self.bands, prefix=file_prefix, hemt_amps = snr_settings['hemt amps'], hemt_freq=snr_settings['hemt freq'])

        return self.output()

    def read_peak_bands(self):

        file_prefix="peak_"+self.filename.replace('../files/','').replace('_optimized_bands.txt','')
        comb_file=self.filename.replace('_optimized_bands','_combined_bands')
        snr_settings=self.snr_settings.copy()
        snr_settings['bolocalc file prefix']=file_prefix

        indices, snrs, areas=np.loadtxt(comb_file, delimiter=',', unpack=True, dtype=np.float64)
        

        snr_max=np.amax(snrs)
        snr_max_ind=np.where(snrs==snr_max)[0][-1]+1

        if snr_settings['hemt amps']==True:
            if self.ratio>0:
                test_comb_1=ratio_initial_bands(self.min_freq_1, self.max_freq_1, self.ratio)
            elif self.equal>0:
                test_comb_1=equal_initial_bands(self.min_freq_1, self.max_freq_1, self.equal)
            else:
                test_comb_1=set_initial_bands(self.min_freq_1, self.max_freq_1, self.dfreq_1, self.freq_step_1)
            
            if self.min_freq_2>0:
                if self.ratio>0:
                    test_comb_2=ratio_initial_bands(self.min_freq_2, self.max_freq_2, self.ratio)
                elif self.equal>0:
                    test_comb_2=equal_initial_bands(self.min_freq_2, self.max_freq_2, self.equal)
                else:
                    test_comb_2=set_initial_bands(self.min_freq_2, self.max_freq_2, self.dfreq_2, self.freq_step_2)
                test_comb=np.concatenate((test_comb_1, test_comb_2))
                print("considering both hemts and bolometers")
            else:
                test_comb=np.copy(test_comb_1)
                print("considering only hemts")
        else:
            if self.ratio>0:
                test_comb_2=ratio_initial_bands(self.min_freq_2, self.max_freq_2, self.ratio)
            elif self.equal>0:
                test_comb_2=equal_initial_bands(self.min_freq_2, self.max_freq_2, self.equal)
            else:
                test_comb_2=set_initial_bands(self.min_freq_2, self.max_freq_2, self.dfreq_2, self.freq_step_2)
            test_comb=np.copy(test_comb_2)
            print("considering only bolometers")

        test_ndet=np.ones(len(test_comb))

        print("maximum snr")
        print(snr_max)
        print("maximum snr ind")
        print(snr_max_ind)

        for ind in indices[:snr_max_ind]:
            new_band, new_ndet=combine_band_right(test_comb,test_ndet,int(ind))
            test_comb=np.copy(new_band)
            test_ndet=np.copy(new_ndet)

        bands_max=np.copy(test_comb)
        dets_max=np.copy(test_ndet)

        filename_peak='../files/'+file_prefix+"_optimized_bands.txt"
        np.savetxt(filename_peak, np.column_stack([bands_max, dets_max]), fmt='%1.1f,%1.1f,%1.1f')

        peak_data=np.loadtxt(filename_peak, delimiter=',')
        self.bands=peak_data[:,:2].reshape(len(peak_data[:,0]),2)
        self.dets=peak_data[:,2]

        self.snr_cost=calc_fisher(snr_settings, self.bands, self.dets)
        self.freqs, self.noise=specter_sens(file_prefix, self.bands, self.dets, hemt_amps=snr_settings['hemt amps'], hemt_freq=snr_settings['hemt freq'])

        if self.raw_noise_calc==True:
            freqs_raw, noise_raw=getnoise_raw(path="../files/", bands=self.bands, prefix=file_prefix, hemt_amps = snr_settings['hemt amps'], hemt_freq=snr_settings['hemt freq'])

        return self.output()

    def output(self):

        self.out_dict={}
        self.out_dict['bands']=self.bands
        self.out_dict['freqs']=self.freqs
        self.out_dict['dets']=self.dets
        self.out_dict['noise']=self.noise
        self.out_dict['snr_cost']=self.snr_cost

        return self.out_dict
