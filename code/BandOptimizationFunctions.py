from HelpFunctions import calc_fisher
import numpy as np
import os
import time
import sys

def set_initial_bands(low_freq_edge,high_freq_edge,initial_dfreq, interval=100.):
    """
    create initial set of narrow bands

    args:
    low_freq_edge: lowest frequency edge (float)
    high_freq_edge: highest frequency edge (float)
    initial_dfreq: narrowest band width (float)
    interval: specify interval after which to increase band width by +initial+dfreq (float)

    output:
    frequency bands (np.array)
    """

    bands=[[low_freq_edge, low_freq_edge+initial_dfreq]] #initiallize first band
    i=0
    edge=low_freq_edge #initialize lowest edge

    while edge < high_freq_edge:

        dnu=initial_dfreq+np.floor(bands[i][1]/interval)*initial_dfreq
        bands.append([bands[i][1],bands[i][1]+dnu])
        i+=1
        edge=bands[-1][1]

    bands[-1][1]=high_freq_edge #specify higest edge

    return bands

def equal_initial_bands(min_freq, max_freq, dfreq):

    bands=[[min_freq, min_freq+dfreq]] #initiallize first band
    edge=bands[0][1] #initialize lowest edge
    i=0
    while edge < max_freq:

        bands.append([bands[i][1],bands[i][1]+dfreq])
        i+=1
        edge=bands[-1][1]
    return bands

def ratio_initial_bands(min_freq, max_freq, ratio=0.3):

    dfreq=min_freq*ratio
    bands=[[min_freq, min_freq+dfreq]]
    edge=bands[0][1]
    i=0
    while edge < max_freq:
        dfreq=bands[i][1]*ratio
        bands.append([bands[i][1],bands[i][1]+dfreq])
        i+=1
        edge=bands[-1][1]

    if edge>max_freq:
        bands[-2][1]=max_freq

    return bands[:-1]

def combine_band_right(bands,dets,ind):
    """
    combine neighbboring bands

    args:
    bands: frequency bands (np.array)
    dets: detector counts (np.array)
    ind: index of the lower frequency band (int) *(detector count is dropped for ind+1)

    output:
    frequency bands (np.array)
    detector counts (np.array)
    """

    comb_bands=np.copy(bands)
    comb_bands[ind]=[bands[ind][0], bands[ind+1][1]]
    newbands=np.delete(comb_bands, ind+1, axis=0)
    newndet=np.delete(dets, ind+1, axis=0)

    return newbands, newndet

def comb_optimize(bolocalc_sensfile_path, output_path, initial_bands, initial_dets, snr_thresh_frac, gen_snr_settings, restart):
    """
    find optimal frequency bands by combining neighboring bands that have the least effect on snr
    saves optimal bands to output_path+fprefix+"_optimized_bands.txt"
    saves each step of the optimization in output_path+fprefix+"_combined_bands.txt"

    args:
    bolocalc_sensfile_path: path to sens.npy in bolocalc directory
    output_path: path for output files
    initial_bands: initial frequency bands [Hz]
    initial_dets: initial detector counts
    snr_thresh_frac: snr threshold after which to stop combining bands [fraction of initial snr]
    gen_snr_settings: general settings for fisher calculations (see calc_fisher in HelpFunctions)
    restart: to start new optimization or to continue previously started optimization [False/True]

    output:
    optimal bands
    optimal detectors

    """

    prefix=gen_snr_settings['bolocalc file prefix']


    #delete any existing cached sensitivities
    if os.path.exists(bolocalc_sensfile_path+prefix+"sens_out.npy")==True:
        os.remove(bolocalc_sensfile_path+prefix+"sens_out.npy")

    #initial bands/dets
    test_comb=np.copy(initial_bands)
    test_ndet=np.copy(initial_dets)


    freq=[] #list of dropped bands
    SNRs=[] #list of corresponding best SNRs
    costs=[] #list of corresponding cost

    #initial SNR (+time it)
    start=time.time()
    snr_initial=calc_fisher(settings_dict=gen_snr_settings, bands=test_comb, dets=test_ndet)
    end=time.time()
    print("time initial:")
    print(end-start)
    print("initial SNR:")
    print(snr_initial[0])

    #set snr threshold
    snr_thresh=snr_thresh_frac*snr_initial[0]

    if restart==True:
        optimized_data=np.loadtxt(output_path+prefix+"_optimized_bands.txt", delimiter=',')
        initial_bands=optimized_data[:,:2].reshape(len(optimized_data[:,0]),2)
        initial_dets=optimized_data[:,2]
        snr_initial=calc_fisher(settings_dict=gen_snr_settings, bands=initial_bands, dets=initial_dets)
        print("restarted from snr:")
        print(snr_initial)
        test_comb=np.copy(initial_bands)
        test_ndet=np.copy(initial_dets)
    else:
        #open file for appending at each iteration
        np.savetxt(output_path+prefix+"_combined_bands.txt",['#ind, snr, cost'], fmt='%s')

    #don't combine if already below threhold
    if snr_initial[0] <= snr_thresh:
        print("didn't combine any bands")
        optimized_bands=np.copy(test_comb)
        optimized_dets=np.copy(test_ndet)
        np.savetxt(output_path+prefix+"_optimized_bands.txt",np.column_stack([optimized_bands, optimized_dets]), fmt='%1.1f,%1.1f,%1.1f')
        return optimized_bands, optimized_dets

    print("starting with the following bands:")
    print(test_comb)

    j=0 #while loop counter
    snr_best=snr_initial[0]
    count_NaNs=0.
    count_all=0.

    with open(output_path+prefix+"_combined_bands.txt", 'a') as f_comb:

        while snr_best > snr_thresh:

            #combine first pair of bands
            start=time.time()
            looptest_comb, looptest_ndet=combine_band_right(test_comb,test_ndet,0)
            SNR=calc_fisher(settings_dict=gen_snr_settings, bands=looptest_comb, dets=looptest_ndet)

            if np.isnan(SNR[0])==True:
                print("NaN config")
                print(count_NaNs+1.)
                count_NaNs+=1.

            snr_best=np.nan_to_num(SNR[0], nan=0.0)
            cost_best=np.nan_to_num(SNR[1], nan=0.0)
            index_best=0
            count_all+=1.
            

            #test combining the rest of the pairs
            for i in range(len(test_comb)-2):

                looptest_comb,looptest_ndet=combine_band_right(test_comb,test_ndet,i+1)

                SNR=calc_fisher(settings_dict=gen_snr_settings,bands=looptest_comb, dets=looptest_ndet)
                count_all+=1.

                if np.isnan(SNR[0])==True:
                    print("NaN configs")
                    print(count_NaNs+1)
                    count_NaNs+=1.


                if np.nan_to_num(SNR[0], nan=0.0) > snr_best:

                    snr_best=np.nan_to_num(SNR[0], nan=0.0)
                    cost_best=np.nan_to_num(SNR[1], nan=0.0)
                    index_best=i+1
                
            end=time.time()
            print("time for one for loop:")
            print(end-start)

            if snr_best > snr_thresh:
                #if haven't reached the threshold yet
                SNRs.append(snr_best)
                freq.append(test_comb[index_best])
                costs.append(cost_best)

                print("remove cached dictionary")
                os.remove(bolocalc_sensfile_path+prefix+"sens_out.npy")

                #create the "best" band combination for the next while loop
                new_test_comb, new_test_ndet=combine_band_right(test_comb,test_ndet,index_best)
                test_comb=np.copy(new_test_comb)
                test_ndet=np.copy(new_test_ndet)
                j+=1
                print("Bands combined so far:")
                print(j)
                print("combined freq: "+str(freq))
                print("SNRs: "+str(SNRs))

                #save information
                np.savetxt(f_comb,np.column_stack([index_best, snr_best, cost_best]),fmt='%1.1f,%1.10f,%1.10f')
                np.savetxt(output_path+prefix+"_optimized_bands.txt",np.column_stack([test_comb, test_ndet]), fmt='%1.1f,%1.1f,%1.1f')


    print("final answer:")
    optimized_bands=np.copy(test_comb)
    optimized_dets=np.copy(test_ndet)

    print(optimized_bands)
    print(optimized_dets)

    snr_final=calc_fisher(settings_dict=gen_snr_settings,bands=optimized_bands, dets=optimized_dets)

    print("final SNR:")
    print(snr_final[0])
    print("final cost:")
    print(snr_final[1])
    print("total NaNs")
    print(count_NaNs)
    print("total fisher calculations")
    print(count_all)

    return optimized_bands, optimized_dets
