import numpy as np
import sys

class AnalyzeDets:

    def __init__(self, fprefix, run=1, nfiles=1):

        self.fprefix=fprefix
        self.run=run
        self.nfiles=nfiles
        self.clean_SNR()
        
    def clean_SNR(self):
        #clean NaNs from grid optimization
        if self.nfiles==1:
            SNR_file = np.loadtxt("../grid_files/"+self.fprefix+"SNRs_run"+str(self.run)+".txt")
        else:
            for i in range(self.nfiles):
                onefile = np.loadtxt("../grid_files/"+self.fprefix+"SNRs_run"+str(self.run)+"_"+str(i)+".txt")
                if i==0:
                    SNR_file=onefile
                else:
                    SNR_file=np.vstack((SNR_file,onefile))

        dets=np.loadtxt("../grid_files/"+self.fprefix+"_ndet_combinations_run"+str(self.run)+".txt")
        mask=np.invert(np.isnan(SNR_file[:,0]))

        print("maximum snr:")
        print(np.amax(np.nan_to_num(SNR_file[:,0], copy=True, nan=0, posinf=0)))
        print("Fraction of non-NaNs") #fraction of values that are not NaN
        print(np.sum(~np.isnan(SNR_file[:,0]))/len(SNR_file[:,0]))

        self.dets=dets[mask]
        self.snrs_area=SNR_file[mask,:]
        
        

    def filter_SNR(self, req_snr, return_all=False):

        if np.ndim(req_snr)==0:
            ind=np.where(self.snrs_area[:,0]>req_snr)
        else:
            ind=np.where((self.snrs_area[:,0]>req_snr[0])&(self.snrs_area[:,0]<req_snr[-1]))

        min_cost=np.amin(self.snrs_area[:,1][ind])

        min_cost_ind=np.where(self.snrs_area[:,1][ind]==min_cost)
        min_cost_snr=self.snrs_area[:,0][ind][min_cost_ind]
        min_cost_det=self.dets[ind][min_cost_ind]

        print("min cost & snr & configs:")
        print(min_cost)
        print(min_cost_snr)
        print(min_cost_det)

        if return_all==True:
            return ind, min_cost, min_cost_snr, min_cost_det
        else:
            return ind
        

    def filter_cost(self, req_cost, return_all=False):

        ind=np.where((np.round(self.snrs_area[:,1],3)>req_cost[0]) & (np.round(self.snrs_area[:,1],3)<=req_cost[1]))

        max_snr=np.amax(self.snrs_area[:,0][ind])
        min_snr=np.amin(self.snrs_area[:,0][ind])

        max_snr_ind=np.where(self.snrs_area[:,0][ind]==max_snr)
        min_snr_ind=np.where(self.snrs_area[:,0][ind]==min_snr)

        max_snr_cost=self.snrs_area[:,1][ind][max_snr_ind]
        min_snr_cost=self.snrs_area[:,1][ind][min_snr_ind]

        max_snr_det=self.dets[ind][max_snr_ind]
        min_snr_det=self.dets[ind][min_snr_ind]

        print("max snr & cost & configs:")
        print(max_snr)
        print(max_snr_cost)
        print(max_snr_det)
        

        if return_all==True:
            return ind, [max_snr_det, min_snr_det], [max_snr_ind, min_snr_ind]
        else:
            return ind

