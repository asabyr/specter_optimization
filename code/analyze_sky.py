import numpy as np

class AnalyzeSkyModel:
    
    def __init__(self, fprefix, run=1, nfiles=1, NaNtozero=False):

        self.fprefix=fprefix
        self.run=run
        self.nfiles=nfiles
        self.NaNtozero=NaNtozero
        self.params_snr,self.snr=self.clean_bias(which='sky_SNR_bias')

    def clean_bias(self, which):
        
        #read either one or many files
        if self.nfiles==1:
            bias_file = np.loadtxt("../grid_files/"+self.fprefix+which+"_run"+str(self.run)+".txt")
        else:
            for i in range(self.nfiles):
                onefile = np.loadtxt("../grid_files/"+self.fprefix+which+"_run"+str(self.run)+"_"+str(i)+".txt", unpack=True)
                
                if i==0:
                    bias_file=onefile.copy()
                else:
                    bias_file=np.append(bias_file, onefile.copy(), axis=None)
        
        
        
        params=np.loadtxt("../grid_files/"+self.fprefix+"_sky_combinations_run"+str(self.run)+".txt")
        
        #count NaNs as zeros or mask out
        if self.NaNtozero==False:
            mask=np.invert(np.isnan(bias_file))
            return params[mask,:],bias_file[mask]
        else:
            clean_bias_file=np.nan_to_num(bias_file, copy=True, nan=0.0)
            return params, clean_bias_file
    
    def sort_results(self):
        indices=np.argsort(-self.snr)
        
        self.params_snr_sort,self.snr_sort=self.params_snr[indices],self.snr[indices]

    def maximum_snr_drop(self):

        ind=np.where(self.snr==np.amin(self.snr))[0]
        print("minimum snr:")
        print(self.snr[ind])
        print(r"parameters for maximum drop: ")
        print(self.params_snr[ind])
        
        return self.snr[ind],self.params_snr[ind]

    def maximum_snr_increase(self):

        ind=np.where(self.snr==np.amax(self.snr))[0]
        print("maximum snr:")
        print(self.snr[ind])
        print(r"parameters for maximum increase: ")
        print(self.params_snr[ind])
    
        return self.snr[ind],self.params_snr[ind]
    