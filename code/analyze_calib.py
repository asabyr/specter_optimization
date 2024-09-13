import numpy as np
import sys

class AnalyzeCalib:

    def __init__(self, fprefix, run=1, nfiles=1, noise_file=""):

        self.fprefix=fprefix
        self.run=run
        self.nfiles=nfiles
        self.noise_file=noise_file
        self.clean_bias()
        
    def clean_bias(self):
        #clean NaNs
        if self.nfiles==1:
            bias_file = np.loadtxt("../calib_files/"+self.fprefix+"Bs_run"+str(self.run)+".txt")

        else:
            for i in range(self.nfiles):
                onefile = np.loadtxt("../calib_files/"+self.fprefix+"Bs_run"+str(self.run)+"_"+str(i)+".txt")

                if i==0:
                    bias_file=onefile
                else:
                    bias_file=np.vstack((bias_file,onefile))

        delta_I=np.loadtxt("../calib_files/"+self.fprefix+"_calib_combinations_run"+str(self.run)+".txt")
        mask=np.invert(np.isnan(bias_file[:,0]))

        print("Fraction of non-NaNs") #fraction of values that are not NaN
        print(np.sum(~np.isnan(bias_file[:,0]))/len(bias_file[:,0]))

        self.delta_I=delta_I[mask]
        self.bias=bias_file[mask,:]
        

    def maximum_bias(self):

        ind=np.where(np.abs(self.bias[:,0])==np.amax(np.abs(self.bias[:,0])))[0]
        print("maximum bias [\%,\sigma]:")
        print(self.bias[ind,0])
        print(self.bias[ind,1])
        
        if len(self.noise_file)>0:
            freq, noise=np.loadtxt("../files/"+self.noise_file, unpack=True)
            print(r"$\Delta I$ for maximum bias [fraction of noise]: ")
            print(self.delta_I[ind]/noise)
        

    def filter_bias(self, bias_range, return_all=False):

        ind=np.where((np.round(self.bias[:,0],1)>bias_range[0]) & (np.round(self.bias[:,0],1)<=bias_range[1]))

        return ind

