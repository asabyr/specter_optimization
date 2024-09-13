import numpy as np
from itertools import product

def Ndet_multichroic(Nmulti,Ntimes):
    '''
    this function outputs a "multichroic" detector combination

    args:
    Nmulti: detector combinations (np.array)
    Ntimes: how many times each element is repeated (np.array)

    output:
    detector combination (np.array)

    e.g.
    det_comb=np.array([1,2,3])
    times=np.array([2,3,4])
    print(Ndet_multichroic(det_comb, times))
    >>[1. 1. 2. 2. 2. 3. 3. 3. 3.]

    '''
    if Nmulti.ndim > 1:
        for i in range(len(Ntimes)):

            if i==0:
                Ndet=np.multiply(Nmulti[:,i],np.ones((len(Nmulti[:,i]),Ntimes[i])).T).T
            else:
                newNdet=np.multiply(Nmulti[:,i],np.ones((len(Nmulti[:,i]),Ntimes[i])).T).T
                Ndet=np.column_stack((Ndet, newNdet))

    else:
        Ndet=[]
        for i in range(len(Ntimes)):
            Ndet=np.hstack((Ndet, Nmulti[i]*np.ones(Ntimes[i])))
    return Ndet

def make_det_combinations(mins, maxs, ns):
    '''
    this function creates all possible detector combinations
    based on a set of inputs

    args:
    mins: min detector count for each band/multichroic band (np.array)
    maxs: max detector count for each band/multichroic band (np.array)
    ns: # of elements for the corresponding linspace (np.array)

    output:
    detector combinations (np.array)
    
    e.g.
    min_arr=np.array([1,1])
    max_arr=np.array([2,2])
    ns_arr=np.array([2,2])
    print(make_det_combinations(min_arr, max_arr, ns_arr))
    >>[[1. 1.]
    [1. 2.]
    [2. 1.]
    [2. 2.]]
    '''
    det_list=[]
    for i in range(len(mins)):

        det_list.append(np.linspace(mins[i], maxs[i], ns[i]))

    return np.asarray(list(product(*det_list)))
