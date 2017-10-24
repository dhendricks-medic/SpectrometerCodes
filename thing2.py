# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:40:02 2017



@author: DHendricks
"""
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan
#
#
##from joblib import Parallel, delayed
#import multiprocessing
#    
## what are your inputs, and what operation do you want to 
## perform on each input. For example...
#
##inputs = range(10) 
##def processInput(i):
##	return i * i
#
#num_cores = multiprocessing.cpu_count()
#
#print(num_cores)    
##results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs
#                   

import multiprocessing as mp

def processInput(i):
        return i * i

if __name__ == '__main__':

    # what are your inputs, and what operation do you want to
    # perform on each input. For example...
    inputs = range(1000)
    #  removing processes argument makes the code run on all available cores
    pool = mp.Pool(processes=4)
    results = pool.map(processInput, inputs)
    print(results)