# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 13:14:48 2017
Read in pickled files, zeropad, and plot the FFT. 

pickled fields are: 

"Amylase," "Calcium," "Chloride," "Creatinine," "Glucose," "Phosphorus," "Potassium,"
"Sodium," "Urea," "Uric Acid," "Specific Gravity,"

"air_1," "air_2"
"sample_1," "sample_2"...."sample_7"
"water_1," "water_2"

"timestamp," "spectrometerNumber," "sessionID," "sampleNumber," 
"resolution," "lighSourceNumber," "lavNumber," "duration," "deviceID," "date," "data"


@author: DHendricks
"""
import scipy as sp
import matplotlib.pyplot as plt
import pickle
#import pandas as pd
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan


#%%
#02_Ideal_Backgrounds\Pickled_Background
with open('Cleaned_Pickles/'+'Glucose_D_norm142_Scans.p', 'rb') as fp:
#with open('Data/2017_09_01/Pickled_Data/'+'A000059_082963.p', 'rb') as fp:
    data = pickle.load(fp)

#bb = data['air']
#cc = data['ideal']
