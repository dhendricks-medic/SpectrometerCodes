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

'x' is optical path difference (in meters)

@author: DHendricks
"""
import scipy as sp
#import matplotlib.pyplot as plt
import pickle
#import pandas as pd
from os import listdir
from os.path import isfile, join, isdir
import spectrometerFunctions as sF
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan
from math import isnan

folderPath = 'Pickles/'
whichLAV = 'A'

target = 'Glucose'
#target = 'Urea'
#target = 'Specific Gravity'
#target = 'Uric Acid'
#target = 'Amylase'
#target = 'Calcium'  
#target = 'Chloride' 
#target = 'Creatinine' ####
#target = 'Phosphorus'
#target = 'Potassium'
#target = 'Sodium'

#
ifSave = True
#ifSave = False
saveName = target+'_'+whichLAV+'_relativeRatioAll'

ratioIndPairs = [
#        [53,68],
#        [95,105],
#        [140,180],
#        [267,272], 
#        [316,330],

#        [40,80],
#        [100,110],
#        [135,200],
#        [260,280], 
#        [303,360],
    [0,450]
             ]
ratioInds = []
for ii in range(len(ratioIndPairs)):
    ratioInds.extend(sp.arange(ratioIndPairs[ii][0],ratioIndPairs[ii][1],2))



def findRatios(yy):
    """ Find ratios between all values passed in through yy"""

    yyPrime = sp.transpose(yy)
    y,yPrime = sp.meshgrid(yy,yyPrime)
    c = y/yPrime
    
    ratios = c[sp.triu_indices(len(yy),1)]
    return(ratios)

allFolders = [f for f in listdir(folderPath) if isdir(join(folderPath, f))]  
print(allFolders)

sampleNames = ['sample_1','sample_2','sample_3']
airNames = ['air_1','air_2']
waterNames = ['water_1','water_2']

idealScanPath = '02_Ideal_Backgrounds/PickledIdeal/'

idealAir,idealBackground = sF.loadIdealFiles(whichLAV)
idealRatio = idealAir/idealBackground # This is what is used to scale the ideal backgroud by air. 

startFolder = {'A':15, 'B':15, 'C':8, 'D':8} ##Note for C, D, only use from folder 8 on (sept 21st) b/c change out pump.## Note for A, B, only use from folder 15 on (oct 2), b/c changed out fiber width.

normalizedScans = []
targetData = []

endishFolder = 10 ## This is just for speed purposes
#for qq in range(startFolder,endishFolder):  ## Step through all folders
for qq in range(startFolder[whichLAV],len(allFolders)):  ## Step through all folders
    singleFolderPath = folderPath+allFolders[qq]+'/'
    allFiles = [f for f in listdir(singleFolderPath) if isfile(join(singleFolderPath, f))]  
    print('folder name is '+singleFolderPath)
    fileNames = []
    for jj in range(len(allFiles)):
        if allFiles[jj].endswith('.p') and allFiles[jj][0] == whichLAV:
            fileNames.append(allFiles[jj])
    numFiles = len(fileNames)

    print('numFiles is '+str(numFiles))
    useFile = []
#    for jj in [2,3]:   #####################################################
    for jj in range(len(fileNames)): ## Step through each file in folder 
#        useFileSingle = True  # This flag determines if I will use the sample.

        print(fileNames[jj])
        with open(singleFolderPath+fileNames[jj], 'rb') as fp:
            data = pickle.load(fp)
            data['x'] = data['x']*100 # Convert to cm. 

        #%% Have the data,now do something to it.
        sampleFFTs = []
        airFFTs = []
        includeScan = True  # If air and water scans checkout, then change to true, and add the data 
        for ii in range(len(sampleNames)):  # Check if scan is liquid
            try:
                fftX,fftY = sF.take_fft(data['x'],data[sampleNames[ii]],pad_x_times = 2)
                fftYInterp = sF.reSample(fftX,fftY)  # Interpolate standard X values
                scanType = sF.isScanWaterOrAir(fftYInterp,idealAir,idealBackground)
                sampleFFTs.append(fftYInterp)
                
                if scanType != 'liquid':
                    includeScan = False
            except:
                includeScan = False
                pass


        if isnan(data[target]) == True:
            includeScan = False

        if includeScan == True:
            ratios = findRatios(sp.mean(sp.asarray(sampleFFTs),axis = 0)[ratioInds])
#            normalizedScans.append(sp.mean(sp.asarray(sampleFFTs),axis = 0) / backgroundPrime)  ## For now, dont scale background.
            normalizedScans.append(ratios)  ## For now, dont scale background.
            targetData.append(data[target])

#%%
if ifSave == True:
    with open('MLModels/'+saveName+'_Scans.p', 'wb') as outfile:
        pickle.dump(sp.asarray(normalizedScans),outfile)
    
    with open('MLModels/'+saveName+'_Vals.p', 'wb') as outfile:
        pickle.dump(sp.asarray(targetData),outfile)
