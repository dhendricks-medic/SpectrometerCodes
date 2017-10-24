# -*- coding: utf-8 -*-
"""
Created on Fri Oct  20 13:14:48 2017
Feature Extraction 

@author: DHendricks
"""
import scipy as sp
import pickle
from os import listdir
from os.path import isfile, join, isdir
import spectrometerFunctions as sF
from math import isnan
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan

def findRatios(yy):
    """ Find ratios between all values passed in through yy"""
    yyPrime = sp.transpose(yy)
    y,yPrime = sp.meshgrid(yy,yyPrime)
    c = y/yPrime
    ratios = c[sp.triu_indices(len(yy),1)]
    return(ratios)


def doThing():
    folderPath = 'Pickles/'
    whichLAV = 'A'
    
    targets = ['Glucose', 'Urea', 'Specific Gravity', 'Uric Acid', 'Amylase', 'Calcium',
               'Chloride', 'Creatinine', 'Phosphorus', 'Potassium', 'Sodium' ]
    target = 'Glucose'
    
    modelID = '_norm142'
    
#    numBadScans = {'A':0, 'B':0, 'C':0, 'D':0} ## Find way to make this part of class....

    aa = extractFeatures(folderPath,whichLAV,target,modelID)
    aa.takeFFTAndCheckScans(ifSave= True)

    whichLAV = 'B'
    bb = extractFeatures(folderPath,whichLAV,target,modelID)
    bb.takeFFTAndCheckScans(ifSave= True)

    whichLAV = 'C'
    cc = extractFeatures(folderPath,whichLAV,target,modelID)
    cc.takeFFTAndCheckScans(ifSave= True)

    whichLAV = 'D'
    dd = extractFeatures(folderPath,whichLAV,target,modelID)
    dd.takeFFTAndCheckScans(ifSave= True)


#    aa.takeFFTAndCheckScans(ifSave= False)

class extractFeatures(object):
    """
    Will simplify process
    """
    def __init__(self,folderPath,whichLAV,target,modelID):
        self.saveName = target+'_'+whichLAV+ modelID
        self.folderPath = folderPath
        self.whichLAV = whichLAV
        self.target = target
        
    def takeFFTAndCheckScans(self,ifSave=False):
        """ This function will walk through every folder and take the FFTs of
        each sample and air scan, and check if the scans are useable. The results
        will be saved in a new pickle
        
        Need to add where it checks if file exists in save location, and if so skips it. only want to 
        do new files. 
        
        should keep track of how many scans I skip, which could tell me about if slot width is good or bad
        run through all LAVS, not just one of them. 
        """

        allFolders = [f for f in listdir(self.folderPath) if isdir(join(self.folderPath, f))]  
        
        sampleNames = ['sample_1','sample_2','sample_3']
        airNames = ['air_1','air_2']
        waterNames = ['water_1','water_2']
        
#        idealScanPath = '02_Ideal_Backgrounds/PickledIdeal/'
        
        self.idealAir,self.idealBackground = sF.loadIdealFiles(self.whichLAV)
        self.idealRatio = self.idealAir/self.idealBackground # This is what is used to scale the ideal backgroud by air. 

        startFolder = {'A':15, 'B':15, 'C':8, 'D':8} ##Note for C, D, only use from folder 8 on (sept 21st) b/c change out pump. For A, B, only use from folder 15 on (oct 2), b/c changed out fiber width.

#        normalizedScans = []
#        targetData = []
        self.badScans = []
        endishFolder = 10 ## This is just for speed purposes
#        #for qq in range(startFolder[self.whichLAV],endishFolder):  ## Step through only a few folders for debugging
        for qq in range(startFolder[self.whichLAV],len(allFolders)):  ## Step through all folders
#        for qq in [-1]: 
            singleFolderPath = self.folderPath+allFolders[qq]+'/'
            allFiles = [f for f in listdir(singleFolderPath) if isfile(join(singleFolderPath, f))]  
            print('folder name is '+singleFolderPath)
            fileNames = []
            for jj in range(len(allFiles)):
                if allFiles[jj].endswith('.p') and allFiles[jj][0] == self.whichLAV:
                    fileNames.append(allFiles[jj])
            numFiles = len(fileNames)
            print('numFiles is '+str(numFiles))

#            for jj in [27,28]:   #####################################################
            for jj in range(len(fileNames)): ## Step through each file in folder 

                print(fileNames[jj])
                with open(singleFolderPath+fileNames[jj], 'rb') as fp:
                    data = pickle.load(fp)
                    data['x'] = data['x']*100 # Convert to cm. 

                includeScan = True  # If air and water scans checkout, then change to true, and add the data 
                sampleFFTs = []
                for ii in range(len(sampleNames)):  # Check if scan is liquid
                    try:
                        fftX,fftY = sF.take_fft(data['x'],data[sampleNames[ii]],pad_x_times = 2)
                        fftYInterp = sF.reSample(fftX,fftY)  # Interpolate standard X values
                        scanType = sF.isScanWaterOrAir(fftYInterp,self.idealAir,self.idealBackground)
                        sampleFFTs.append(fftYInterp)
                        if scanType != 'liquid':
                            includeScan = False
                            self.badScans.append(fileNames[jj])
                    except:
                        includeScan = False
                        self.badScans.append(fileNames[jj])
                        pass

                if includeScan == False: ## If I know sample is bad, then don't waste time on checking air. I already can't use it. 
                    pass
                else:
                    airFFTs= []
                    for ii in range(len(airNames)):  # Check if air is air
                        try:
                            fftX,fftY = sF.take_fft(data['x'],data[airNames[ii]],pad_x_times = 2)
                            fftYInterp = sF.reSample(fftX,fftY)  # Interpolate standard X values
                            scanType = sF.isScanWaterOrAir(fftYInterp,self.idealAir,self.idealBackground)
                            airFFTs.append(fftYInterp)
                            if scanType != 'air':
                                self.badScans.append(fileNames[jj])
                                includeScan = False
                        except:
                            self.badScans.append(fileNames[jj])
                            includeScan = False
                            pass

                if includeScan == True: ################## FIX THIS HERE. WHAT DO I SAVE? NEW DICTIONARY. Crate new Dictionary with all files. 
                    data['isDataGood'] = True
                    for ii in range(len(sampleNames)):  # 
                        data[sampleNames[ii]+'_FFT'] = sampleFFTs[ii]
                    for ii in range(len(airNames)):  # 
                        data[airNames[ii]+'_FFT'] = sampleFFTs[ii]
                else:
                    data['isDataGood'] = False

#                if ifSave == True: # only save if trigger is set to save.
#                    if includeScan == True: # Only save the good scans....
#                        with open('Cleaned_Pickles/'+fileNames[jj][:-2]+'_FFT'+'.p', 'wb') as outfile:
#                            pickle.dump(data,outfile)
                if ifSave == True: # only save if trigger is set to save.
                    with open('All_Cleaned_Pickles/'+fileNames[jj][:-2]+'_FFT'+'.p', 'wb') as outfile:
                        pickle.dump(data,outfile)

        print('Number of bad scans are: '+str(len(self.badScans)))

if __name__ == '__main__':
    doThing()


""" This is code I want to keep for future. 

                #### I do need to save the FFTd values so I can save them as new pickle.
##                pickle.dump(sp.asarray(normalizedScans),outfile)
                pickle.dump(data,outfile)

##            with open('MLModels/'+saveName+'_Vals.p', 'wb') as outfile:
###                pickle.dump(sp.asarray(targetData),outfile)
##                pickle.dump(data,outfile)




                    #                    # First stab is to do a frequency by frequency scaling. 
#                    airPrime = sp.mean(sp.asarray(airFFTs),axis = 0)
#                    backgroundPrime = airPrime/idealRatio
#                    backgroundDoublePrime = idealBackground/idealBackground[142]*(sp.mean(sp.asarray(sampleFFTs),axis = 0)[142])
#        #            amplitudeAir = max(sp.mean(sp.asarray(airFFTs),axis = 0))
#        #            normalizedScans.append(sp.mean(sp.asarray(sampleFFTs),axis = 0) / idealBackground)  ## For now, dont scale background.
#        
#        #            normalizedScans.append(sp.mean(sp.asarray(sampleFFTs),axis = 0) / backgroundPrime)  ## For now, dont scale background.
#                    normalizedScans.append(sp.mean(sp.asarray(sampleFFTs),axis = 0) / backgroundDoublePrime)  ## Scale idealbackground by index 142, and see what happens. 
#                    targetData.append(data[target])


#    #%%

#                if isnan(data[target]) == True:
#                    includeScan = False


"""