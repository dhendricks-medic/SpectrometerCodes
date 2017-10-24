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
import pickle
#import matplotlib.pyplot as plt
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan


def take_fft(x,y,pad_x_times=[]):
    """
    Python has FFT functions, but I want things in certain formats, so this function
    was written to cater to my needs.
    
    Parameters
    ----------
    
    :x: 1-d array or list, the x data from the raw interferrogram, in cm
    :y: 1-d array or list, the y data from the raw interferrogram
    :pad_x_times: int, the number of times to zero pad the interferrogram. The
        default is nothing, ie, pad_x_times = []
        
    Returns
    -------
    :fftX: numpy array, The wave numbers associated with the FFT'd data.
    :fftY: numpy array, The Y coeeficients associated with the FFT'd data.
    
    """
    len_array = len(x)

    if int(sp.log2(len_array)) == sp.log2(len_array):  # check if already is power of 2. If not power of 2, then num Samples needs to be next highest power of 2
        pow2 = int(sp.log2(len_array))
        numSamples = 2**pow2
        y_pad = y*1
    else:
        pow2 = int(sp.log2(len_array))+1  # Find next largest power of 2, so I can truncate data at a power of 2
        numSamples = 2**pow2
        len_to_pad = numSamples-len_array
#            self.y_pad = sp.pad(self.y,int(len_to_pad/2),mode = 'edge')### Currently non zero edges. Look at adding zeros, or at tapering off. 
        y_pad = sp.pad(y,int(len_to_pad/2),mode = 'linear_ramp')### Currently non zero edges. Look at adding zeros, or at tapering off. 
    
    if pad_x_times != []:
        len_to_pad = int(len(y_pad)/2*pad_x_times)
#            self.y_pad = sp.pad(self.y,int(len_to_pad/2),mode = 'edge')### Currently non zero edges. Look at adding zeros, or at tapering off. 
        y_pad = sp.pad(y,int(len_to_pad/2),mode = 'linear_ramp')### Currently non zero edges. Look at adding zeros, or at tapering off. 
        numSamples = len(y_pad)  

    wave_num = 1/(sp.mean(sp.diff(x))) # sample resolution

    fft_raw = sp.fft(y_pad)
    fft_y_full = 2.0/numSamples*abs(fft_raw[:int(numSamples/2)])

    fft_wavenum = sp.linspace(1e-12,1.0/(2.0*(1/wave_num)),numSamples/2)#/1e6

    x1 = 3500  # Must be << min('standardXValues') to avoid interp problems.
    x1_ind = sp.argmin(sp.absolute(fft_wavenum-x1))
    x2 = 9900 # Must be >> max('standardXValues') to avoid interp problems.
    x2_ind = sp.argmin(sp.absolute(fft_wavenum-x2))

    fftX = fft_wavenum[x1_ind:x2_ind]
    fftY = fft_y_full[x1_ind:x2_ind]

    return(fftX,fftY)


def isScanWaterOrAir(fftY, idealAir, idealBackground, waterCorrThresh = 0.02, airCorrThresh = 0.02):
    """
    Scale the fft to the range -1 => 1, then correlation with water signal, and with air signal
    and see which is higher.
    
    Will need to convolve with perfect air and water signals from the particular LAV
    in question, because of minor setup differences. 
    
    Parameters
    ----------
    :waterCorrThresh: float, 0-1. How close must the cross correlation between 
        the sample scan and the ideal water scan be to the auto-correlation of the 
        ideal water scan, in order for me to believe a good sample was taken. 
        default is 2%. This is arbitrary at the moment, and needs to be. A value
        0f 0.02 means the difference between the two correlations is < +- 2%.

    :AirCorrThresh: float, 0-1. How close must the cross correlation between 
        the sample scan and the ideal air scan be to the auto-correlation of the 
        ideal air scan, in order for me to believe a good air scan was taken. 
        default is 2%. This is arbitrary at the moment, and needs to be figured.
        A value 0f 0.02 means the difference between the two correlations is < +- 2%.
    Returns
    -------
    :fftType: str, tells whether we beleive the fft is a good scan of liquitd,
        air, or if it was a bad scan. The returns are "liquid", "air', or "bad". 
    """
    yMax = max(fftY)
    scaledFFTY = fftY/yMax*2-1 # set to range -1 to 1.
    
    yMaxIdealAir = max(idealAir)
    scaledIdealAir = idealAir/yMaxIdealAir*2-1
    
    yMaxIdealBackground = max(idealBackground)
    scaledIdealBackground = idealBackground/yMaxIdealBackground*2-1

    testWater = sp.correlate(scaledFFTY,scaledIdealBackground,mode = 'valid')
    testAir = sp.correlate(scaledFFTY,scaledIdealAir,mode = 'valid')
    
    perfectWater = sp.correlate(scaledIdealBackground,scaledIdealBackground,mode = 'valid')
    waterCorr = perfectWater*(waterCorrThresh)

    perfectAir = sp.correlate(scaledIdealAir,scaledIdealAir,mode = 'valid')
    airCorr = perfectAir*(airCorrThresh)
    

    if sp.absolute(testWater - perfectWater) < waterCorr:
        fftType = 'liquid'
    elif sp.absolute(testAir - perfectAir) < airCorr:
        fftType = 'air'
    else:
        fftType = 'bad'
    return(fftType)


def loadIdealFiles(whichLAV):
    """
    Given which LAV is being used, load the "ideal air" and "ideal background"
    to use for normalizing scans. 
    """
    initPath = '02_Ideal_Backgrounds/PickledIdeal/'
    if whichLAV == 'D':
        idealPath = initPath+'1503018_k116130474_1.5_D.p'
    elif whichLAV == 'C':
        idealPath = initPath+'1611034_k1161330442_1_C.p'
    elif whichLAV == 'B':
        idealPath = initPath+'1612006_k116130565_1.5_B.p'
    elif whichLAV == 'A':
        idealPath = initPath+'1611038_k116130587_1_A.p'
    else:
        raise ValueError('Unknown LAV sent to "loadIdealFiles"')
        
    with open(idealPath, 'rb') as fp:
        data = pickle.load(fp)
    idealAir = data['air']
    idealBackground = data['ideal']
    return(idealAir,idealBackground)

def compareScans(scan1,scan2):
    """
    scan1 and scan2 must be arrays, containing the y axis data from each scan
    """
    diff = scan1-scan2
    mean_diff = sp.mean(diff)
    print('mean_diff = '+str(mean_diff))

def reSample(fftX,fftY):
    standardXValues = sp.linspace(3600,9500,500) ## This is what all files will be resampled to. 

    from scipy.interpolate import splrep, splev
    interpFunction = splrep(fftX,fftY)
    newY = splev(standardXValues, interpFunction)
    return(newY)

