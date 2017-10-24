# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 11:22:17 2017

load the pickled data and try a model. Dont expect this to do well. Need to 
scale by air scan. 

@author: DHendricks
"""
import scipy as sp
import matplotlib.pyplot as plt
import pickle
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan

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

LAV = 'A'

#saveName = target + '_'+LAV+'_airNorm'
#saveName = target + '_'+LAV+'_relativeRatio0'
saveName = target + '_'+LAV+'_norm142'


with open('MLModels/'+saveName+'_Scans.p', 'rb') as fp:
    data = pickle.load(fp)

with open('MLModels/'+saveName+'_Vals.p', 'rb') as fp:
    target = pickle.load(fp)

print('total number of samples is '+str(len(data)))
#%% Learn

X_train, X_test, y_train, y_test = train_test_split(data, target, test_size=0.33, random_state=42)
regr = LinearRegression()

scores = cross_val_score(regr,data,target,cv = 3)
otherPredict = cross_val_predict(regr,data,target,cv=3)

print(scores)
regr.fit(X_train,y_train)
y_prediction = regr.predict(X_test)

coeefs = regr.coef_
r2 = sp.round_(r2_score(y_test,y_prediction),3)
print('r2 score is ' + str(r2))

maxY = max(y_test)
minY = min(y_test)


#%% PLOT

plt.figure(1,facecolor = 'w',figsize = [10,8])
plt.clf()
plt.subplot(211)
#plt.plot(y_prediction,y_test,'k.',ms=8)
plt.plot(otherPredict,target,'k.',ms=8)


plt.plot([minY,maxY],[minY,maxY],'k--',lw=1)
plt.text(minY,(maxY-minY)*.95+minY,'R Squared value is '+str(r2),fontsize = 15)
plt.grid()
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.xlabel('Predicted Values',fontsize = 20)
plt.ylabel('AU480 Values',fontsize = 20)
plt.title(saveName,fontsize = 25)
#plt.savefig('MLModels/figs/'+saveName+'.png')
#

#plt.figure(2,facecolor = 'w',figsize = [10,8])
#plt.clf()
plt.subplot(212)
plt.plot(coeefs,'b',lw=2)
plt.grid()
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

#plt.ylim(0,2000)
#plt.xlim(-1000,2000)


"""
multipleRuns = False
if multipleRuns== True:
    all_coeefs = []
    rscore = []
    for ii in range(10):
        X_train, X_test, y_train, y_test = train_test_split(data, target, test_size=0.2,)
        
        regr = LinearRegression()
        regr.fit(X_train,y_train)
        y_prediction = regr.predict(X_test)
        
        coeefs = regr.coef_
        all_coeefs.append(coeefs)    
        r2 = sp.round_(r2_score(y_test,y_prediction),3)
        rscore.append(r2)
    
    ave_coeefs = sp.mean(sp.asarray(all_coeefs),axis = 0)
    print('ave_rscore is '+str(sp.mean(rscore)))


"""
