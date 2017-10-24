# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 15:38:05 2017



@author: DHendricks
"""
import scipy as sp
import matplotlib.pyplot as plt
import pickle
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict
from sklearn.linear_model import LinearRegression 
from sklearn.metrics import r2_score
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan



x = sp.linspace(0,10,100)
data = []
num_samples = 500


from  scipy import random
random.seed(seed = 3)
targets = sp.randn(num_samples)*100
#targets = random.random()*100
#print(targets)


random.seed(seed=None)
for ii in range(num_samples):

#    y.append(sin(x)*(targets[ii]+sp.rand()*10))
    ynew = sin(x)+sp.rand()*10
    ynew[4] += targets[ii]+sp.rand()*10
#    ynew[5] += targets[ii]+sp.rand()*10
#    ynew[6] += targets[ii]+sp.rand()*10

    data.append(ynew) #*(targets[ii]+sp.rand()*10))

data = sp.asarray(data)
X_train, X_test, y_train, y_test = train_test_split(data, targets, test_size=0.33, random_state=42)
regr = LinearRegression()

otherPredict = cross_val_predict(regr,data,targets,cv=4)

#scores = cross_val_score(regr,data,targets,cv = 4)
#print(scores)
regr.fit(X_train,y_train)
y_prediction = regr.predict(X_test)
#
coeefs = regr.coef_
r2 = sp.round_(r2_score(y_test,y_prediction),3)
print('r2 score is ' + str(r2))
#
maxY = max(y_test)
minY = min(y_test)


plt.figure(3,facecolor = 'w')
plt.clf()
plt.subplot(211)
#plt.plot(y_test,y_prediction,'k.')
plt.plot(targets,otherPredict,'k.')
plt.subplot(212)
plt.plot(coeefs)

