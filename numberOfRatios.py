# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:36:14 2017



@author: DHendricks
"""
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
pi,sin,cos,tan = sp.pi,sp.sin,sp.cos,sp.tan


#aa = sp.array([1,2,2,2,3])
#aa = sp.array([1,2,3,4,5])
aa = sp.rand(43)
bb = sp.transpose(aa)
a,b = sp.meshgrid(aa,bb)

c = a*b

#dd =sp.triu(c,1)
ans = c[sp.triu_indices(len(aa),1)]
print(len(ans))

