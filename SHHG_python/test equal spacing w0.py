# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 13:50:10 2019

@author: ouille
"""

from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

l = 790*1e-3 #µm
E0 = 1.2398/l
harmonics = (E0)
lambd = l
for i, k in enumerate (np.arange(2,30)):
    harmonics= np.append(harmonics, E0*k)
    lambd = np.append(lambd, l/k)
B = - np.arcsin( lambd/600000 - np.sin(85.3*np.pi/180) )

p = np.arange(1,30)
plt.figure()
plt.plot(p,B*180/np.pi, 'o')
#plt.ylim(10,32)
#plt.xlim(0.47699843206843770, 0.47699860321964116)
plt.ylabel('angle B = arcsin ($\lambda$*$\sigma_0$ - sin A) ')
plt.xlabel('harmonic number')


pos = 469/np.tan(B) #mm
plt.figure()
plt.plot(p,pos , 'o')        
#plt.ylim(10,32)  
#plt.xlim(907.5044, 907.5048)      
plt.ylabel('postion x = 469 mm / tan B')
plt.xlabel('harmonic number')


#file = r'C:\Users\ouille\Desktop\SHHG analysis\HHG calibration\\1402-24fs-refspec.txt'
#test = np.loadtxt (file)
#plt.figure()
#plt.plot(test)
#plt.title ('Ref for stefan calibration')
