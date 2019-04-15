# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:10:32 2019

@author: ouille
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage


filepath = r'C:\Users\ouille\Desktop\HHG_python'
filename = 'tir237.dat'
outdir = r'C:\Users\ouille\Desktop\HHG_python'
im = np.loadtxt (filepath + '\\' + filename)
plt.figure()
plt.imshow(ndimage.rotate(im[100:650,50:380], 90), cmap='jet')
plt.ylabel('vertical direction')
plt.xlabel('photon energy')


