# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:02:22 2019

@author: ouille




# step 1 : execute this script
# step 2 : use the output file for calibration with 'calibration_stefan.py'
# step 3 : Use the calibration file to plot calibrated spectra

"""


#%%
#import libraries
import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage



#%% Entrées à changer :

#input = raw data from the MCP camera
#filepath = r'Z:\SHHG\Stefan_and_Marie\2nd-campaign-May2019\20190605\harmonics'
#filepath = r'Z:\SHHG\LILU\20190829\harmonics'
#filename = 'tir310.dat'
#filename = 'tir66.dat'
#
#filepath = r'Z:\SHHG\LILU\20191202\harmonics'
#filename = 'tir171.dat'
#outdir = r'C:\Users\ouille\Desktop\Python_travail_en_cours'

#Après le confinement, avec Jaismeen
#filepath = r'Z:\SHHG\Jaismeen_Marie\20200612\harmonics'
#filename = 'tir460.dat'
#outdir = r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files'

#test pour le scan 20190605, avec chirp
filepath = r'Z:\SHHG\Stefan_and_Marie\2nd-campaign-May2019\20190605\harmonics'
filename = 'tir38.dat'
outdir = r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files'




#%% load and plot :
    
im = np.loadtxt (filepath + '\\' + filename)

#look at the HHG spectrum as recorded
plt.figure()
#plt.imshow(ndimage.rotate(im[20:560,200:380], 90), cmap='jet')
plt.imshow(ndimage.rotate(im, 90), cmap='jet')
plt.ylabel('vertical direction')
plt.xlabel('photon energy')

#integrate vertically :
plt.figure()
im_v = np.zeros(np.size(im[:,0]))
#im_v = np.zeros(560-20)
for i, j in enumerate (im[:,0]):
    im_v[i] = sum (im[i])
#for i, j in enumerate (im[20:560,0]):
#    im_v[i] = sum (im[20+i, 200:380])
plt.plot(im_v)


#save txt
np.savetxt( (outdir+ '\\' + filename[0:-4]+'_VerticallyIntegrated.txt'), im_v) 