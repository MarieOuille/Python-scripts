# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:54:56 2020

@author: ouille

This script generates a .png figure showing the temporal and spectral phases and intensity profiles.
It requires 3 .txt files 
These 3 files are directly generated from the Wizzler software 
(!! It can be tricky : don't write any filename when saving from the Wizzler and point a folder only !!
 Otherwise the software returns an error and crashes)
"""

#%% librairies
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


#%% filepath, to be changed :
filepath = r'C:\Users\ouille\Desktop\converted'
filename = '2020-05-29_14-41-25.51_SRSI_'
outdir = r'C:\Users\ouille\Desktop\converted'   #output directory, to save the figure as png


#%% infos
nums = (filepath + '\\' + filename + 'nums.txt')
data = open(nums, "r")
print(data.read())
data = open(nums, "r")
lines = data.readlines()


#%% figure parameters
fig = plt.figure(figsize=(12,5))
s=12 #labels size for the figures 
fig.text(0.39,0.95,'Wizzler measurement', size=s+2)
fig.text(0.39,0.92, filename[0:22] , size=s-3)
fig.text(0.12,0.83, str(lines[2][0:-1]), size=s+1)
gs = gridspec.GridSpec(1, 2, left=0.05, right=0.95, wspace = 0.3, top = 0.80)


#%% temporal
time = np.loadtxt (filepath + '\\' + filename + 'time.txt', skiprows=1)
t = time[:,0]
t_int = time[:,1]
t_pha = time[:,2]

ax11 = fig.add_subplot(gs[0])
ax11.plot(t, t_int, color='navy')
ax11.set_ylabel('Intensity', color='navy', size=s)
ax11.set_xlabel('Time (fs)', size=s)
ax11.set_xlim(-100,100)
ax11.set_ylim(0,1)
ax11.tick_params(axis='y', color='navy', labelcolor='navy', labelsize=s)
ax11.tick_params(axis='x',labelsize=s)

ax12 = ax11.twinx()
ax12.plot(t, t_pha, color='red')
ax12.set_ylabel('Phase', color='red', size=s)
ax12.set_ylim(-5,5)
ax12.tick_params(axis='y', color='red', labelcolor='red', labelsize=s)



#%%Spectral
spectral = np.loadtxt (filepath + '\\' + filename + 'spectral.txt', skiprows=1)
f = spectral[:,0]
f_int = spectral[:,1]
f_pha = spectral[:,2]
f_fit = spectral[:,3]

#fig2,ax21 = plt.subplots()
ax21 = fig.add_subplot(gs[1])
ax21.plot(3e8/f*1e-3, f_int, color='navy')
ax21.set_ylabel('Intensity', color='navy', size=s)
ax21.set_xlabel('Wavelength (nm)', size=s)
ax21.set_xlim(730,850)
ax21.set_ylim(0,1)
ax21.tick_params(axis='y', color='navy', labelcolor='navy',labelsize=s)
ax21.tick_params(axis='x',labelsize=s)

ax22 = ax21.twinx()
ax22.plot(3e8/f*1e-3, f_pha, color='red',label='phase')
ax22.plot(3e8/f*1e-3, f_fit,  '--', color='pink', label='fit')
ax22.tick_params(axis='y', color='red', labelcolor='red', labelsize=s)
ax22.set_ylabel('Phase', color='red',size=s)
ax22.set_ylim(-3,5)
ax22.legend(fontsize=s)


#%% save figure
fig = plt.savefig(outdir + '\\' + filename + 'figure.png', bboxinches='tight', dpi=300)