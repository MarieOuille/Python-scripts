# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:00:51 2018

@author: ouille
"""

#Librairies used
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import interpolate

# info on where and how to save the image
#outdir = r'C:\Users\ouille\Desktop\codes_python_divers'
outdir = r'Z:\Laser\Compte-rendus\chirped_mirrors_PC1332_test'


#Load files
filepath = r'Z:\Laser\D-Scan\HCF\20191219'
filepath2 =r'Z:\Laser\D-Scan\HCF\20191217'
filepath3 = r'Z:\Laser\D-Scan\HCF\20191217'
filename = '1300mbar-4_35mJ-bis'    #PC70 + PC1332
filename2 = 'PC1332_1300mbar_4_45mJ_3_18fs'   #PC1332
#filename3 =  'PC70_chirped_mirrors_1200mbar_4-45mJ'    #PC70
f_spectrum  =  filepath + '\\' + filename + '_spectrum.csv'
f_pulse     =  filepath +  '\\' + filename + '_retrieved_pulse.csv'
f_spectrum2  =  filepath2 + '\\' + filename2 + '_spectrum.csv'
f_pulse2     =  filepath2 +  '\\' + filename2 + '_retrieved_pulse.csv'
#f_spectrum3  =  filepath3 + '\\' + filename3 + '_spectrum.csv'
#f_pulse3     =  filepath3 +  '\\' + filename3 + '_retrieved_pulse.csv'

#load spectral intensity and phase
raw_data = np.genfromtxt(f_spectrum, delimiter=',', dtype=float)
spectrum_wl=raw_data[:,0]
spectrum_I = raw_data[:,1]
spectrum_phase = raw_data[:,2]

#load temporal intensity and phase
raw_data = np.genfromtxt(f_pulse, delimiter=',', dtype=float)
pulse_T = raw_data[:,0]
pulse_Env = raw_data[:,3]
pulse_Car = raw_data[:,4] 

#load spectral intensity and phase
raw_data2 = np.genfromtxt(f_spectrum2, delimiter=',', dtype=float)
spectrum_wl2=raw_data2[:,0]
spectrum_I2 = raw_data2[:,1]
spectrum_phase2 = raw_data2[:,2]

#load temporal intensity and phase
raw_data2 = np.genfromtxt(f_pulse2, delimiter=',', dtype=float)
pulse_T2 = raw_data2[:,0]
pulse_Env2 = raw_data2[:,3]
pulse_Car2 = raw_data2[:,4] 

##load spectral intensity and phase
#raw_data3 = np.genfromtxt(f_spectrum3, delimiter=',', dtype=float)
#spectrum_wl3=raw_data3[:,0]
#spectrum_I3 = raw_data3[:,1]
#spectrum_phase3 = raw_data3[:,2]
#
##load temporal intensity and phase
#raw_data3 = np.genfromtxt(f_pulse3, delimiter=',', dtype=float)
#pulse_T3 = raw_data3[:,0]
#pulse_Env3 = raw_data3[:,3]
#pulse_Car3 = raw_data3[:,4] 
#
#


#General plot parameters
fig = plt.figure(figsize=(12,8))
spectral = fig.add_subplot(111)
spectral.plot(spectrum_wl, spectrum_I, color='navy',label = 'PC70 & PC1332')
spectral.plot(spectrum_wl2, spectrum_I2, color='cyan',label = 'PC1332')
#spectral.plot(spectrum_wl3, spectrum_I3, color='black',label = 'PC70')
spectral.set_ylim([0,1.1])
spectral.set_xlabel('Wavelength (nm)')
spectral.set_ylabel('Spectral Intensity (a.u.)')
spectral.tick_params(axis='y', colors='navy')
spectral.set_xlim([450,1050])
spectral.legend()      
fig = plt.savefig(str(outdir) + '\\comp_spec_1300mbar.png', dpi=200)  


fig3 = plt.figure(figsize=(12,8))
spectral1 = fig3.add_subplot(111)
spectral1.plot(spectrum_wl, spectrum_phase, 'navy', label = 'PC70 & PC1332')
spectral1.plot(spectrum_wl2, spectrum_phase2, 'cyan', label = 'PC1332')
#spectral1.plot(spectrum_wl3, spectrum_phase3, 'black', label = 'PC70')
spectral1.set_ylim([-4,4])
spectral1.set_ylabel('Spectral phase (rad)')
spectral1.set_xlabel('Wavelength (nm)')
spectral1.set_xlim([450,1050])
spectral1.legend(loc='upper right')
fig3 = plt.savefig(str(outdir) + '\\comp_spectral_phase_1300mbar.png', dpi=200)


fig3zoom = plt.figure(figsize=(12,8))
spectral1zoom = fig3zoom.add_subplot(111)
spectral1zoom.plot(spectrum_wl, spectrum_phase, 'navy', label = 'PC70 & PC1332')
spectral1zoom.plot(spectrum_wl2, spectrum_phase2, 'cyan', label = 'PC1332')
#spectral1zoom.plot(spectrum_wl3, spectrum_phase3, 'black', label = 'PC70')
spectral1zoom.set_ylim([-8,6])
spectral1zoom.set_ylabel('Spectral phase (rad)')
spectral1zoom.set_xlabel('Wavelength (nm)')
spectral1zoom.set_xlim([450,550])
spectral1zoom.legend(loc='upper right')  
fig3zoom = plt.savefig(str(outdir) + '\\comp_spectral_phase_zoom_1300mbar.png', dpi=200)

                     
#SUBPLOT : Temporal intensity and phase
fig2 = plt.figure(figsize=(12,8))
temp = plt.plot(pulse_T, pulse_Car**2, color= 'navy', label = 'PC70 & PC1332')
plt.plot(pulse_T2-2, pulse_Car2**2, color= 'cyan', label = 'PC1332')
#plt.plot(pulse_T3, pulse_Car3**2, color= 'black', label = 'PC70')
temp = plt.xlabel('Time (fs)')
temp=plt.xlim([-40,40])
temp = plt.ylim([0,1.05])
temp = plt.title('Retrieved temporal shape', weight='bold')
temp = plt.ylabel('Temporal Intensity (a.u.)')
temp = plt.legend() 
fig2 = plt.savefig(str(outdir) + '\\comp_temp_1300mbar.png', dpi=200)
                   
  

