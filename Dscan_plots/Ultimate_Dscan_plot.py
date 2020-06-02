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
outdir = r'C:\Users\ouille\Desktop\Generated_Files'
P = '1200' #mbar
E = '4' #mJ


#Load files
filepath = r'Z:\Laser\D-Scan\HCF\20200303'
filename = '1200mbar-3_4mJThorlabs'



f_spectrum  =  filepath + '\\' + filename + '_spectrum.csv'
f_pulse     =  filepath +  '\\' + filename + '_retrieved_pulse.csv'
f_retrieved = filepath +  '\\' + filename + '_retrieved_dscan.csv'
f_measured  = filepath +  '\\' + filename + '_measured_dscan.csv'

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

#load retrieved trace
data=  np.genfromtxt(f_retrieved, dtype=float)
data=data[2:np.size(data)]
fheader = np.genfromtxt(f_retrieved, dtype=None)
nWL = int(fheader[0][4:7])
nI= int(fheader[1][3:6])
ret_SHG_wl = data[0:nWL]
ret_insertions = data[nWL:(nWL+nI)]
ret_dscan = np.reshape(data[nWL+nI:np.size(data)],(nI, nWL))

#Load measured trace
data = np.genfromtxt(f_measured, dtype=float)
data=data[2:np.size(data)]
fheader = np.genfromtxt(f_retrieved, dtype=None)
nWL = int(fheader[0][4:7])
nI= int(fheader[1][3:6])
meas_SHG_wl = data[0:nWL]
meas_insertions = data[nWL:(nWL+nI)]
meas_dscan = np.reshape(data[nWL+nI:np.size(data)],(nI, nWL))






#General plot parameters
fig = plt.figure(figsize=(9,8))
gs = GridSpec(2, 2, wspace=0.4, hspace=0.4)

#SUBPLOT : Spectral intensity and phase
spectral = plt.subplot(gs[1, 1])
spectral.plot(spectrum_wl, spectrum_I, color='navy')
spectral.set_ylim([0,1])
spectral.set_xlabel('Wavelength (nm)')
spectral.set_ylabel('Spectral Intensity (a.u.)', color='navy')
spectral.tick_params(axis='y', colors='navy')
spectral1 = spectral.twinx()
spectral1.plot(spectrum_wl, spectrum_phase, 'r--')
spectral1.set_ylim([-4,4])
spectral1.tick_params(colors='r')
spectral1.set_ylabel('Spectral phase (rad)', color='r')
spectral1.set_xlim([450,1050])
                     
#SUBPLOT : Temporal intensity and phase
temp = plt.subplot(gs[1,0])
temp = plt.plot(pulse_T, pulse_Car**2, color= 'navy')
temp = plt.xlabel('Time (fs)')
temp=plt.xlim([-40,40])
plt.ylim([0,1.05])
temp = plt.title('Retrieved temporal shape', weight='bold')
temp = plt.ylabel('Temporal Intensity (a.u.)')
     
#SUBPLOT : measured dscan trace                
measured = plt.subplot(gs[0,0])
measured = plt.pcolor(meas_SHG_wl, meas_insertions, meas_dscan/np.max(meas_dscan), cmap = 'jet')
measured = plt.xlabel('Wavelength (nm)');
measured = plt.ylabel('Insertion (mm)');
measured = plt.title('Measured', weight='bold');

                   
#SUBPLOT : retireved dscan trace                    
retrieved = plt.subplot(gs[0,1])
retrieved = plt.pcolor(ret_SHG_wl, ret_insertions, ret_dscan/np.max(ret_dscan), cmap='jet')
retrieved = plt.xlabel('Wavelength (nm)');
retrieved = plt.ylabel('Insertion (mm)');
retrieved = plt.title('Retrieved', weight='bold');
              
                     
                     
           
                     
                     
#FWHM of the temporal intensity profile
#Function to obtain the FWHM :
def FWHM (x, y):
    ymax = np.max(y)
    print('\n', ymax*0.5 , 'max*0.5')
    for i, value in enumerate (y):
        if value == ymax:
            index_max = i
         #   print(index_max)
    for i, value in enumerate ( y[0:index_max] ):
        if value >= 0.5*y[index_max]:
            low = i
            break
    print(y[low] , 'y[low]')
    for i, value in enumerate (y):
        if i > index_max+1:
            if value <= 0.5*y[index_max]:
                high = i
                break    
    print(y[high] , 'y[high]', '\n')
    return (x[high] - x[low])
#interpolation for more precision 
X_new = np.arange (-40,40, 0.0001)
pulse_Car = np.nan_to_num(pulse_Car)
pulse_T = np.nan_to_num(pulse_T)
f = interpolate.interp1d(pulse_T[1:np.size(pulse_T)-2], pulse_Car[1:np.size(pulse_T)-2]**2)
T = FWHM(X_new, f(X_new))
print(str(T) + ' fs FWHM')                
fig.suptitle( str(round(float(T),1)) + 'fs FWHM - P=' + str(P) + 'mbar - E=' + str(E)+ 'mJ  --  date = ' + str(filepath[len(filepath)-8:len(filepath)]) )                     
  

                   
#☺save figure
plt.savefig(str(outdir) + '\\' +str(filepath[len(filepath)-8:len(filepath)]) + '__' +  str(filename) + '.png', dpi=200)