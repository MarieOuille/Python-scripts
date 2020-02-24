# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:33:54 2020

@author: ouille

Fringeezz very basic plot
"""



"""
 STRUCTURE OF THE HDF FILES :    
 The hdf5 file is divided into 3 parts : CEP values ; time stamps ; configuration
 Each carrier-enveloppe phase value is recorded as an integer : you must divide by 10 000 to get it in radians
 phase values are recorded in an array
 Each line of this array corresponds to a 20second-long measurement (parameter 'phase points' defined through the Fringeezz software : if phase points = 20 000, each line corresponds to 20sec)
 timestamps : list of times (sec) when each measurement line ended 
 printing timestamps[1] - timestamps[0] should therefore return 20 as there are 20 seconds between two lines of measurements

"""





# libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import welch # PSD using Welch method



## FIRST : EXTRACT DATA 
# option 1 = ask for file location and filename
#print ('Where is the file located?')
#filepath = input('The path of the file is :')  # e.g. : C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\Fringeezz
#print ('What is the name of the hdf5 file?')
#filename = input('The name of the file is :')   # e.g. :  'phase_log_20170512_110754' 

#option 2 : enter it manually here :
filepath = r'Z:\Laser\CEP\20200218'
filename = 'phase_log_20200218_164551'




file= str(filename) + '.hdf5' 

#output directory = new directory created where "fringeezz.py" is located
path = os.getcwd()
if not os.path.exists('analyzed_data'):
    os.mkdir('analyzed_data')
output_dir=path+'\\analyzed_data'

def fringeezz_load_h5_data_log(file, filepath, print_structure=True):
    """ load data log file in HDF5 format from the fringeezz from default directory.
    input: file name, filepath if not default, print_structure=True if you want a print of the hdf5 file structure, else = False
    return: fringeezz configuration as string, timestamps array, phase matrix (rad)
    """
    with h5py.File(os.path.join(filepath,file)) as h5f:
        if print_structure:
            for name in h5f:
                print(h5f[name])
        config = np.array(h5f['configuration'], dtype = '|S1545')
        if 'phases' in h5f:
            phases = np.array(h5f['phases'], dtype = int)
        else:
            phases = float('nan')
        timestamps = np.array(h5f['timestamps'], dtype = int).ravel()
    return config[0].decode('UTF-8'), timestamps, phases
config, timestamps, phases = fringeezz_load_h5_data_log(file,filepath) #calls the function defined above to extract data from your hdf5 file 
print(config) #To get info on the configuration of the Fringeezz at the time of the measurement :

# so now, let's deal with these raw data and turn it into nice basic lists :
    
phases = (phases.ravel())/4096*2*np.pi
time = np.arange (0, np.size(phases)/1000, 0.001)  #divided by 1000 (number of points per second since we have 1kHz rep rate, and with a 0.001 sec interval)

plt.figure(figsize=(10,4))
plt.plot(time/60, phases,  '.', markersize = 1, color = 'navy')
plt.xlabel('time(min)')
plt.ylabel ('phase(rad)')
#plt.xlim(0,np.max(time/60))
plt.xlim(0,0.8)
plt.ylim(np.min(phases),np.max(phases))
plt.title(str(filename))