# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 12:30:50 2018

@author: ouille
"""

#FIGURE 3
#CEP stab + Energy stab + Spectrum stab

import numpy as np
from matplotlib import pyplot as plt

################################################# SPECTRA
filepathSPEC = r'Z:\Laser\Spectra\HCF\2018\05\23'
filenameSPEC ='4.2mJThorlabs-1200mbar0'

#100 mesures, chacune intégrée sur 10 ms donc 10 pulses.
fileSPEC = str(filenameSPEC) + '.dat'
output_dirSPEC =filepathSPEC
filesSPEC = str(filepathSPEC) + "\\" + str(fileSPEC)


dataSPEC = open(filesSPEC).read().split('\n') # probe timing (ns)
wavelength = dataSPEC[8].split('\t')
wavelength = wavelength[0:np.size(wavelength)-1]
correction_factor=dataSPEC[7].split('\t')
correction_factor = correction_factor[0:np.size(correction_factor)-1]
for index, item in enumerate(wavelength):
    wavelength[index] = float(item)
    correction_factor[index] = float(correction_factor[index])
wavelength = np.transpose (wavelength)
correction_factor=np.transpose(correction_factor)
meanSPEC = dataSPEC[10].split('\t')
meanSPEC = meanSPEC[0:np.size(meanSPEC)-1]
for index, item in enumerate(meanSPEC):
    meanSPEC[index] = float(item)
meanSPEC = np.transpose (meanSPEC)
line9 = dataSPEC[9].split('\t')
nb_of_meas = float(line9[0])

alldatas = np.loadtxt(filesSPEC, dtype=float, skiprows=11)
background = alldatas[0]

    
                         
                         ###############FIGURE : plot the mean spectrum :
plt.figure(figsize=(10,6))
s=14
spectrumstab = plt.plot(wavelength, (meanSPEC-background)*correction_factor/np.max((meanSPEC-background)*correction_factor), linewidth = 1, color='navy')
spectrumstab = plt.xlim([450,1050])
Spectrumstab = plt.ylabel('intensity (a.u.)', size=s)
Spectrumstab = plt.xlabel('wavelength (nm)', size=s)
spectrumstab = plt.tick_params(labelsize=s)
#plt.savefig(r'Z:\Laser\Spectra\HCF\2019\01\15' + '\\' + str(filenameSPEC) + '.png', dpi=800)