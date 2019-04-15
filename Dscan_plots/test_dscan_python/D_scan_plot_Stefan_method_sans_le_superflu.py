# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:00:51 2018

@author: ouille
"""

#Librairies used
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy as sc

#Load data
#indicate the file path :
filepath = r'Z:\Laser\D-Scan\HCF\20171110'
#what's the name of the file?
filename= '1200mbar_3.8mJ'
#NOTE:the figure will be saved in the same folder


f_spectrum  =  filepath + '\\' + filename + '_spectrum.csv'
f_pulse     =  filepath +  '\\' + filename + '_retrieved_pulse.csv'
f_retrieved = filepath +  '\\' + filename + '_retrieved_dscan.csv'
f_measured  = filepath +  '\\' + filename + '_measured_dscan.csv'


#Spectral intensity and phase
raw_data = np.genfromtxt(f_spectrum, delimiter=',', dtype=float)
spectrum_wl= np.round(raw_data[:,0],4)   #list of wavelengths
spectrum_I = np.round(raw_data[:,1],4)  #spectral intensity
spectrum_phase = np.round( raw_data[:,2] ,4) #spectral phase


#Temporal intensity and phase
raw_data = np.genfromtxt(f_pulse, delimiter=',', dtype=float)  
pulse_T = raw_data[:,0] #list of times
pulse_Env = raw_data[:,1]  #Electric field envelope
pulse_Car = raw_data[:,2]   #Electric field (carrier)



######## ADD TOD

plt.figure()

####ftAxis function
def ftAxis(nPoints, nuMax):
    deltaNu = 2*nuMax/nPoints
    deltaT = 1/(2*nuMax)
    nu = np.arange(-nuMax,nuMax-(2*nuMax/nPoints),2*nuMax/nPoints)
    t = np.arange(-nPoints/2*deltaT,(nPoints/2-1)*deltaT,deltaT)
    return [nu, t]

###
nPoints = 2**12
[nu,t] = ftAxis(nPoints,5)


c0=299.792458 #nm/fs
##########
lowwl=450
highwl=1050
nu2 =  c0/lowwl 
nu1 =  c0/highwl 



####### interpolate
end= np.size(nu)-1
spectrum_I_nu=[]
spectrum_I_nu = np.interp(nu[nPoints/2 + 1 : end], c0/spectrum_wl[::-1] , spectrum_I[::-1]*(spectrum_wl[::-1]**2))
spectrum_phase_nu = np.interp( nu[nPoints/2 +1:end], c0/spectrum_wl[::-1], spectrum_phase[::-1])



##########
lambdac = 780. #nm
nuc = c0/lambdac  # fs-1
GD= 0.
GDD= 0.
TOD= 5.
FOD= 0.
correctionphase = GD*2*np.pi*(nu[nPoints/2 + 1:end]-nuc) + GDD/2*(2*np.pi*(nu[nPoints/2 +1:end]-nuc))**2 + TOD/6*(2*np.pi*(nu[nPoints/2 +1:end]-nuc))**3 + FOD/24*(2*np.pi*(nu[nPoints/2 +1:end]-nuc))**4


################
fig = plt.figure(figsize=(8,8))
gs = GridSpec(1, 1, wspace=0, hspace=0)

s = plt.subplot(gs[0,0])
s.plot(spectrum_wl,spectrum_I, 'navy')
s.set_xlabel('Wavelength (nm)')
s.set_ylabel('Spectral Intensity (arb.u.)')
s.set_xlim([lowwl,highwl])
s.set_ylim([0,1])

ph=s.twinx()
ph=plt.plot(c0/nu[nPoints/2+1:end],spectrum_phase_nu,'purple', label ='initial phase')
ph=plt.plot(c0/nu[nPoints/2 + 1:end],spectrum_phase_nu-correctionphase,'r', label='corrected phase')
ph=plt.plot(c0/nu[nPoints/2 + 1:end],correctionphase,'--', color='orange', label='phase correction applied')
ph=plt.ylim([-4,4])
ph=plt.xlim([lowwl,highwl])
ph = plt.legend(loc='upper right', fontsize='large')


#### put together symmetric spectrum over both negative and positive frequencies
spectrum_I_nu = np.append(np.zeros(nPoints/2),spectrum_I_nu)
spectrum_phase_nu = np.append(np.zeros(nPoints/2),spectrum_phase_nu)
spectrum_phase_nu_corr = spectrum_phase_nu - np.append(np.zeros(nPoints/2),correctionphase)

##def fonction ft :
def ft(data):
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(data)))

spectrum_I_nu = np.nan_to_num(spectrum_I_nu)
spectrum_phase_nu = np.nan_to_num(spectrum_phase_nu)
pulse_t = ft( np.sqrt(spectrum_I_nu)*np.exp((-1j)*(spectrum_phase_nu-np.pi/2)) )
pulse_I = (pulse_t *np.conj(pulse_t)).real
maxI = max(pulse_I)
maxind = 2039
pulse_tau = np.diff(t[abs(np.diff(np.sign(pulse_I - 0.5*(max(pulse_I)))))==2])
print('Measured pulse duration = ' + str(pulse_tau) + ' fs FWHM.')
spectrum_phase_nu_corr = np.nan_to_num(spectrum_phase_nu_corr)
pulse_t_corr = ft( np.sqrt(spectrum_I_nu)*np.exp(-1j*(spectrum_phase_nu_corr-np.pi/2)) )
pulse_I_corr = (pulse_t_corr *np.conj(pulse_t_corr)).real
pulse_tau_corr = np.diff(t[abs(np.diff(np.sign(pulse_I_corr - 0.5*(max(pulse_I_corr)))))==2])
print('Pulse duration with corrected phase = '+str(pulse_tau_corr)+' fs FWHM,')
print('with peak intensity gain of a factor = '+str(max(pulse_I_corr)/max(pulse_I))+'.')


pulse_t_FTL = ft( np.sqrt(spectrum_I_nu)*np.exp(-1j*(-np.pi/2)) )
pulse_I_FTL = (pulse_t_FTL *np.conj(pulse_t_FTL)).real
pulse_tau_FTL = np.diff(t[abs(np.diff(np.sign(pulse_I_FTL - 0.5*(max(pulse_I_FTL)))))==2])
print('FTL duration = '+str(pulse_tau_FTL)+' fs FWHM.')
print('with peak intensity gain of a factor = '+str(max(pulse_I_FTL)/max(pulse_I))+'.')


plt.figure()
plt.plot(t[2:np.size(t)],pulse_I,'navy', label='initial profile')
plt.plot(t[2:np.size(t)],pulse_I_corr,'r', label='corrected phase')
plt.plot(t[2:np.size(t)],pulse_I_FTL, 'pink', label='with a flat phase')
plt.xlim([-30,30])
plt.legend(loc='upper left', fontsize='large')