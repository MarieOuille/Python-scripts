# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 12:37:38 2018

@author: ouille
"""


#idée = donner des paramètres de dispersion et voir les spectrogrammes


import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import matplotlib.gridspec as gridspec

#useful functions :
def ift (x):
    return (np.fft.fftshift(np.fft.ifft(np.fft.ifftshift( x ))))
def ft(x):
    return (np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x))))
def ftAxis(nPoints, nuMax):
    deltaT = 1/(2*nuMax)
    nu = np.arange(-nuMax,nuMax-(2*nuMax/nPoints),2*nuMax/nPoints)
    t = np.arange(-nPoints/2*deltaT,(nPoints/2-1)*deltaT,deltaT)
    return [nu, t]
def FWHM (x, y):
    ymax = np.max(y)
    #print('\n', ymax*0.5 , 'max*0.5')
    for i, value in enumerate (y):
        if value == ymax:
            index_max = i
   #         print(index_max)
    for i, value in enumerate ( y[0:index_max] ):
        if value >= 0.5*y[index_max]:
            low = i
            break
    #print(y[low] , 'y[low]')
    for i, value in enumerate (y):
        if i > index_max:
            if value <= 0.5*y[index_max]:
                high = i
                break    
   # print(y[high] , 'y[high]', '\n')
    return (x[high] - x[low])






#pulse duration and central wevelegnth :
tau_0 = 35 #input pulse duration FWHM in fs    3.6fs initially
lambda_0 = 800 #input carrier wavelength in nm
omega_0 = 2*np.pi*299792458 / lambda_0 *1e-6  # in PHz
save = 0

#spectral phase :
k00 = 0 # constant phase = CEP   
k10 = 0 # delay in fs
k20 = 275 # input pulse GDD fs^2 (1mm of fused silica = 36.163 fs² at 800nm, cf refractive index.info)
k30 = 0# input pulse TOD fs^3
k40 = 0 #FOD



#fourier transform configuration
nPoints = 2**14   #2**13
[nu,t] = ftAxis(nPoints,4) #%M. Joffre's ftAxis(nPoints, nuMax) % nu up to 4 PHz
                            #%(to sample well the field oscillations, just so you can make pretty plots)
                            #%(and to resolve steep pulsefronts due to SS)    
omega = (2*np.pi)*nu     

#champ E(w)
a = tau_0/(2*np.sqrt(2*np.log(2))) #sigma de la gaussienne
phase_w = k00 + k10*(omega-omega_0) + k20/2*(omega-omega_0)**2 + k30/6*(omega-omega_0)**3 + k40/24*(omega-omega_0)**4 
field_w = np.abs(np.exp((-a**2)/2*(omega-omega_0)**2)) * np.exp(-1j*phase_w)
#field_w = np.abs(np.exp(-(omega-omega_0)**2/(2*a**2)) * np.exp(-1j*phase_w))



#FIGURE
#Frequency domain :
s=13
fig=plt.figure(figsize=(13,4))
gs = gridspec.GridSpec(1, 2, wspace = 0.3)
Ew = plt.subplot(gs[0,0])
Ew.plot(omega, (np.abs(field_w))**2, color = 'navy',label = '|E($\omega$)|')
Ew.set_xlabel('$\omega$ (Phz)', size=s)
Ew.set_ylabel('|E($\omega$)|', color='navy', size=s)
Ew.legend(loc='upper left', fontsize=s)
Ew.tick_params(axis='y', labelcolor='navy',color='navy',labelsize=s)
Ew.tick_params(axis='x', labelsize=s)
Ew.set_title('Frequency domain',fontsize=s+2)

Phiw = Ew.twinx()
Phiw.plot(omega, phase_w, '--', color = 'red', label='$\Phi (\omega)$')
Phiw.set_ylabel('$\Phi (\omega)$ (rad)', color = 'red', size=s)
Phiw.tick_params(axis='y', labelcolor='red',color='red', labelsize=s)
for i, value in enumerate(np.abs(field_w)):
    if value > 1e-15:
        lowi = i
        break
for i, value in enumerate(np.abs(field_w)):
    if value == np.max(np.abs(field_w)):
        highi = i
        break
Phiw.set_xlim([omega[lowi], omega[2*highi - lowi]])
Phiw.set_ylim([np.min(phase_w[lowi : 2*highi-lowi])-5,np.max(phase_w[lowi : 2*highi-lowi])+5])
Phiw.legend(loc = 'upper right', fontsize=s)


#E(t)
field_t = ift(field_w)
Time = plt.subplot(gs[0,1])
plt.plot(t, np.abs(field_t), color='navy',label='|E(t)|')
plt.yticks([])
plt.title('Time domain', fontsize=s+2)
plt.tick_params(axis='both',labelsize=s)
Time = plt.plot(t, field_t, color='grey', label='E(t) = $FT^{-1}$(E($\omega$))')
for i, value in enumerate (np.abs(field_t)):
    if value > 1e-5:
        lowi2 = i 
        break
for i, value in enumerate (np.abs(field_t)):
    if value == np.max(np.abs(field_t)):
        maxxi = i
        break
for i, value in enumerate (np.abs(field_t)):
    if i > maxxi:
        if value < 9.9*1e-6:
            maxi = i 
            break    
plt.xlim([t[lowi2],t[maxi]])
plt.xlabel('t(fs)', size=s)
plt.legend(loc='lower right', fontsize=s)
plt.text(FWHM(t, np.abs(field_t)),np.max(np.abs(field_t))/2, '|E(t)|² is ' + str(round(FWHM(t, (np.abs(field_t))**2),1)) + ' fs FWHM')




if save ==1:
    fig.savefig(r'C:\Users\ouille\Desktop\\' + 'freq_time_' +  str(tau_0) + 'fsFWHM_phase' +  str(round(k00,2)) + '_' + str(k10) +'fs_' + str(k20) +'fs2_' + str(k30) + 'fs3_'+ str(k40) + 'fs4.png', bbox_inches='tight', dpi = 300)



# would be nice to have a spectrogram but it doesnt work yet
#plt.figure()
#ff,tt,Sxx = signal.spectrogram(field_t, fs=int(1/(t[4]-t[3])),window = np.hanning(2**3), nperseg = 2**3, noverlap = 2**2, nfft = 2**5)
#plt.pcolormesh(tt,ff, np.log10(Sxx/np.max(Sxx)), cmap='jet')
#plt.xlim([400,600])
#plt.ylim([-4,0])