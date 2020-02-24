# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 18:18:43 2020

@author: ouille


modélisation franges CEP

@ faire :
    - utiliser les spectres mesurés
    - utiliser la phase spectrale mesurée par le dscan
    - prendre en compte l'efficacité spectrale du BBO selon l'épaisseur
    



"""

#%% libraries
import numpy as np
from matplotlib import pyplot as plt
import sys
if not r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' in sys.path : 
    sys.path.insert(0,r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' )
from FourierTransforms import ft, ift, ftAxis
from FWHM import FWHM






#%% input pulse parameters
c = 3e8
tau_0 = 4 #input pulse duration FWHM in fs    3.6fs initially
lambda_0 = 790 #input carrier wavelength in nm
omega_0 = 2*np.pi*299792458 / lambda_0 *1e-6  # in PHz



#spectral phase :
k00 = 0.5*np.pi# constant phase = CEP   
k10 = 0 # delay in fs
k20 = 36# input pulse GDD fs^2 (1mm of fused silica = 36.163 fs² at 800nm, cf refractive index.info)
k30 = 0# input pulse TOD fs^3
k40 = 0 #FOD

#fourier transform configuration
nPoints = 2**14   #2**13
[nu,t] = ftAxis(nPoints,4) #%M. Joffre's ftAxis(nPoints, nuMax) % nu up to 4 PHz
                            #%(to sample well the field oscillations, just so you can make pretty plots)
                            #%(and to resolve steep pulsefronts due to SS)    
omega = (2*np.pi)*nu     



#%% spectral domain        
#field E(w)
def a(t0):   #sigma de la gaussienne pour une durée t0 (en fs)
    return t0/(2*np.sqrt(2*np.log(2)))
def phase(w,w0):
    return k00 + k10*(w-w0) + k20/2*(w-w0)**2 + k30/6*(w-w0)**3 + k40/24*(w-w0)**4 
def field(w,w0,t0):
    return np.abs(np.exp((-(a(t0))**2)/2*(w-w0)**2)) * np.exp(-1j*phase(w,w0))



#plt.figure()
#plt.title('fundamental spectrum')
#plt.plot(2*np.pi*c/(omega*1e6), (np.abs(field(omega,omega_0,4)))**2)    #x en nm
#plt.xlim(600,1200)

#plt.figure()
#plt.title('initial spectral phase')
#plt.plot(2*np.pi*c/(omega*1e6), phase(omega,omega_0) , color = 'red') 
#plt.xlim(600,1200)
#plt.ylim(-10,50)


#fundamental spectrum, 3fs duration, centered around 790nm
simple = field(omega,omega_0,3)
#simple = np.exp((-(a(3))**2)/2*(omega-omega_0)**2) * np.exp(-1j*phase(1*omega,omega_0))


#%% frequency doubling
#twice the spectral phase
omega_d = 2*np.pi*299792458 / 480 *1e-6  #enveloppe centered around 480nm
#spectral bandwidth depends on the BBO efficiency (thickness, orientation)... parameter a(t0) to be tweaked below... width of the gaussian
att =0.08 #attenuation induced by rotating the PBS
double = np.sqrt(att) * np.exp((-(a(10))**2)/2*(omega-omega_d)**2) * np.exp(-2j*phase(omega,omega_0))



#%% interferences
interf = np.abs((simple + double)**2)

     
         
#fig = plt.figure(figsize=(10,7))
#plt.plot(2*np.pi*c/(omega*1e6), (np.abs(simple))**2, color = 'red', label = '$I_{fund}$') 
#plt.plot(2*np.pi*c/(omega*1e6), (np.abs(double))**2 , color= '#00a6ff' , label = '$I_{SHG}$') 
#plt.plot(2*np.pi*c/(omega*1e6), (np.abs(double))**2   + (np.abs(simple))**2 , '--', color= '#3c9900' , alpha = 0.3, label = '$I_{fund} + I_{SHG}$') 
#plt.plot(2*np.pi*c/(omega*1e6), interf , color = '#3c9900', label = '|$E_{SHG}+E_{fund}$|²') 
#plt.xlim(425,650)
#plt.ylim(0,0.3)
#plt.xlabel('Wavelength (nm)', fontsize=14)
#plt.ylabel('Spectral intensity (a.u.)', fontsize=14)
#plt.legend(fontsize='14')

         
fig = plt.figure(figsize=(10,7))
plt.plot((omega), (np.abs(simple))**2, color = 'red', label = '$I_{fund}$') 
plt.plot((omega), (np.abs(double))**2 , color= '#00a6ff' , label = '$I_{SHG}$') 
plt.plot((omega), (np.abs(double))**2   + (np.abs(simple))**2 , '--', color= '#3c9900' , alpha = 0.3, label = '$I_{fund} + I_{SHG}$') 
plt.plot((omega), interf , color = '#3c9900', label = '|$E_{SHG}+E_{fund}$|²') 
plt.xlim(omega_0,4.5)
plt.ylim(0,1)
plt.xticks(fontsize=14) ; plt.yticks(fontsize=14)
plt.xlabel('Frequency $\omega$ (PHz)', fontsize=14)
plt.ylabel('Spectral intensity (a.u.)', fontsize=14)
plt.legend(fontsize='14')





fig.suptitle('CEP = ' + str(round(k00/np.pi,1)) + '$\pi$', fontsize = 16)
fig = plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\CEP_fringes\GDD36fs2\\' + 'GDD' +str(k20)+   'fs2_CEP' + str(round(k00/np.pi,1))[0] + '-' + str(round(k00/np.pi,1))[2] + 'pi_480nm.png')