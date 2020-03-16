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
from matplotlib.gridspec import GridSpec






#%% input pulse parameters
c = 3e8
tau_0 = 4 #input pulse duration FWHM in fs    3.6fs initially
lambda_0 = 790 #input carrier wavelength in nm
omega_0 = 2*np.pi*299792458 / lambda_0 *1e-6  # in PHz



#spectral phase :
k00 = 0.0*np.pi# constant phase = CEP   
k10 = 30 # delay in fs
k20 = 0# input pulse GDD fs^2 (1mm of fused silica = 36.163 fs² at 800nm, cf refractive index.info)
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

#frequency doubling
#twice the spectral phase at w/2
omega_d = 2*np.pi*299792458 / 480 *1e-6  #enveloppe centered around 480nm
#spectral bandwidth depends on the BBO efficiency (thickness, orientation)... parameter a(t0) to be tweaked below... width of the gaussian
att =0.08 #attenuation induced by rotating the PBS
double = np.sqrt(att) * np.exp((-(a(10))**2)/2*(omega-omega_d)**2) * np.exp(-2j*phase(omega/2,omega_0))

# interferences
interf = (np.abs(simple + double))**2

     

#%% PLOTS
############## spectrum seen on the APS800 or the oscilloscope of the fringeezz         
fig1 = plt.figure(figsize=(10,7))
plt.plot((omega), (np.abs(simple))**2, color = 'red', label = '$I_{fund}$') 
plt.plot((omega), (np.abs(double))**2 , color= '#00a6ff' , label = '$I_{SHG}$') 
plt.plot((omega), (np.abs(double))**2   + (np.abs(simple))**2 , '--', color= '#3c9900' , alpha = 0.3, label = '$I_{fund} + I_{SHG}$') 
plt.plot((omega), interf , color = '#3c9900', label = '|$E_{SHG}+E_{fund}$|²') 
plt.xlim(3.1,4.3)
plt.ylim(0,0.6)
plt.xticks(fontsize=14) ; plt.yticks(fontsize=14)
plt.xlabel('Frequency $\omega$ (PHz)', fontsize=14)
plt.ylabel('Spectral intensity (a.u.)', fontsize=14)
plt.legend(fontsize='14')
plt.axvline(x=3.73, color='purple')
plt.axvline(x=4.14, color='purple')
fig1.suptitle('CEP = ' + str(round(k00/np.pi,1)) + '$\pi$', fontsize = 16)
#fig1 = plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\CEP_fringes\GDD100fs2\\' + 'GDD' +str(k20)+   'fs2_CEP' + str(round(k00/np.pi,1))[0] + '-' + str(round(k00/np.pi,1))[2] + 'pi_480nm.png')




############## visualize the spectral phases and the difference
fig2 = plt.figure(figsize=(10,6))
plt.plot(omega, phase(omega, omega_0), c = 'red', label = '$\phi_{fund} (\omega)$')
plt.plot(omega, 2*phase(omega/2, omega_0), c = '#00a6ff', label = '$\phi_{SHG} (\omega) = 2\phi_{fund} (\omega/2)$')
plt.xlim(0,6)
plt.ylim(-100,50)
#plt.axvline(omega_0, color='black', alpha = 0.2)
#plt.axvline(2*omega_0, color='black', alpha = 0.2)
delta_phi = phase(omega, omega_0) -2*phase(omega/2, omega_0)
plt.plot( omega,  delta_phi, color='#3c9900'  , label = '$\Delta \phi = \phi_{fund} - \phi_{SHG}$')
plt.axvline(x=3.73, color='purple', alpha = 0.3,label ='Fringeezz')
plt.axvline(x=4.14, color='purple', alpha = 0.3)
plt.legend(fontsize='13')
plt.xlabel('frequency (PHz)')
plt.ylabel('phase(rad)')
plt.title( 'CEP=' + str(round(k00/np.pi,1)) + '$\pi$, GD=' + str(k10) +'fs, GDD=' + str(k20) + 'fs², TOD=' + str(k30) + '$fs^3$, FOD=' + str(k40) + '$fs^4$' )  



################### cos(Delta phi)
fig3 = plt.figure(figsize=(6,4))
plt.plot(omega, np.cos(delta_phi), color='#3c9900')
plt.xlabel('frequency $\omega$')
plt.ylabel('$cos(\Delta\phi)$')
plt.xlim(3.1,4.3)
plt.ylim(-1.3,1.3)
plt.axvline(x=3.73, color='purple', alpha = 0.6, label ='Fringeezz')
plt.axvline(x=4.14, color='purple', alpha = 0.6)
plt.legend()



















#%% D-scan spectrum and spectral phase

#Load data
filepath = r'Z:\Laser\D-Scan\HCF\20191217'
filename= 'PC1332_1300mbar_4_45mJ_3_18fs'


f_spectrum  =  filepath + '\\' + filename + '_spectrum.csv'
f_pulse     =  filepath +  '\\' + filename + '_retrieved_pulse.csv'
f_retrieved = filepath +  '\\' + filename + '_retrieved_dscan.csv'
f_measured  = filepath +  '\\' + filename + '_measured_dscan.csv'


#Spectral intensity and phase
raw_data = np.genfromtxt(f_spectrum, delimiter=',', dtype=float)
spectrum_wl=raw_data[:,0]   #list of wavelengths
spectrum_I = raw_data[:,1]  #spectral intensity
spectrum_phase = raw_data[:,3]  #spectral phase


#Temporal intensity and phase
raw_data = np.genfromtxt(f_pulse, delimiter=',', dtype=float)  
pulse_T = raw_data[:,0] #list of times
pulse_Env = raw_data[:,3]  #Electric field envelope
pulse_Car = raw_data[:,4]   #Electric field (carrier)


#add spectral phase
c0=299.792458 #nm/fs
lowwl=450
highwl=1050
nu2 =  c0/lowwl 
nu1 =  c0/highwl 
####### interpolate
end= np.size(nu)-1
spectrum_I_nu=[]
spectrum_I_nu = np.interp(nu[nPoints/2 + 1 : end], c0/spectrum_wl[::-1] , spectrum_I[::-1]*(spectrum_wl[::-1]**2))    #I(w) interpolé à partir de I(lambda)
spectrum_phase_nu = np.interp( nu[nPoints/2 +1:end], c0/spectrum_wl[::-1], spectrum_phase[::-1])           #Phi(w) interpolée à partir de Phi(lambda)
##########
lambdac = 780. #nm
nuc = c0/lambdac  # fs-1
CEP_d = 1.8*np.pi
GD= 0.
GDD= 0.
TOD= 0.
FOD= 0.
correctionphase =CEP_d + GD*2*np.pi*(nu[nPoints/2 + 1:end]-nuc) + GDD/2*(2*np.pi*(nu[nPoints/2 +1:end]-nuc))**2 + TOD/6*(2*np.pi*(nu[nPoints/2 +1:end]-nuc))**3 + FOD/24*(2*np.pi*(nu[nPoints/2 +1:end]-nuc))**4
################



#### put together symmetric spectrum over both negative and positive frequencies
spectrum_I_nu = np.append(np.zeros(nPoints/2),spectrum_I_nu)
spectrum_phase_nu = np.append(np.zeros(nPoints/2),spectrum_phase_nu)
spectrum_phase_nu_corr = spectrum_phase_nu - np.append(np.zeros(nPoints/2),correctionphase)

spectrum_I_nu = np.nan_to_num(spectrum_I_nu)
spectrum_phase_nu = np.nan_to_num(spectrum_phase_nu)
pulse_t = ft( np.sqrt(spectrum_I_nu)*np.exp((-1j)*(spectrum_phase_nu-np.pi/2)) )
pulse_I = (pulse_t *np.conj(pulse_t)).real
          
maxI = max(pulse_I)
maxind = 2039
pulse_tau = np.diff(t[abs(np.diff(np.sign(pulse_I - 0.5*(max(pulse_I)))))==2])
spectrum_phase_nu_corr = np.nan_to_num(spectrum_phase_nu_corr)

pulse_t_corr = ft( np.sqrt(spectrum_I_nu)*np.exp(-1j*(spectrum_phase_nu_corr-np.pi/2)) )
pulse_I_corr = (pulse_t_corr *np.conj(pulse_t_corr)).real
pulse_tau_corr = np.diff(t[abs(np.diff(np.sign(pulse_I_corr - 0.5*(max(pulse_I_corr)))))==2])

pulse_t_FTL = ft( np.sqrt(spectrum_I_nu)*np.exp(-1j*(-np.pi/2)) )
pulse_I_FTL = (pulse_t_FTL *np.conj(pulse_t_FTL)).real
pulse_tau_FTL = np.diff(t[abs(np.diff(np.sign(pulse_I_FTL - 0.5*(max(pulse_I_FTL)))))==2])























#%% CEP fringes with dscan meas
fig4 = plt.figure(figsize=(10,7))

#fundmanetal spectrum
SpecNu = fig4.add_subplot(111)
SpecNu.set_xlim(3.6,2*np.pi*nu2)
SpecNu.set_ylim(ymin=0, ymax=0.15)
SpecNu.set_xlabel('Frequency $\omega$ (PHz)', fontsize=14)
SpecNu.set_ylabel('Spectral intensity (a.u.)', fontsize=14)
SpecNu.tick_params(labelsize=14) 
fund_d = np.sqrt(spectrum_I_nu)/np.max(np.sqrt(spectrum_I_nu))   *np.exp(-1j*(spectrum_phase_nu_corr)) 
SpecNu.plot(2*np.pi*nu[1:-1] ,   (np.conj(fund_d) * fund_d).real  , color='red'   ,  label = '$I_{fund}$')

# SHG
att_d = 0.04
double_d_phase=np.zeros(np.size(spectrum_phase_nu_corr))
i1 = round(np.size(nu)/2 )
SHG_ws = np.arange(3.59872, 4.18163, (4.18163-3.59872)/(9555-9365) )  #double_d_phase entre les indices 9171 et 10148 doit contenir 2*spectrum_phase_nu_corr entre les indices 8681 et 9171
#plt.figure()
#plt.plot(2*np.pi*nu[1:-1], spectrum_phase_nu_corr)
#plt.xlim(1.8, 4.18)
#plt.ylim(-15,5)
SHG_phase_interp = np.interp(SHG_ws/2, 2*np.pi*nu[1:-1], 2*spectrum_phase_nu_corr)
#plt.plot(SHG_ws, SHG_phase_interp, '--')
double_d_phase = np.append(   np.append( ( np.zeros(9365) ) , SHG_phase_interp)    , np.zeros( np.size(nu) - 9555)  )
#plt.plot( 2*np.pi*nu,  double_d_phase)
double_d = np.sqrt(att_d) * np.exp((-(a(10))**2)/2*(2*np.pi*nu[i1:]-omega_d)**2) *np.exp(-1j*double_d_phase[i1:] ) 
SpecNu.plot(2*np.pi*nu[i1:], (np.abs(double_d))**2 , color= '#00a6ff'  , label = '$I_{SHG}$')

#sum
sum_d =  (np.conj(fund_d[i1-2:]) * fund_d[i1-2:]).real +  (np.abs(double_d))**2
SpecNu.plot(2*np.pi*nu[i1:]  ,  sum_d,  '--', color= '#3c9900' , alpha = 0.3, label = '$I_{fund} + I_{SHG}$')

#Interf
interf_d = fund_d[i1-2:] + double_d
SpecNu.plot(2*np.pi*nu[i1:]   ,  (np.conj(interf_d) * interf_d).real    ,  color = '#3c9900', label = '|$E_{SHG}+E_{fund}$|²'  )
SpecNu.legend( fontsize = 14)
SpecNu.axvline(x=3.73, color='purple')
SpecNu.axvline(x=4.14, color='purple')

fig4.suptitle('CEP = ' + str( round(CEP_d/np.pi,1) ) + '$\pi$' , fontsize = 16 )
fig4.text(0.15,-0.02,'additional phase : GD=' +str(GD) + 'fs, GDD=' + str(GDD) +'fs², TOD=' + str(TOD) +'fs^3, FOD=' + str(FOD) + 'fs^4'    , fontsize=14)
#fig4 = plt.savefig(     r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\CEP_fringes\dscan_GDD100fs2\\' + 'GDD' +str(GDD)+   'fs2_CEP' + str   (round(CEP_d/np.pi,1)) + 'pi_480nm.png' ,  bbox_inches='tight')





















#%% SHG, fund and Delta Phi
fig5 = plt.figure(figsize=(10,7))

delta_phase_dscan = spectrum_phase_nu_corr[i1-2:]    -  double_d_phase[i1:]
delta = fig5.add_subplot(111)
delta.plot(2*np.pi*nu[i1:], delta_phase_dscan, color= '#3c9900', label = '$\phi_{fund} (\omega)$' )
delta.plot(2*np.pi*nu[i1:], spectrum_phase_nu_corr[i1-2:]   , color = 'red' , label='$\phi_{fund} - \phi_{SHG}$'  )
delta.plot(2*np.pi*nu[i1:], double_d_phase[i1:]  , color = '#00a6ff'   ,  label = '$\phi_{SHG} (\omega) = 2\phi_{fund} (\omega/2)$')
delta.set_xlim(2*np.pi*nu1,2*np.pi*nu2)
delta.set_ylim(-145,4)
delta.legend(fontsize=14)
delta.axvline(x=3.73, color='purple', alpha = 0.3)
delta.axvline(x=4.14, color='purple', alpha = 0.3)

fig5.text(0.15,-0.02,'additional phase : GD=' +str(GD) + 'fs, GDD=' + str(GDD) +'fs², TOD=' + str(TOD) +'fs^3, FOD=' + str(FOD) + 'fs^4'    , fontsize=14)




























#%% spectrum and spectral phase from dscan + correction
fig6 = plt.figure(figsize=(10,7))

Phase_nu = fig6.add_subplot(111)
Phase_nu.plot(2*np.pi*nu[1:-1],-spectrum_phase_nu_corr, color='red' )
Phase_nu.set_ylim(-10,160)
Phase_nu.set_ylabel('spectral phase', color='red', fontsize=14)
Phase_nu.set_xlabel('Frequency $\omega$ (PHz)',fontsize=14)
Phase_nu.tick_params( axis = 'both', labelsize=14)
Phase_nu.tick_params( axis='y',  labelcolor='red')

Int_nu = Phase_nu.twinx()
Int_nu.plot(2*np.pi*nu[1:-1] ,   (np.conj(fund_d) * fund_d).real , color='navy')
Int_nu.set_xlim(2*np.pi*nu1, 2*np.pi*nu2)
Int_nu.set_ylim(0,1)
Int_nu.set_ylabel('Spectral intensity (a.u.)', color='navy',fontsize=14)
Int_nu.tick_params(labelsize=14, labelcolor='navy')
Int_nu.axvline(x=3.73,  ls = '--', color='purple', alpha= 0.6)
Int_nu.axvline(x=4.14,  ls = '--', color='purple', alpha = 0.6)

fig6.suptitle('d-scan measurement, additional phase : \n GD=' +str(GD) + 'fs, GDD=' + str(GDD) +'fs², TOD=' + str(TOD) +'fs^3, FOD=' + str(FOD) + 'fs^4'    , fontsize=16)