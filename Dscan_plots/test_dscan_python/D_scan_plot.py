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
spectrum_wl=raw_data[:,0]   #list of wavelengths
spectrum_I = raw_data[:,1]  #spectral intensity
spectrum_phase = raw_data[:,2]  #spectral phase


#Temporal intensity and phase
raw_data = np.genfromtxt(f_pulse, delimiter=',', dtype=float)  
pulse_T = raw_data[:,0] #list of times
pulse_Env = raw_data[:,1]  #Electric field envelope
pulse_Car = raw_data[:,2]   #Electric field (carrier)


#Retrieved trace
data=  np.genfromtxt(f_retrieved, dtype=float)
data=data[2:np.size(data)]
fheader = np.genfromtxt(f_retrieved, dtype=None)
nWL = float(fheader[0][4:7])  #number of wavelengths
nI= float(fheader[1][3:6])  #number of insertion values
ret_SHG_wl = data[0:nWL] #list of SHG wavelengths(x-axis)
ret_insertions = data[nWL:(nWL+nI)]   #list of insertion values(y-axis)
ret_dscan = np.reshape(data[nWL+nI:np.size(data)],(nI, nWL))   #matrix where each line corresponds to the spectrum for a given insertion


#Measured trace
data = np.genfromtxt(f_measured, dtype=float)
data=data[2:np.size(data)]
fheader = np.genfromtxt(f_retrieved, dtype=None)
nWL = float(fheader[0][4:7]) #number of wavelengths
nI= float(fheader[1][3:6])  #number of insertion values
meas_SHG_wl = data[0:nWL]  #list of SHG wavelengths(x-axis)
meas_insertions = data[nWL:(nWL+nI)]  #list of insertion values(y-axis)
meas_dscan = np.reshape(data[nWL+nI:np.size(data)],(nI, nWL))  #matrix where each line corresponds to the spectrum for a given insertion





##PLOT THE FIGURE:
##Structure of the big figure :
fig = plt.figure(figsize=(9,8))
gs = GridSpec(2, 2, wspace=0.4, hspace=0.4)

#SUBPLOT : Spectral intensity and phase
spectral = plt.subplot(gs[1, 1])
spectral.plot(spectrum_wl, spectrum_I,  color='navy')
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
temp = plt.plot(pulse_T, pulse_Car**2,  color='navy')
temp = plt.xlabel('Time (fs)')
temp = plt.xlim([-40,40])
temp = plt.title('Retrieved temporal shape', weight='bold')
temp = plt.ylabel('Temporal Intensity (a.u.)')
for k in np.arange (0,np.size(pulse_Car)):
    if round((pulse_Car[k])**2,2) >= 0.50:
        k1 = k
        break
for k in np.arange (k1+5,np.size(pulse_Car)):
    if round((pulse_Car[k])**2,2) <= 0.50:
        k2 = k 
        break
FWHM = pulse_T[k2]-pulse_T[k1]
temp = plt.axvspan(pulse_T[k1], pulse_T[k2], facecolor='grey', alpha=0.5)
temp = plt.text(4,0.5, str(np.round(FWHM, decimals = 1)) + ' fs FWHM', color='black')


#SUBPLOT : measured d-scan trace                    
measured = plt.subplot(gs[0,0])
measured = plt.pcolor(meas_SHG_wl, meas_insertions, meas_dscan, cmap='jet')
measured = plt.xlabel('Wavelength (nm)');
measured = plt.ylabel('Wedges insertion (mm)');
measured = plt.title('Measured', weight='bold');

                    
#SUBPLOT : retrieved d-scan trace           
retrieved = plt.subplot(gs[0,1])
retrieved = plt.pcolor(ret_SHG_wl, ret_insertions, ret_dscan, cmap='jet')
retrieved = plt.xlabel('Wavelength (nm)');
retrieved = plt.ylabel('Wedges insertion (mm)');
retrieved = plt.title('Retrieved', weight='bold');
                     
                     
#save figure
plt.savefig(str(filepath) + '\\' + str(filename) + '_dscan_analyzed'+ '.png')




















#### AJOUT TOD
#plt.close()
#plt.figure()
##domaine fréquentiel :
#c0=300 #nm/fs
#def freq_conv(lam):
#    return 2 * np.pi * c0 / lam
#spectrum_freq = freq_conv(spectrum_wl)
##phase_freq = spectrum_phase[::-1]
#plt.plot(spectrum_freq, spectrum_phase)
#plt.xlim([1.88,3.43])
##plt.xlim([2.4,5])
#plt.ylim([-2,3])
#plt.xlabel('w (rad/fs)')
##add TOD:
#TOD=-10 #fs3
#w0=freq_conv(800)
#def add_TOD(w,ph):
#    return ph + TOD/6*((w-w0)**3)
#plt.plot(spectrum_freq, add_TOD(spectrum_freq,spectrum_phase), color='r')
##FFT, retrieve duration
#def E(I,ph):
#    return np.sqrt(I) * np.exp(1j*ph)
#E_w = E(spectrum_I[600:1900], add_TOD(spectrum_freq[600:1900],spectrum_phase[600:1900]))
##E_w, w : domaine spectral
##FT:
#E_t = np.fft.fft(E_w)
#
#
#plt.figure()
#E_t_reconstructed = E_t[np.size(E_t)/2 : np.size(E_t)]
#E_t_reconstructed = np.append (E_t_reconstructed, (E_t[0:np.size(E_t)/2]))
#plt.plot((abs(E_t_reconstructed))**2, color='r')
##1/spectrum_phase[600:1900], 
##plt.xlim([640,660])
#plt.ylim(0,300000)
#plt.xlim([640,660])
##plt.close()
#
##sans correction :
#E_temoin = E (spectrum_I[600:1900], spectrum_phase[600:1900])
#E_temoin_fft = np.fft.fft (E_temoin)
#E_temoin_reconstructed = E_temoin_fft[np.size(E_temoin_fft)/2 : np.size(E_temoin_fft)]
#E_temoin_reconstructed = np.append(E_temoin_reconstructed, E_temoin_fft[0:np.size(E_temoin_fft)/2])
#plt.plot((abs(E_temoin_reconstructed))**2)
##plt.xlim([620,680])

















###deal with fft x axis to get the duration

plt.figure()
#domaine fréquentiel :
c0=300 #nm/fs
def freq_conv(lam):
    return 2 * np.pi * c0 / lam
spectrum_freq = freq_conv(spectrum_wl)  #en fs-1
spec_freq = np.append(spectrum_freq, np.arange(spectrum_freq[np.size(spectrum_freq)-1], 0, round((spectrum_freq[2048] -spectrum_freq[2047]),4)))
spec_phase = np.append(spectrum_phase, np.zeros(np.size(spec_freq)-np.size(spectrum_freq)))
plt.plot(spec_freq, spec_phase)
plt.xlim([1.88,3.43])
#plt.xlim([2.4,5])
plt.ylim([-2,3])
plt.xlabel('w (fs-1)')
plt.ylabel ('Spectral phase (rad)')
#add TOD:
TOD=4 #fs3
w0=freq_conv(800)
def add_TOD(w,ph):
    return ph + TOD/6*((w-w0)**3)
plt.plot(spec_freq, add_TOD(spec_freq,spec_phase), color='r')


#check phi(lambda)
#plt.figure()
#plt.plot(spectrum_wl, spectrum_phase)
#plt.xlim([550,1000])
#plt.ylim([-2,3])
#plt.xlabel('Wavelength (nm)')
#plt.ylabel ('Spectral phase (rad)')

#Electric field in spectral domain 
def E(I,ph):
    return np.sqrt(I) * np.exp(1j*ph)
spec_I = np.append(spectrum_I, np.zeros(np.size(spec_freq)-np.size(spectrum_freq)))
E_w = E(spec_I, add_TOD(spec_freq,spec_phase))


#échantilloner E_w
x=np.arange(0,spec_freq[1],0.0001)
E_w[0]=0

#   #PROBLEME ECHANTILLONAGE !!
#Ew_ech =  np.interp(x,spec_freq[::-1], abs(E_w[::-1]))
#plt.figure()
#plt.plot(spec_freq, abs(E_w))
#plt.plot(x, Ew_ech, color='r')
#plt.xlabel('w (fs-1)')
#plt.ylabel ( 'I(w)  rouge : echantilloné / bleu : initial' )
#plt.xlim([1.5,4.5])
##plt.xlim([-2,10])
#
#plt.figure()
#plt.plot(spectrum_wl, spectrum_I)
#plt.xlabel('wavelength (nm)')
#plt.xlim([400,1100])
#plt.ylabel('spectral intensity')
#
##FFT :
#Et_ech = np.fft.fft(Ew_ech)
#
##x-axis
#plt.figure()
#Et_rec = Et_ech[np.size(Et_ech)/2 : np.size(Et_ech)]
#Et_rec = np.append (Et_rec, (Et_ech[0:np.size(Et_ech)/2]))
#time_ax = np.arange(-5000,5000, 10000/np.size(x))
#plt.plot(time_ax, (abs(Et_rec[0:np.size(Et_rec)]))**2, color='r')
#plt.xlim([-5,5])
#plt.xlabel('time(fs)')
#
#
###sans correction :
#E_temoin = E (spec_I, spec_phase)
#E_temoin_ech = np.interp(x,spec_freq[::-1], abs(E_temoin[::-1]))
#E_temoin_fft = np.fft.fft (E_temoin_ech)
#E_temoin_reconstructed = E_temoin_fft[np.size(E_temoin_fft)/2 : np.size(E_temoin_fft)]
#E_temoin_reconstructed = np.append(E_temoin_reconstructed, E_temoin_fft[0:np.size(E_temoin_fft)/2])
#plt.plot(time_ax,(abs(E_temoin_reconstructed))**2)
#
#
#plt.figure()
#plt.plot((abs(np.fft.fft(E_w)[0:100]))**2)
#plt.plot((abs(np.fft.fft(E_temoin)[0:100]))**2)
##plt.plot(E_temoin)
#


###########################TENTATIVE NUMERO 2
   #PROBLEME ECHANTILLONAGE !! ça devient réel !! 
#   ### et problème avec l'échelle de temps (un facteur 10 on dirait)
Ew_real = np.zeros(np.size(E_w))
for i, value in enumerate (E_w):
    Ew_real[i]=value.real

Ew_ech =  np.interp(x,spec_freq[::-1], Ew_real[::-1])

#FFT :
Et_ech = np.fft.fft(Ew_ech)

#x-axis
plt.figure()
Et_rec = Et_ech[np.size(Et_ech)/2 : np.size(Et_ech)]
Et_rec = np.append (Et_rec, (Et_ech[0:np.size(Et_ech)/2]))
time_ax = np.arange(-5000,5000, 10000/np.size(x))
plt.plot(time_ax, (abs(Et_rec[0:np.size(Et_rec)]))**2, color='r')
plt.xlim([-5,5])
plt.xlabel('time(fs)')

#
###sans correction :
E_temoin = E (spec_I, spec_phase)
Etem_real = np.zeros(np.size(E_w))
for i, value in enumerate (E_temoin):
    Etem_real[i]=value.real
E_temoin_ech = np.interp(x,spec_freq[::-1], (Etem_real[::-1]))
E_temoin_fft = np.fft.fft (E_temoin_ech)
E_temoin_reconstructed = E_temoin_fft[np.size(E_temoin_fft)/2 : np.size(E_temoin_fft)]
E_temoin_reconstructed = np.append(E_temoin_reconstructed, E_temoin_fft[0:np.size(E_temoin_fft)/2])
plt.plot(time_ax,(abs(E_temoin_reconstructed))**2)

#
#plt.figure()
#plt.plot((abs(np.fft.fft(E_w)[0:100]))**2)
#plt.plot((abs(np.fft.fft(E_temoin)[0:100]))**2)
#plt.plot(E_temoin)























########ATTEMP TO RETRIEVE THE TOD :
    
######################SPECTRAL_PHASE_POLYNOMIAL_FIT
########polynomial fit
#x_polyfit=np.arange(0,280)
#indices = np.arange(0,280)
#y_polyfit=np.arange(0,280)
#for i, value in enumerate(x_polyfit):
#    x_polyfit[i]=500+i*2
#    for k, valeur in enumerate(spectrum_wl):
#        if round(valeur) == x_polyfit[i]:
#            indices[i]=k
#            y_polyfit[i]=spectrum_phase[k]
#

#
#my_fit = np.polyfit(x_polyfit,y_polyfit, 7)
#my_fit_func = np.poly1d(my_fit)
#########plot
#plt.plot(spectrum_wl, my_fit_func(spectrum_wl))
#plt.plot(spectrum_wl, spectrum_phase, color='r')
#plt.xlim([400,1100])
#plt.ylim([-3,4])
#
















############ATTEMPT : Fit with a Taylor series expansion, order = 3
#plt.close()
#inf = 550
#sup = 1000
#step = 4
#x_polyfit = np.arange(inf, sup, step)
#indices = np.arange(inf, sup, step)
#y_polyfit = np.arange(inf, sup, step)
#
#for i, value in enumerate(x_polyfit):
#    for k, valeur in enumerate(spectrum_wl):
#        if value == round(valeur):
#            indices[i] = k
#            y_polyfit[i] = spectrum_phase[k]
#l_0 = 800   #800nm central wavelength
#c0=300   #speed of light (nm/fs)
#
#def T (l,a,b,c,d):
#    return a + b*2*np.pi*c0*(1/l-1/l_0) + c/2*((2*np.pi*c0)**2)*(1/l-1/l_0)**2 + d/6*((2*np.pi*c0)**3)*(1/l-1/l_0)**3 
#
#initial_guesses=[-1.93764783e-02,  0,   0, -100]
#func_fit = sc.optimize.curve_fit(T, x_polyfit, y_polyfit, initial_guesses)
##func_fit = sc.optimize.curve_fit(T, x_polyfit, y_polyfit)
#
#
#a=func_fit[0][0]
#b=func_fit[0][1]
#c=func_fit[0][2]
#d=func_fit[0][3]
#
#plt.plot(spectrum_wl, T(spectrum_wl,a,b,c,d))
#plt.plot(spectrum_wl, spectrum_phase)
#plt.ylim([-2,2.5])
#plt.xlim([inf,sup])
#
#print('GDD = ' + str(func_fit[0][2]) +' fs2 and TOD = ' + str(func_fit[0][3]) + ' fs3')
#
#











#
##############ATTEMPT : Fit with a Taylor series expansion, order = 4

inf = 520
sup = 1000
step = 1
x_polyfit = np.arange(inf, sup, step)
indices = np.arange(inf, sup, step)
y_polyfit = np.arange(inf, sup, step)

for i, value in enumerate(x_polyfit):
    for k, valeur in enumerate(spectrum_wl):
        if value == round(valeur):
            indices[i] = k
            y_polyfit[i] = spectrum_phase[k]
l_0 = 800   #800nm central wavelength
c0=300   #speed of light (nm/fs)

def T (l,a,b,c,d,e):
    return a + b*2*np.pi*c0*(1/l-1/l_0) + c/2*((2*np.pi*c0)**2)*(1/l-1/l_0)**2 + d/6*((2*np.pi*c0)**3)*(1/l-1/l_0)**3 +  e/24*((2*np.pi*c0)**4)*(1/l-1/l_0)**4 

initial_guesses=[-1.93764783e-02,  -6.25325242e-01,   4.02739828e-01, -10.12931757e+01,  0]
func_fit = sc.optimize.curve_fit(T, x_polyfit, y_polyfit, initial_guesses)
#func_fit = sc.optimize.curve_fit(T, x_polyfit, y_polyfit)


a=func_fit[0][0]
b=func_fit[0][1]
c=func_fit[0][2]
d=func_fit[0][3]
e=func_fit[0][4]
plt.figure()
plt.plot(spectrum_wl, T(spectrum_wl,a,b,c,d,e))
plt.plot(spectrum_wl, spectrum_phase)
plt.ylim([-2,2.5])
plt.xlim([inf,sup])

print('GDD = ' + str(func_fit[0][2]) +' fs2 and TOD = ' + str(func_fit[0][3]) + ' fs3')




















#############ATTEMPT : Fit with a Taylor series expansion, order=6
#plt.close()
#inf = 520
#sup = 1000
#step = 1
#x_polyfit = np.arange(inf, sup, step)
#indices = np.arange(inf, sup, step)
#y_polyfit = np.arange(inf, sup, step)
#
#for i, value in enumerate(x_polyfit):
#    for k, valeur in enumerate(spectrum_wl):
#        if value == round(valeur):
#            indices[i] = k
#            y_polyfit[i] = spectrum_phase[k]
#l_0 = 800   #800nm central wavelength
#c0=300   #speed of light (nm/fs)
#
#
#
#def T (l,a,b,c,d,e,f,g):
#    return a + b*2*np.pi*c0*(1/l-1/l_0) + c/2*((2*np.pi*c0)**2)*(1/l-1/l_0)**2 + d/6*((2*np.pi*c0)**3)*(1/l-1/l_0)**3 +  e/24*((2*np.pi*c0)**4)*(1/l-1/l_0)**4 + f/120*((2*np.pi*c0)**5)*(1/l-1/l_0)**5  + g/720*((2*np.pi*c0)**6)*(1/l-1/l_0)**6 
#
#initial_guesses=[-1.00438667e-01,   1.92438992e-01,   6.22149225e+00,
#         -3.68560826e+01,  -2.18716631e+02,   0,
#         0]
#func_fit = sc.optimize.curve_fit(T, x_polyfit, y_polyfit, initial_guesses)
##func_fit = sc.optimize.curve_fit(T, x_polyfit, y_polyfit)
#
#
#a=func_fit[0][0]
#b=func_fit[0][1]
#c=func_fit[0][2]
#d=func_fit[0][3]
#e=func_fit[0][4]
#f=func_fit[0][5]
#g=func_fit[0][6]
#
#plt.plot(spectrum_wl, T(spectrum_wl,a,b,c,d,e,f,g))
#plt.plot(spectrum_wl, spectrum_phase)
#plt.ylim([-2,2.5])
#plt.xlim([inf,sup])
#
#print('GDD = ' + str(func_fit[0][2]) +' fs2 and TOD = ' + str(func_fit[0][3]) + ' fs3')
#
#
#[xf,yf]=plt.ginput(10)
#plt.show()
#
