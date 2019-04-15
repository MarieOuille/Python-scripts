# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:28:57 2019

@author: ouille
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


# 0 ) Je définis la transformation de Fourier :
def ft(data):
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(data)))
def ftAxis(nPoints, nuMax):
    deltaNu = 2*nuMax/nPoints
    deltaT = 1/(2*nuMax)
    nu = np.arange(-nuMax,nuMax-(2*nuMax/nPoints),2*nuMax/nPoints)
    t = np.arange(-nPoints/2*deltaT,(nPoints/2-1)*deltaT,deltaT)
    return [nu, t]
nPoints = 2**18
c0=299.792458 #nm/fs
lambdac = 800. #nm
nuc = c0/lambdac  # fs-1




# 1) je définis mon impulsion gaussienne
n = 60  #nombre de bounces sur miroirs TOD 
dw = 0.03
GD= -7900 ; GDD=-2500; TOD=600000 ; FOD=0   ## n=60 design 1901
#GD= -11110 ; GDD=0; TOD=600000 ; FOD=0   ## n=60 design 1805
#GD= 0 ; GDD=0; TOD=10000 ; FOD=0  ##n=1
def A(w) :
    return np.exp(-4*np.log(2)*((w-nuc)/dw)**2)
def P(w) :
    return GD* 2*np.pi*(w-nuc) + GDD/2*( 2*np.pi*(w-nuc))**2 + TOD/6*( 2*np.pi*(w-nuc))**3 + FOD/24*( 2*np.pi*(w-nuc))**4



# 2) je charge mon fichier de phase spectrale des miroirs TOD et je l'échantillonne :
wvl = np.loadtxt(r'C:\Users\ouille\Desktop\DEVIS_COMMANDES_MISSIONS\TOD mirrors\TOD1901simple.csv' , skiprows=1, delimiter=';', usecols=(0,))
R = np.loadtxt(r'C:\Users\ouille\Desktop\DEVIS_COMMANDES_MISSIONS\TOD mirrors\TOD1901simple.csv' , skiprows=1, delimiter=';', usecols=(1,))
GDDt = np.loadtxt(r'C:\Users\ouille\Desktop\DEVIS_COMMANDES_MISSIONS\TOD mirrors\TOD1901simple.csv' , skiprows=1, delimiter=';', usecols=(3,))


nulist = np.arange(c0/wvl[np.size(wvl)-1],c0/wvl[0], 0.0000001)
[nu,t] = ftAxis(nPoints,nulist.max())
GDD_nu = np.interp( nulist, c0/wvl[::-1] , GDDt[::-1])


# integrate GDD :
domega = 2*np.pi*(nulist[2]-nulist[1])
GD_TOD = np.zeros (np.size(nulist))
phase_TOD = np.zeros (np.size(nulist))
for i in np.arange (1, np.size(nulist)):
    GD_TOD [i] = GD_TOD[i-1] + GDD_nu[i]*domega
    phase_TOD[i] = phase_TOD [i-1] + GD_TOD [i] * domega
plt.figure()
plt.plot(nulist, GDD_nu)
plt.title ('GDD = f(w)')
plt.figure()
plt.plot(nulist, GD_TOD)
plt.title ('GD = f(w)  ---   GDD intégrée')
plt.figure()
plt.plot(nulist, phase_TOD)
plt.title ('$\Phi$ = f(w)  ---   GD intégré')

             
#on rajoute des zéros à gauche:
nulistz = np.append ( np.arange(0,nulist[0]-0.01, 0.01) , nulist )
phase_TODz = np.append ( np.zeros(np.size(nulistz) - np.size(GDD_nu)) , phase_TOD  )  

##interpolation et TF :
phase_interp = np.interp(nu, nulistz, phase_TODz)
pulse_t = ft( np.sqrt(A(nu)**2)*np.exp((1j)*(n*phase_interp+P(nu))) )
pulse_tbis = ft( np.sqrt(A(nu)**2) )


#plot figures :
fig = plt.figure(figsize=(13,4))
gs = gridspec.GridSpec(1, 2, wspace = 0.5)
ax1 = plt.subplot(gs[0,0])

ax1.plot(2*np.pi*nu, A(nu)**2,color='navy')

ax1.set_ylabel('I ($\omega$)', color='navy')
ax2=ax1.twinx()

ax2.plot(2*np.pi*nu, -n*phase_interp - P(nu) , color='red')

ax2.set_ylabel('$\Phi(\omega)$', color='red')
ax3=plt.subplot(gs[0,1])

ax3.semilogy(t,(pulse_t *np.conj(pulse_t)).real)
ax3.semilogy(t,(pulse_tbis *np.conj(pulse_tbis)).real)

ax3.set_title ('|E(t)|² with and without $\Phi(\omega)$')
ax3.set_xlim(-500,500)
ax3.set_ylim(ymin = 1e2, ymax=2e8)
ax3.set_xlabel('t(fs)')





plt.figure()
plt.plot(t,(pulse_t *np.conj(pulse_t)).real)
plt.plot(t,(pulse_tbis *np.conj(pulse_tbis)).real)
plt.xlim(-400,400)
plt.title ('|E(t)|² with and without $\Phi(\omega)$')



plt.figure()
plt.plot(2*np.pi*nu, -n*phase_interp - P(nu) , color='red')
plt.xlim(2.15,2.6)
plt.ylim(-850,-880)