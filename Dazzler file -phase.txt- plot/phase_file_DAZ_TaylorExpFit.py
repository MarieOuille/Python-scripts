# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 13:43:25 2018

@author: ouille
"""


import numpy as np
from matplotlib import pyplot as plt
import scipy as sc

### ordres 1,2,3,4 du dazzler numero 2
filepath = r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\Dazzler file -phase.txt- plot'
filename ='phase-testDazSN321percent'
file= str(filename) + '.txt'
output_dir=filepath

files = str(filepath) + "\\" + str(file)
lam = np.loadtxt(files, dtype='float', usecols = (0,)) # probe timing (ns)
phase = np.loadtxt(files, dtype='float', usecols = (1,))

#plt.figure()
#plt.plot(lam, phase)
#plt.title('Dazzler 2 phase from the phase.txt file')
#plt.ylabel('phase (rad)')
#plt.xlabel('wavelength (nm)')


### convert to omega
omega = 2 * np.pi * 299.792458 / lam
#plt.figure()
#plt.plot(omega, phase)
#plt.xlabel('frequency (rad/fs)')
#plt.ylabel('phase (rad)')
#



###fit to retrieve the different orders
###Taylor series expansion, order = 5
l_0 = 792   #800nm central wavelength
c0=299.792458   #speed of light (nm/fs)

def T (l,a,b,c,d,e,f):
    return a + b*2*np.pi*c0*(1/l-1/l_0) + c/2*((2*np.pi*c0)**2)*(1/l-1/l_0)**2 + d/6*((2*np.pi*c0)**3)*(1/l-1/l_0)**3 +  e/24*((2*np.pi*c0)**4)*(1/l-1/l_0)**4 + f/120*((2*np.pi*c0)**5)*(1/l-1/l_0)**5 

initial_guesses=[-20,  0,   3000, 3000, 0, 0]
func_fit = sc.optimize.curve_fit(T, lam,phase, initial_guesses)
#func_fit = sc.optimize.curve_fit(T, lam,phase)

a=func_fit[0][0]
b=func_fit[0][1]
c=func_fit[0][2]
d=func_fit[0][3]
e=func_fit[0][4]
f=func_fit[0][5]
plt.figure()
plt.plot(lam, T(lam,a,b,c,d,e,f), label='Taylor series exp fit')
plt.plot(lam, phase, color='r', label='Dazzler phase.txt file')
plt.ylabel('phase (rad)')
plt.xlabel('wavelength (nm)')
plt.legend(loc='best', fontsize='large')
plt.title ('GDD = ' + str(round(func_fit[0][2])) +' fs2 \n TOD = ' + str(round(func_fit[0][3])) + ' fs3  \n  FOD = ' + str(round(func_fit[0][4]))  + ' fs4   \n  5th order = ' + str(round(func_fit[0][5])) + ' fs5 \n')
plt.savefig( r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\Dazzler file -phase.txt- plot\\Taylor series expansion fit of the phase file obtained from Wizzler.png',dpi=600,  bbox_inches='tight')