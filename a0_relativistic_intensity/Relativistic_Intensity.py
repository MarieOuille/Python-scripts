# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 11:37:26 2018

@author: ouille
"""

import numpy as np

c = 299792458 #m/s
ep0 = 8.85418782*10**-12  #m-3 kg-1 s4 A2
me = 9.10938356*10**-31 #kg
n =  1
e = 1.60217662*10**-19 #C

def I (l):
    return c**5 * n * ep0 * me**2 * 4 *np.pi**2 / ( 2* e**2* (l*10**-9)**2  ) *10**-4
                                                
l = int(input("What's the wavelength in nm? : "))
print ( 'The relativistic intensity is ' + str("%.1e"%(I(l))) + ' W/cm²')


####################

def a(Intensity, wl):
    return (Intensity/I(wl))**0.5
Int = float( input("What's the intensity in W/cm²? (e.g. : 3.4e+18) : "))
print( 'a0 is equal to ' + str(round(a(Int,l),2)) )