# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 15:11:34 2019

@author: haessler, adapted from Matlab to Python by Marie
"""
#import libraries
import numpy as np
import scipy as sc
import sys
from matplotlib import pyplot as plt
if not 'C:\\Users\\ouille\\Desktop\\GitHub_users\\GitHub_MarieOuille\\python_functions' in sys.path : # put the path of the files : "smooth.py" and "peakdet.py" here 
    sys.path.insert(0,'C:\\Users\\ouille\\Desktop\\GitHub_users\\GitHub_MarieOuille\\python_functions')
from phys_constants import ct
from scipy import optimize


#define the function
def fit_calibration_stefan(nn,xn,deltaE,l,inc_angle,pxsize,gratinglines):
    

    plt.figure()
    s=16
    # --- Plot data that was originally in data set "xn vs. nn"
    plt.plot(nn,xn, 'o', label = 'identified peaks')
    plt.title('b) fit', fontsize=s)
    plt.xlabel('harmonic number $n$ (with an offset $n_0$)', fontsize=s)
    plt.ylabel('position (pixel number)', fontsize=s)

    
    const1 = ct.hbar /ct.e *2*np.pi*ct.c *gratinglines    # hbar/e * 2*pi*c /d
    def expr (n, B,Coff):
        return -l/pxsize / np.tan( np.arcsin( const1/((n+ Coff)*deltaE) - np.sin(inc_angle/180*np.pi) )  ) - B/pxsize
    [cf_,gof] = sc.optimize.curve_fit(expr, nn, xn, bounds = ([0,5], [0.2,12] ) )  #cf_ = [B, C]



    plt.plot( np.arange(nn[-1], nn[0], 0.01) , expr( np.arange(nn[-1], nn[0], 0.01), cf_[0] ,cf_[1] ) , color= 'red' , label = 'fit')
    plt.legend(fontsize=s)



    return [cf_, gof]