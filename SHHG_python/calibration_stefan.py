# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:50:12 2019

@author : Haessler, adapted from Matlab to Python 


description of the script :
This script uses an harmonics spectrum to calibrate the correspondance between position <-> photons energies
First it finds the coordinates of the harmonics
then it calls a function that fits those coordinates using the expression position = expr (harmonics number ,depending on parameters B and C  + constants)
Once B(position of the detector) and C(harmonics number offset) are found, we can convert the x axis into photon energies (eV)



"""



#%%
#libraries
import sys
import numpy as np
from matplotlib import pyplot as plt
if not 'C:\\Users\\ouille\\Desktop\\GitHub_users\\GitHub_MarieOuille\\python_functions' in sys.path : # put the path of the files : "smooth.py" and "peakdet.py" here 
    sys.path.insert(0,'C:\\Users\\ouille\\Desktop\\GitHub_users\\GitHub_MarieOuille\\python_functions')
from smooth import smooth
from peakdet import peakdet
if not 'C:\\Users\\ouille\\Desktop\\GitHub_users\\GitHub_MarieOuille\\Python-scripts\\SHHG_python' in sys.path :  # put the path of the file : "fit_calibration_stefan.py" here
    sys.path.insert(0,'C:\\Users\\ouille\\Desktop\\GitHub_users\\GitHub_MarieOuille\\Python-scripts\\SHHG_python' )
from fit_calibration_stefan import fit_calibration_stefan
from phys_constants import ct
import scipy as sc



#%%

# FILENAMES
#filepath = r'C:\Users\ouille\Desktop\SHHG analysis\HHG calibration'
#filename  ='1402-24fs-refspec.txt';  #vertically integrated harmonics spectrum
#filepath = r'C:\Users\ouille\Desktop\Python_travail_en_cours'
#filename  ='tir310_VerticallyIntegrated.txt';  #vertically integrated harmonics spectrum
#filename  ='tir66_VerticallyIntegrated.txt';  #vertically integrated harmonics spectrum
#filename  ='tir171_VerticallyIntegrated.txt';  #vertically integrated harmonics spectrum

##After lockdown, Jaismeen and Marie :
filepath = r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files'
#filename = 'tir601_VerticallyIntegrated.txt'
filename = 'tir38_VerticallyIntegrated.txt'



# parameters to tweak
lambd = 780 # laser central wavelength in nanometer
deltaE =ct.hbar*2*np.pi*ct.c /(lambd*1e-9) /ct.e  # photon energy spacing, if odd *and even* harmonics present
smoothwidth = 7 # over how many pixels to make a running-average-smoothing of the spectra # this makes peak-finding easier in noisy data.
peaksearch_contrast = 20000 #how many counts difference you require between a 'minimum' and a 'maximum' in the peak search... used to be 1
pxwidth = 1 #number of pixels left and right of the found maxima to include in the center-of-weight search
gratinglines=600e3 #600 lines per mm
cam_dist = l = 482e-3 #m       #cam_dist = l = 469e-3 #m   
inc_angle = 85.3 # degres
pxsize = 1.65e-4
const1 = ct.hbar/ct.e *2*np.pi*ct.c *gratinglines  


#save the txt calibration file or not ? 1 for yes, 0 for no
save = 1  
outdir = r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files'




#%% Load data :
    
data0=np.loadtxt(filepath + '\\' + filename)   #spectre d'harmoniques intégré verticalement
data=np.zeros((640,2))
data[:,1]=data0[::-1] #spectre
data[:,0]=np.arange(1,np.size(data0)+1)  #liste de 1 à 640 par pas de 1


    
    
#%% Find harmonics peaks :
    
spc_max=smooth(data[:,1],smoothwidth)
[maxtab_o, mintab_o]=peakdet(spc_max,peaksearch_contrast, data[:,0])
maxtab=maxtab_o[0:-1,:]
mintab=mintab_o


s=16
plt.figure()
plt.plot(data[:,0], spc_max, color='Navy')   
plt.xlabel('position (pixel number)', fontsize=s)
plt.ylabel('signal amplitude (a.u.)', fontsize=s)
plt.title('a) vertically integrated raw data', fontsize=s)
plt.xlim(0,640)

cowtab = np.zeros((len(maxtab),2))
for i in (np.arange(0,len(maxtab))) :
    if maxtab[i,0] < pxwidth:
        ppxwidth=maxtab[i,0]-1
    elif maxtab[i,0]+pxwidth > len(spc_max) :
        ppxwidth=len(spc_max)-maxtab[i,0]
    else :
        ppxwidth = pxwidth
    cowtab[i,:] = ([ sum( spc_max[int((maxtab[i,0]-ppxwidth-1)):int((maxtab[i,0]+ppxwidth)) ] * data[ int( (maxtab[i,0]-ppxwidth-1)): int(( maxtab[i,0]+ppxwidth )) ,0 ] ) / sum( spc_max[int ((maxtab[i,0]-ppxwidth-1)) : int((maxtab[i,0]+ppxwidth)) ] )    ,  int( maxtab[i,1])  ])
    if i == 0:
        plt.plot(cowtab[i,0],cowtab[i,1], 'o', label='identified peaks')
    else :
        plt.plot(cowtab[i,0],cowtab[i,1], 'o')
        
plt.legend(fontsize=s)
#    area(data(maxtab(i)-ppxwidth:maxtab(i)+ppxwidth,1), spc_max(maxtab(i)-ppxwidth:maxtab(i)+ppxwidth))
nn=np.arange(0,len(maxtab))[::-1]   #harmonics numbers with an offset
xn=cowtab[:,0]   #x coordinates of the peaks 


         
#%% FIT and rescale :   
         
# Make a fit using (xn, nn) to find B and C :
[cf_,gof] = fit_calibration_stefan(nn,xn,deltaE,l,inc_angle,pxsize,gratinglines)
B = cf_ [0]
Coff = cf_[1]

print('************** The fit gives: ***************')
print('\n','For a constant imposed harmonic spacing of ',str(deltaE),' eV,')
print('\n','the lowest multiple of that spacing" = ',str(Coff),',')
print('\n','and the detector "blue edge" at x_0 = ',str(B*1000),' mm.', '\n')


#convert pixels to energies (eV) :  
Ev = const1 / (    np.sin( np.arctan(-l/(data[:,0]*pxsize+B)) )   +    np.sin(inc_angle/180*np.pi)     )    #photon energy

##save the text file used later for calibration of all the data :
if save == 1:
    np.savetxt((outdir + '\\NEW_Ev-calib_'+filename[0:-4]+'.txt'), Ev, newline=',') 
    print ('a new calib file has been created and saved in' + str(outdir))

#converts the intensity
rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
rescale = np.append(rescale, rescale[-1])  #just to make the list 1 element longer

                   
         
                   
                   
                   

#%% PLOT FIGURE           
                   
plt.figure()
plt.plot(Ev, rescale*spc_max, color='navy')  #plots the spectrum with x AND y axis recalibrated
plt.xlabel('Energy [eV]', fontsize=s)
plt.ylabel('Signal level (a.u.)', fontsize=s)
#plt.title(filename[0:-4])
plt.title('c) calibrated and rescaled XUV spectrum', fontsize=s)

harmpos = (np.arange(Coff,Coff+15,1))*deltaE
for i in np.arange (0,len(harmpos)) :
    if i == 0 :
        plt.plot([harmpos[i],harmpos[i]], [0,1e9],color='0.8', label='lines spaced by $\Delta E$')    #vertical lines at harmonics positions
    else :
        plt.plot([harmpos[i],harmpos[i]], [0,1e9],color='0.8')      
plt.legend(fontsize=s-1)
plt.xlim(Ev.min(), 35)     #plt.xlim([min(Ev), Ev[round(cowtab[0,0])]+5])
plt.ylim([0, max(rescale*spc_max)])