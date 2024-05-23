# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 17:13:37 2019

@author: ouille



TEST, copie de HHG_electrons_ions_CEPscan.py
pour faire une belle figure pour l'article CEP
avec les harmoniques en full

"""



#%%   #import libraries
import numpy as np
import re
from PIL import Image
import sys
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
if not r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' in sys.path : 
    sys.path.insert(0,r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' )
from phys_constants import ct
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['#BF360C','#F4511E', '#f3c1bb', 'white' ,'#d9e1f9' , '#a0b4f0','#6d8ce8','royalblue', '#2250d9' ,'#0000e6','darkblue', 'midnightblue','black' ])
import scipy as sc 
from pathlib import Path
#import scipy.signal



#%%  ##indicate file paths and load data.... things to change :
remoteDrive = r'Z:'        #file path    Z: or C: drive 
dataPath = 'SHHG\Jaismeen_Marie'      #dataPath = e.g : 'SHHG\Stefan_and_Marie\\2nd-campaign-May2019'  
outdir = r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python'



# data indications input :
# Scan parameter can be : 'CEP' or 'Gradient'
# Indications should be given in this format : dateCode = '20200630' ; file1 = 447 ; file_last = 489;  harm = 1 ; elec = 1  ; ions = 1  ; espectro =1 ; (Espec_ref = 197 ;) scan_param = 'Gradient' ;

                                                                                                                                                                                                                                                                                                           
#dateCode = '20200824'; file1 = 326 ; file_last = 351 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'CEP'
dateCode = '20200824'; file1 = 501 ; file_last = 526 ; harm = 1 ; elec = 0  ; ions = 0  ; espectro =0 ; scan_param = 'CEP'
#dateCode = '20200624' ; file1 = 280 ; file_last = 315 ; harm = 1 ; elec = 1  ; ions = 1  ; espectro =0 ; scan_param = 'Gradient' ;
#dateCode = '20200630' ; file1 = 447 ; file_last = 489;  harm = 1 ; elec = 1  ; ions = 1  ; espectro =1 ; Espec_ref = 197 ; scan_param = 'Gradient' ;

plasma_speed = 7 #nm/ps


save = 0        # save=0 : don't save, save = 1 : save figures






#%% Load data
fileNumbers = np.arange(file1, file_last + 1,1)

DataLog = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\\DataLog.txt' 
                                  
ShotNo = np.loadtxt(DataLog, skiprows = (2), usecols=(0,) )
if ShotNo[0] >0:   #offset in shot num.....
    offset_shot=int(ShotNo[0])
PulseEnergy = np.loadtxt(DataLog, skiprows = (2), usecols=(1,) )  #mJ
PulseDuration = np.loadtxt(DataLog, skiprows = (2), usecols=(2,) )  #fs
Chirp = np.loadtxt(DataLog, skiprows = (2), usecols=(3,) )  #fs²
WedgeInsert = np.loadtxt(DataLog, skiprows = (2), usecols=(4,) )  #mm
ParabolaPos = np.loadtxt(DataLog, skiprows = (2), usecols=(5,) )  #µm
PrepulseDelay = np.loadtxt(DataLog, skiprows = (2), usecols=(6,) ) #ps
SeqLen = np.loadtxt(DataLog, skiprows = (2), usecols=(7,) )   #ms
coms=np.loadtxt(DataLog, dtype= str, comments='CEP=',delimiter='\t',skiprows = (2), usecols=(11,) )


#reduce to the shot numbers we're interested in
ScanDelays = -PrepulseDelay[fileNumbers-offset_shot]

#extract the CEP value from the comments line :
CEPs=np.zeros(np.size(ShotNo))
for r, t in enumerate (np.arange(0,np.size(ShotNo)))  :
   # a = re.findall("[+-]\d+", coms[r])
    a=re.findall("[+-]?\d+\.\d+", coms[r])
    if np.shape(a)[0] == 0 :
        a_temp = (re.findall("[+-]?\d", coms[r]))
        if a_temp :
            a.append(a_temp[0])
            
    if np.shape(a)[0] == 0 :
        CEPs[r] = 0
    else :
        CEPs[r] = float(a[0])
CEPs = CEPs[fileNumbers-offset_shot] # List of CEP values per shot number

           
           
if scan_param == 'CEP' :        
    Params = CEPs     
    param_uniq = np.unique(CEPs)  # list of unique CEP values
    nb_param = np.unique(CEPs, return_counts = True)[1]  # number of occurences for each CEP value

elif scan_param == 'Gradient' :
    Params = ScanDelays
    param_uniq = np.unique(ScanDelays)   
    nb_param = np.unique(ScanDelays, return_counts = True)[1]           
                  
                  
                  
                  
    
    
    
    
    
    
                  
                  
                   
#%%  XUV HARMONICS SPECTRA                
if harm == 1 : 
    
    HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
              
                                               
    ## CALIBRATION ## 
    # harmonics ROI :
    if dateCode == '20200824':
        y_low = 100 ; y_high = 375
    elif dateCode == '20200630' :
        y_low = 100 ; y_high = 375
    elif dateCode == '20200624' :
        y_low = 107 ; y_high = 384 
    else :
        print('Be careful, you should check the HHG ROI')
        #%   #Check ROI ?
        plt.figure()
        plt.pcolor(np.loadtxt (HarmonicsFolder + 'tir' + str(file1) + '.dat' ) , cmap='inferno')
        plt.title('Select the ROI (2 points, first left then right)')
        ginput_HHG = plt.ginput(2)
        y_low = int(ginput_HHG[0][0]) ; y_high = int(ginput_HHG[1][0])
        plt.close()
        print('you chose : ' + str(y_low) + ' and ' + str(y_high) + ' for dateCode = ' + str(dateCode) + '\n' + 'you should add that in the script' )
    ##load calibration file for harmonics :
    Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_VIntegratedharmonics_tir31.txt', delimiter=',', dtype='float')




#%%

   #load calibration for harmonics energies :
    Ev = Ev[0:np.size(Ev)]
    rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
    rescale = np.append(rescale, rescale[-1]) #1element longer
    rescale=rescale[::-1]


    fig = plt.figure(figsize=(3,8))
    gs = gridspec.GridSpec(13, 1, left=0.06, top=0.88, wspace=0, hspace=0, figure=fig)

    #CWE cutoff/ max palsma frequency for a given wavelenght
    wp = np.sqrt (ct.e**2*ct.ne / (ct.me*ct.e0)) 
    lp = 2*np.pi*ct.c/wp
    Ev_max = ct.h*ct.c/lp/ct.e  
    
    for n, p in enumerate (param_uniq):
        for j, k in enumerate (fileNumbers):
            if CEPs [j]== param_uniq[n]:
                HarmonicsFile = np.loadtxt (HarmonicsFolder + 'tir' + str(k) + '.dat' ) 
                fig.add_subplot(gs[n, 0])
                plt.pcolor(Ev[::-1],np.arange(0,175), np.rot90(HarmonicsFile[:,y_low+100:y_high]), cmap = 'inferno')
                plt.text(35,80,'CEP='+str(p), color='white', fontsize=8)
                plt.xlim(xmin=11,xmax=40)
                plt.axvline(Ev_max, color='white')

            if n==12:
                plt.xlabel('Photon energy (eV)', fontsize=12)

                 
    plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\inf_20200824-tirs500qq_CEP_harm.png')


