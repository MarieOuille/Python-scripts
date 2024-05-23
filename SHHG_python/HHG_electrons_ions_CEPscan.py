# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 17:13:37 2019

@author: ouille
This is a (almost) a copy of the matlab scripts for SHHG analysis



THINGS TO DO : 
    
- ajouter plasma speed/gradient length
- peut-être créer un fichier pour tous les paramètres/calib
- créer dossier date pour sauvegarde des données s'il n'existe pas deja dans outdir
- check that there aren't too many hot pixels in LANEX, otherwise --> display error or increase from 100counts difference to 150 or 200



optional stuff :
    
- TRY : CEP scan with ions
- TRY : CEP scan with espectro
- TRY : z-scan
- TRY : random gradient scan

- add text to explain how to find alpha, a1 and d for lanex
- add text to explain how to calibrate the harmonics spectrum
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
#custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['#BF360C','#F4511E', '#f3c1bb', 'white' ,'#d9e1f9' , '#a0b4f0','#6d8ce8','royalblue', '#2250d9' ,'#0000e6','darkblue', 'midnightblue','black' ])
import scipy as sc 
from pathlib import Path
#import scipy.signal


#%% custom cmap
#Create your own color map :
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'blue', 'green', 'yellow', 'orange', 'red', 'darkred'])





#%%  ##indicate file paths and load data.... things to change :
remoteDrive = r'Z:'        #file path    Z: or C: drive 
dataPath = 'SHHG\Jaismeen_Marie'      #dataPath = e.g : 'SHHG\Stefan_and_Marie\\2nd-campaign-May2019' 
#dataPath='SHHG\Stefan_and_Marie\\2nd-campaign-May2019'  
#dataPath='SHHG\LILU'
#outdir = r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python'
outdir = r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python'



# data indications input :
# Scan parameter can be : 'CEP' or 'Gradient'
# Indications should be given in this format : dateCode = '20200630' ; file1 = 447 ; file_last = 489;  harm = 1 ; elec = 1  ; ions = 1  ; espectro =1 ; (Espec_ref = 197 ;) scan_param = 'Gradient' ;

                                                                                                                                                                                                                                                                                                           
#dateCode = '20200824'; file1 = 326 ; file_last = 351 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'CEP'
#dateCode = '20200824'; file1 = 501 ; file_last = 526 ; harm = 0 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 486 ; scan_param = 'CEP'

#dateCode = '20200824'; file1 = 149 ; file_last = 169 ; harm = 1 ; elec = 1  ; ions = 1  ; espectro =1 ; Espec_ref = 147 ; scan_param = 'Gradient'



#dateCode = '20200624' ; file1 = 280 ; file_last = 315 ; harm = 1 ; elec = 1  ; ions = 1  ; espectro =0 ; scan_param = 'Gradient' ;
#dateCode = '20200630' ; file1 = 447 ; file_last = 489;  harm = 1 ; elec = 1  ; ions = 1  ; espectro =1 ; Espec_ref = 197 ; scan_param = 'Gradient' ;
#dateCode = '20200824'; file1 = 580 ; file_last = 607 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 486 ; scan_param = 'CEP'
#dateCode='20200923'; file1 = 288 ; file_last = 313 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'CEP'


## ESPECTRO CEP 2eme scan
#dateCode = '20200904'; file1 = 326 ; file_last = 351 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 316 ; scan_param = 'CEP'


#dateCode = '20200902'; file1 = 253 ; file_last = 280 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 219 ; scan_param = 'CEP'

#CEP locked gradient scans
#dateCode = '20200902'; file1 = 298 ; file_last = 310 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 219 ; scan_param = 'Gradient'
#dateCode = '20200902'; file1 = 312 ; file_last = 324 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 219 ; scan_param = 'Gradient'


##1ER SCAN ESPEXCTRO CEP THESE
#dateCode = '20200902'; file1 = 221 ; file_last = 246 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 219 ; scan_param = 'CEP'




##scan graident
#dateCode = '20200910'; file1 = 138 ; file_last = 148 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 219 ; scan_param = 'Gradient'






#dateCode = '20200702'; file1 = 485 ; file_last = 549 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 458 ; scan_param = 'CEP'

#dateCode = '20200701'; file1 = 134 ; file_last = 216 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 215 ; scan_param = 'Gradient'

#LWFA 6fs spectro
#dateCode = '20200701'; file1 = 219 ; file_last = 283 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 215 ; scan_param = 'Gradient'


#dateCode = '20200824'; file1 = 528 ; file_last = 539 ; harm = 1 ; elec = 1  ; ions = 1  ; espectro =1 ; Espec_ref = 486 ; scan_param = 'Parabola_pos'
#dateCode = '20200824'; file1 = 488 ; file_last = 499 ; harm = 1 ; elec = 0  ; ions = 0  ; espectro =0 ; Espec_ref = 486 ; scan_param = 'Parabola_pos'
#dateCode = '20200824'; file1 = 488 ; file_last = 499 ; harm = 1 ; elec = 0  ; ions = 0  ; espectro =1 ; Espec_ref = 486 ; scan_param = 'Parabola_pos'
#dateCode = '20200824'; file1 = 444 ; file_last = 455 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 486 ; scan_param = 'Parabola_pos'

#dateCode = '20200824'; file1 = 567 ; file_last = 578 ; harm = 1 ; elec = 0  ; ions = 0  ; espectro =0 ; Espec_ref = 486 ; scan_param = 'Parabola_pos'
#dateCode = '20200902'; file1 = 325 ; file_last = 350 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 219 ; scan_param = 'CEP'
#dateCode = '20200824'; file1 = 315 ; file_last = 323 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 486 ; scan_param = 'Parabola_pos'



#dateCode = '20200824'; file1 = 353 ; file_last = 362 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Parabola_pos'
dateCode = '20200824'; file1 = 391 ; file_last = 406 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Parabola_pos'


#
#dateCode = '20200619' ; file1 = 26 ; file_last = 46 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ;  scan_param = 'Gradient' ; #ions=1 also possible
#plasma_speed = 19 #2020.06.26 SDI meas #1001 - 1098
#

#•pinhole
#dateCode = '20200619' ; file1 = 457 ; file_last = 500 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;
#dateCode = '20200619' ; file1 = 457 ; file_last = 550 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;
#plasma_speed=27 #8.8/np.cos(55*np.pi/180)**2
#dateCode = '20200619' ; file1 = 120 ; file_last = 172 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;
#plasma_speed =21 # 6.8/np.cos(55*np.pi/180)**2
#20200619, no pinhole
#dateCode = '20200619' ; file1 = 382 ; file_last = 396 ; harm = 1 ; elec = 1  ; ions = 1  ; espectro =0 ; scan_param = 'Gradient' ;
#last = 446 if long gradients


#fileNumbers = np.append( np.arange(123,149), np.arange(162,173) ) 

#thèse, scans haute intensité
#dateCode = '20200626' ; file1 = 25 ; file_last = 69 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;   #30 fs
#plasma_speed = 19 #2020.06.26 SDI meas #1001 - 1098
#dateCode = '20200619' ; file1 = 743 ; file_last = 764 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;   #7 fs
#dateCode = '20200624' ; file1 = 280 ; file_last = 315 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;   #4 fs
#dateCode = '20200904' ; file1 = 240 ; file_last = 261 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;   #4 fs
#dateCode = '20200904' ; file1 = 70 ; file_last = 94 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;   #6.1 fs
#plasma_speed=18


#thèse LWFA 3fs
#dateCode = '20200904' ; file1 = 356 ; file_last = 406 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =1 ; Espec_ref = 316; scan_param = 'Gradient' ;   #6.1 fs
#dateCode = '20200904' ; file1 = 411 ; file_last = 461 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; Espec_ref = 316; scan_param = 'Gradient' ;   #6.1 fs




#24fs et 24+500fs2 chirp, comme dans le papier de Stefan... a revoir, calib HHG (Ev file à refaire) + Lanex angles etc
#dateCode = '20190605' ; file1 = 28 ; file_last = 48 ; harm = 1 ; elec = 1  ; ions = 1  ; espectro =0 ; scan_param = 'Gradient' ;   #6.1 fs
##dateCode = '20190605' ; file1 = 186 ; file_last = 206 ; harm = 1 ; elec = 1  ; ions = 1  ; espectro =0 ; scan_param = 'Gradient' ;   #6.1 fs
#plasma_speed = 24 #SDI on that day
#



#dateCode = '20191202' ; file1 = 229 ; file_last =265 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ; scan_param = 'Gradient' ;   #6.1 fs


#dateCode = '20200629'; file1 = 262 ; file_last = 263 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ;scan_param = 'CEP'   #+1ps

#dateCode = '20200626'; file1 = 466 ; file_last = 550 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ;scan_param = 'CEP'   #+1ps
#dateCode = '20200626'; file1 = 557 ; file_last = 643 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ;scan_param = 'CEP'   #-1ps
#dateCode = '20200626'; file1 = 649 ; file_last = 738 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ;scan_param = 'CEP'   #-2ps
#dateCode = '20200626'; file1 = 746 ; file_last = 827 ; harm = 1 ; elec = 1  ; ions = 0  ; espectro =0 ;scan_param = 'CEP'   #-3ps




#PUT THISWITH SPATIAL CALIB VS DATE
#Espec = np.loadtxt(remoteDrive + '\\' + dataPath + '\dsdE_20200630.txt');  flip =0
#Espec = np.loadtxt(remoteDrive + '\\' + dataPath + '\dsdE_20200902.txt')   ;  flip = 1
#Espec = np.loadtxt(remoteDrive + '\\' + dataPath + '\dsdE_20200701_until_20200706.txt');  flip =0
#Espec = np.loadtxt(remoteDrive + '\\' + dataPath + '\dsdE_Aug2020.txt')   ;  flip =0



#plasma_speed =  24#17#8  # 8.8 #7.1  #9.2#7.1#8.2# 8#8.8 #6.8  # 7 #5.1 for 30fs, 7 for 4fs, and let's say 6 for 7fs? #nm/ps
#17.1 on 20190521
#13.3 on 20190605


save = 1        # save=0 : don't save, save = 1 : save figures



#
#flip=1

#%% Load data
fileNumbers = np.arange(file1, file_last + 1,1)
#fileNumbers = np.array([580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,607])

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
    #a = re.findall("[+-]\d+", coms[r]) #[56:-1])
    a=re.findall("[+-]?\d+\.\d+", coms[r])
    if np.shape(a)[0] == 0 :
        a_temp = (re.findall("[+-]?\d", coms[r])) #[56:-1]))
        if a_temp :
            a.append(a_temp[0])
            
    if np.shape(a)[0] == 0 :
        CEPs[r] = 0
    else :
        CEPs[r] = float(a[0])
CEPs = CEPs[fileNumbers-offset_shot] # List of CEP values per shot number


###### essai décaler 0 de la CEP pour 20200626
#for i, j in enumerate (CEPs):
#    if j <= 0.5:   ###-2ps
#        CEPs[i]=CEPs[i]+2.5
#    elif j > 0.5:
#        CEPs [i] = CEPs[i]-4
#    if j >= -2:  #######--1ps
#        CEPs[i]=CEPs[i]-1
#    elif j < -2:
#        CEPs [i] = CEPs[i]+5.5          
           
           
if scan_param == 'CEP' :        
    Params = CEPs     
    param_uniq = np.unique(CEPs)  # list of unique CEP values
    nb_param = np.unique(CEPs, return_counts = True)[1]  # number of occurences for each CEP value

elif scan_param == 'Gradient' :
    Params = np.round(ScanDelays,1)
    param_uniq = np.unique(np.round(ScanDelays,1))   
    nb_param = np.unique(np.round(ScanDelays,1), return_counts = True)[1]           
                  
                  
elif scan_param == 'Parabola_pos' :
    Params = ParabolaPos[fileNumbers-offset_shot]
    param_uniq = np.unique(Params)   
    nb_param = np.unique(Params, return_counts = True)[1] 
    
    
    
    
                  
                  
                   
#%%  XUV HARMONICS SPECTRA                
if harm == 1 : 
    
    HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
              
                                               
    ## CALIBRATION ## 
    # harmonics ROI :
    if dateCode == '20200824':
        y_low = 100 ; y_high = 375
    elif dateCode == '20200630' :
        y_low = 100 ; y_high = 375
    elif dateCode == '20200626' :
        y_low = 110 ; y_high = 385      ##checked
    elif dateCode == '20200629' :
        y_low = 110 ; y_high = 385     
    elif dateCode == '20200624' :
        y_low = 107 ; y_high = 384 
    elif dateCode == '20200923' :
        y_low =104 ; y_high = 378
    elif dateCode >= '20200902' and dateCode <='20200910' :
        y_low =99 ; y_high = 376
    elif dateCode == '20200702' :
        y_low = 110 ; y_high = 382
    elif dateCode =='20200701':
        y_low = 106; y_high = 383
    elif dateCode == '20200619' :
        y_low=106 ; y_high = 385
    elif dateCode == '20191202' :
        y_low=78 ; y_high = 403
    elif dateCode == '20190605' :
        y_low=39 ; y_high=359
        Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_tir38_VerticallyIntegrated.txt', delimiter=',', dtype='float')

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






   #load calibration for harmonics energies :
    Ev = Ev[0:np.size(Ev)]
    rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
    rescale = np.append(rescale, rescale[-1]) #1element longer
    rescale=rescale[::-1]
    
    HarmonicsScan = []
    for j, k in enumerate (fileNumbers):
        HarmonicsFile = np.loadtxt (HarmonicsFolder + 'tir' + str(k) + '.dat' ) 
        HarmonicsScan.append(list (np.sum(HarmonicsFile[:,y_low:y_high],1)) * rescale)  #integrated and put in an array
    
    #flip 90degree to get the delays on the horizontal axis, more convenient
    HarmScan_r = np.zeros ((np.shape(HarmonicsScan)[1],np.shape(HarmonicsScan)[0]))
    for k, osef in enumerate (np.zeros(np.shape(HarmonicsScan)[0])):
        for j, osef_aussi in enumerate (np.zeros(np.shape(HarmonicsScan)[1])):
            HarmScan_r[j,k] = HarmonicsScan[k][j]
            
            
            
            
    # organize data to average if the same parameter value occurs more than once :  
              
    Bharm = np.zeros((640,np.size(param_uniq)))
    for j, valuej in enumerate (param_uniq):
        Aharm = np.zeros((640,nb_param[j]))
        count = 0
        for i, value in enumerate (fileNumbers):
            if Params [i] == valuej :
                Aharm[:,count] = HarmScan_r[:,i]
                if count < nb_param[j] :
                    count = count +1 
        Bharm[:,j]=np.sum(Aharm,1)/nb_param[j]    #average over all the measurements for the same CEP value
      #  Bharm[:,j]=Bharm[:,j]-Bharm[:,-1]      #substract the last HHG spectrumto remove BG
                 
    #CWE cutoff/ max palsma frequency for a given wavelenght
    wp = np.sqrt (ct.e**2*ct.ne / (ct.me*ct.e0)) 
    lp = 2*np.pi*ct.c/wp
    Ev_max = ct.h*ct.c/lp/ct.e                



#
#    for i, j in enumerate(Bharm[0,:]):
#        Bharm[:,i] = Bharm[:,i]/np.max(Bharm[:,i])
#
#










#%% electrons 
ElecScan = []
if elec == 1:
    # data folder 
    SpaElectronsFolder = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\electrons\\'                     
                                 

    # LANEX calibration : 
    if int(dateCode) >= 20200611 and int(dateCode) <= 20200624 :             #validated
        px_ll_lanex = [193,522] #lower left corner of the LANEX : y,x
        px_ur_lanex = [777,1063] #upper right corner of the LANEX : y,x
        a1 = 0.106 ; d = 0.208 #m
        alpha = np.deg2rad(31.6) #rad
    elif int(dateCode) == 20200625 :                                         #validated
        px_ll_lanex = [193,416] #lower left corner of the LANEX : y,x
        px_ur_lanex = [777,945] #upper right corner of the LANEX : y,x
        a1 = 0.113 ; d = 0.236 #m
        alpha = np.deg2rad(38) #rad
    elif int(dateCode) == 20190605 :                                         #validated
        px_ll_lanex = [273,479] #lower left corner of the LANEX : y,x
        px_ur_lanex = [791,965] #upper right corner of the LANEX : y,x
        a1 = 0.096 ; d = 0.236 #m checked
        alpha = np.deg2rad(36) #rad   checked
    elif int(dateCode) >= 20200626 and int(dateCode)<=20200630 :             #validated
        px_ll_lanex = [193,417] #lower left corner of the LANEX : y,x
        px_ur_lanex = [777,965] #upper right corner of the LANEX : y,x
        a1 = 0.117 ; d = 0.235 #m
        alpha = np.deg2rad(38) #rad
    elif int(dateCode) >= 20200701 and int(dateCode)<=20200707 :             #validated
        px_ll_lanex = [188,330] #lower left corner of the LANEX : y,x
        px_ur_lanex = [777,885] #upper right corner of the LANEX : y,x
        a1 = 0.151 ; d = 0.291 #m
        alpha = np.deg2rad(41.5) #rad
    elif int(dateCode) == 20200824 :                                         #validated
        px_ll_lanex = [185,330] #lower left corner of the LANEX : y,x
        px_ur_lanex = [775,885] #upper right corner of the LANEX : y,x
        a1 = 0.161 ; d = 0.301 #m
        alpha = np.deg2rad(41.1) #rad
    elif int(dateCode) == 20191202 :                                         #checked
        px_ll_lanex = [208,446] #lower left corner of the LANEX : y,x
        px_ur_lanex = [799,1007] #upper right corner of the LANEX : y,x
        a1 = 0.051 ; d = 0.085 #m
        alpha = np.deg2rad(18) #rad
#    elif int(dateCode) == 20200904:                                         #validated
#        px_ll_lanex = [460,343] #lower left corner of the LANEX : y,x
#        px_ur_lanex = [1015,953] #upper right corner of the LANEX : y,x
#        a1 =  0.1384 ; d = 0.254  #m
#        alpha = np.deg2rad(36.54) #rad
    elif int(dateCode) >= 20200902 and int(dateCode) <= 20200922  :                                         #validated
        px_ll_lanex = [350,450] #lower left corner of the LANEX : y,x
        px_ur_lanex = [925,1000] #upper right corner of the LANEX : y,x
        a1 = 0.127 ; d = 0.254 #m
        alpha = np.deg2rad(36.54) #rad
    elif int(dateCode) == 20200923:                                         #validated
        px_ll_lanex = [340,260] #lower left corner of the LANEX : y,x
        px_ur_lanex = [950,850] #upper right corner of the LANEX : y,x
        a1 = 0.0275 ; d =  0.112#m
        alpha = np.deg2rad(36.6) #rad
    else :
        print('\n Be careful, there is no LANEX calibration')
        TEST_lanex = np.array( Image.open (SpaElectronsFolder + 'tir' + str(file1) + '.tiff' )  )                          
        plt.figurefigsize=(10,10)
        plt.pcolor(TEST_lanex) #,  vmin =1, vmax=100)
        plt.title('Select 2 points : 1) lower left & 2) upper right corner')
        text_ROI = plt.ginput(2)
        plt.close()
        print('For dateCode = ' + str(dateCode) + ', ')
        print('px_ll_lanex =' + str( int(text_ROI[0][1])) + ', ' + str( int(text_ROI[0][0]) )   )
        print('px_ur_lanex =' + str( int(text_ROI[1][1])) + ', ' + str( int(text_ROI[1][0]) )   )
        print('The execution will be aborted, you still need alpha, a1 and d \n')
        sys.exit()

    LANEX_width = 0.065 #m
    
    
    
    
    
    
                                 ## --------------------------------- ## 
                 
    yll =  px_ll_lanex[0] ; xll =  px_ll_lanex[1]     # ROI  
    yur =  px_ur_lanex[0] ; xur =  px_ur_lanex[1]     # ROI  
    
    s = round(LANEX_width/(xur-xll),7)  #m/pixel
    
    # coordinates list, downsized by dsz to make it faster
    #dsz = 2      
    if espectro ==1 :
        dsz=1
    elif espectro==0:
        dsz =2
    
    #coordinates changes :
    px_temp = np.linspace(0,int(xur-xll)-1,xur-xll)[0:-1:dsz]
    py_temp = np.linspace(0,yur-yll-1,yur-yll)[0:-1:dsz]
    #create coordinate mesh :
    [px,py]=np.meshgrid(px_temp,py_temp)
    #with real dimensions in the LANEX plane, origin = intersection with target normal
    A = a1+px*s
    B = -(yur-yll)*s/2 + py*s
    
         

    
    
                                 ## --------------------------------- ##       
                                                
    if espectro == 1 :
        
        ## FIRST : check the zero deflection point with the reference file (ie only with the pinhole, no magnets)
        #x, y, and values in px
        if flip == 1 :
            LANEXpx_ref = np.flipud (np.array( Image.open (SpaElectronsFolder + 'tir' + str(Espec_ref) + '.tiff' )  )[yll:yur, xll:xur] )
        else :
            LANEXpx_ref = np.array( Image.open (SpaElectronsFolder + 'tir' + str(Espec_ref) + '.tiff' )  )[yll:yur, xll:xur]
            
 
        #we need to get rid of hot pixels
        temp = sorted(LANEXpx_ref.flatten())[::-1]    #array -> list, then organize from highest to lowest value
        n_idx=[]
        for i, value in enumerate (temp[0:20]):  #search only within the 10 highest values for >= 100 counts differences
            if temp[i]-temp[i+1] >= 10:
                n_idx=np.append(n_idx,i)
        if np.size(n_idx) >0 :
            idx = np.where(LANEXpx_ref >= temp[int(n_idx.max())])
            for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
                LANEXpx_ref[idx[0][j], idx[1][j]] = LANEXpx_ref[idx[0][j]-1, idx[1][j]]  
            
        LANEXpx_ref = LANEXpx_ref[0:-1:dsz, 0:-1:dsz]           # downsize ie take one every dsz points
    
    
        #with real dimensions in the LANEX plane, origin = intersection with target normal
        LANEXm_ref = LANEXpx_ref / (s**2)  #now in Coulomb/m²    
        
        # next : find the center of the spot                      
        plt.figure(figsize=(10,10))
        plt.title('Select the zero deflexion point', fontsize=16)
        plt.pcolor(A,B,LANEXm_ref, cmap = 'viridis')
        plt.xlabel('a (m)')
        plt.ylabel('b (m)') 
        plt.xlim( A.min() + (A.max()-A.min())/2-0.25*(A.max()-A.min())    , A.min() + (A.max()-A.min())/2+0.25*(A.max()-A.min())    )
        plt.ylim( B.min() + (B.max()-B.min())/2-0.25*(B.max()-B.min())    , B.min() + (B.max()-B.min())/2+0.25*(B.max()-B.min())  )
        [(aref, bref)] = plt.ginput(1, show_clicks = True)
        plt.close()

    
        ## NOW let's treat the data
        # first load the spectrometer calibration :

            
        EE = Espec[:,0]
        dSdE= Espec[:,1]
        S = Espec[:,2]
        teta = Espec[:,3]
        dEE = Espec[:,4]
        

        # START LOOP HERE
        for j, k in enumerate (fileNumbers):
  
        # load file
            SpaElectronsFile = np.array (   Image.open (SpaElectronsFolder + 'tir' + str(k) + '.tiff' )    )
            
            #x, y, and values in 
            if flip == 1 :
                LANEXpx = np.flipud(SpaElectronsFile[yll:yur,xll:xur]  ) #only take the ROI, in Coulomb/px   
            else :
                LANEXpx = SpaElectronsFile[yll:yur,xll:xur]  
        
        
        

            
        
            #test : downsize ie take one every dsz points
            LANEXpx = LANEXpx[0:-1:dsz, 0:-1:dsz]
    
        
            #with real dimensions in the LANEX plane, origin = intersection with target normal
            LANEXm = LANEXpx / (s**2)  #now in Coulomb/m²  
#            LANEXm = LANEXm - 0.95*LANEXm_ref #test
#            for ff, fff in enumerate (LANEXm[0,:]):
#                for vv, vvv in enumerate( LANEXm[:,0]):
#                    if LANEXm[vv,ff]<0:
#                        LANEXm[vv,ff]=0
##            
##            
#            
            
            
                        #we need to get rid of hot pixels
            temp = sorted(LANEXm.flatten())[::-1]    #array -> list, then organize from highest to lowest value
            n_idx=[]
            for i, value in enumerate (temp[0:100]):  #search only within the 10 highest values for >= 100 counts differences
                if temp[i]-temp[i+1] >= 100 and temp[i]>20003830107:
                    n_idx=np.append(n_idx,i)
            idx = np.where(LANEXm >= temp[int(n_idx.max())])
            for i,w in enumerate (np.arange(0,int(np.size(idx)/2))):
                LANEXm[idx[0][w], idx[1][w]] = LANEXm[idx[0][w]-1, idx[1][w]]  
            
            
            
            
            
            
            
            
            
            
            
            
    
    
            if j ==0:
                # define the integration area                  
                plt.figure(figsize=(10,10))
                plt.title('Select the signal area by clicking 2 points : \n left and right', fontsize=16)
                plt.pcolor(LANEXm, cmap = 'viridis', vmin=0, vmax=50)
                zone = plt.ginput(2, show_clicks = True)
                plt.close()
                
                alim1 = int(zone[0][0])
                alim2 = int(zone[1][0])
    
          

            
            Positions = np.round( (B[:,0]-bref)*1000  ,1)              
            Spectrum = np.sum ( LANEXm[:,alim1:alim2], axis=1 ) 
#            Spectrum_ref = np.sum ( LANEXm_ref[:,alim1:alim2], axis=1 ) 
#            Spectrum=Spectrum-Spectrum_ref
          #  Spectrum = Spectrum - np.mean(Spectrum[0:30]) # test remove background
            
            idx0 = np.where(Positions >= np.round(S.min(),1))[0][0]    +1    
    
            E_interp = np.interp(Positions[idx0::], S[::-1], EE[::-1] ) 
                
                                
            if j==0:
                ElecScan_r=np.zeros((np.size(E_interp), np.size(fileNumbers)))
                                
            dBdE_interp = np.interp(Positions[idx0::], S[::-1], dSdE[::-1])
            

            
            #ElecScan_r[:,j] = -Spectrum[idx0::]*dBdE_interp
            ElecScan_r[:,j] = -Spectrum[idx0::]*dBdE_interp
#           ###20200904


#################### pour avoir les raw image des spectros electrons
##            i1=8
##            i2=2     
##### 20200902
#            i1=1
#            i2=8
#                            
#            if j==i1:
#                plt.figure(figsize=(5,4))
#                plt.pcolor( np.linspace(-1,1,np.size(LANEXm[0,alim1:alim2])), Positions, LANEXm[:,alim1:alim2], cmap=custom_cmap)
#                plt.title('CEP =' +str(param_uniq[i1]) +' rad')
#                plt.ylim(ymin=0)
#                plt.colorbar()
#                plt.savefig(r'G:\THESE_TRAVAIL\CEP_espectro\\i1_cb.png', bbox_inches='tight')
#                
#                
#                plt.figure(figsize=(1,4))
#                plt.pcolor( np.linspace(-1,1,np.size(LANEXm[0,alim1:alim2])), Positions[idx0::], LANEXm[idx0::,alim1:alim2], cmap=custom_cmap)
#                #plt.pcolor( np.linspace(-1,1,np.size(LANEXm[0,alim1:alim2])), E_interp, LANEXm[idx0::,alim1:alim2], cmap=custom_cmap)
#                plt.title('CEP =' +str(param_uniq[i1]) +' rad')
#                plt.ylim(ymin=0, ymax=np.max(Positions))
#                plt.ylabel('Position (mm)')
#                plt.xlabel('Position (a.u.)')
#                plt.yticks((0,10,20,30,40))
##                plt.colorbar()
##                
#                plt.twinx()
#                plt.ylim(ymin=0, ymax=np.max(Positions))
#                plt.yticks((Positions[idx0::][0],Positions[idx0::][50], Positions[idx0::][100],Positions[idx0::][150],Positions[idx0::][200],Positions[idx0::][250]), (str(np.round(E_interp[0],2)),str(np.round(E_interp[50],2)),str(np.round(E_interp[100],2)),str(np.round(E_interp[150],2)), str(np.round(E_interp[200],2)), str(np.round(E_interp[250],2))  ))
#                #plt.yticklabels( (str(E_interp[50]),str(E_interp[100]),str(E_interp[150])) )
#                #plt.twinx()
#                #plt.plot(Positions[idx0::], E_interp)
##                
##                
#                plt.ylabel('E (Mev)')
#                plt.savefig(r'G:\THESE_TRAVAIL\CEP_espectro\\i1.png', bbox_inches='tight')
#                
#                
#                
#                
#            if j==i2:
#                plt.figure(figsize=(5,4))
#                plt.pcolor( np.linspace(-1,1,np.size(LANEXm[0,alim1:alim2])), Positions, LANEXm[:,alim1:alim2], cmap=custom_cmap)
#                plt.title('CEP =' +str(param_uniq[i2]) +' rad')
#                plt.ylim(ymin=0)
#                plt.colorbar()
#                plt.savefig(r'G:\THESE_TRAVAIL\CEP_espectro\\i2_cb.png', bbox_inches='tight')
#
#                plt.figure(figsize=(1,4))
#                plt.pcolor( np.linspace(-1,1,np.size(LANEXm[0,alim1:alim2])), Positions[idx0::], LANEXm[idx0::,alim1:alim2], cmap=custom_cmap)
#                #plt.pcolor( np.linspace(-1,1,np.size(LANEXm[0,alim1:alim2])), E_interp, LANEXm[idx0::,alim1:alim2], cmap=custom_cmap)
#                plt.title('CEP =' +str(param_uniq[i2]) +' rad')
#                plt.ylim(ymin=0, ymax=np.max(Positions))
#                plt.ylabel('Position (mm)')
#                plt.xlabel('Position (a.u.)')
#                plt.yticks((0,10,20,30,40))
##                plt.colorbar()
##                
#                plt.twinx()
#                plt.ylim(ymin=0, ymax=np.max(Positions))
#                plt.yticks((Positions[idx0::][0],Positions[idx0::][50], Positions[idx0::][100],Positions[idx0::][150],Positions[idx0::][200],Positions[idx0::][250]), (str(np.round(E_interp[0],2)),str(np.round(E_interp[50],2)),str(np.round(E_interp[100],2)),str(np.round(E_interp[150],2)), str(np.round(E_interp[200],2)), str(np.round(E_interp[250],2))  ))
#                #plt.yticklabels( (str(E_interp[50]),str(E_interp[100]),str(E_interp[150])) )
#                #plt.twinx()
#                #plt.plot(Positions[idx0::], E_interp)
#                
#                
#                plt.ylabel('E (Mev)')
#                plt.savefig(r'G:\THESE_TRAVAIL\CEP_espectro\\i2.png', bbox_inches='tight')
#                
#                
#
##○TEST pinhole a executer dans la console
#SpaElectronsFile = np.array(  Image.open (SpaElectronsFolder + 'tir' + str(316) + '.tiff' ) )
#LANEXpx = SpaElectronsFile[yll:yur,xll:xur]
#LANEXpx = LANEXpx[0:-1:dsz, 0:-1:dsz]
#plt.pcolor(np.rad2deg(temp_phi),np.rad2deg(temp_theta), LANEXpx)
#


                 ## --------------------------------- ##   


    elif espectro == 0:  
        #cartesian coordinates x,y,z : no scaling needed here for LANEXm
        x = A*np.cos(alpha)                #along the target plane
        y = d - A*np.sin(alpha)            #along the target normal
        z = B                              #vertical       
    
        #spherical 
        r=np.zeros(np.shape(x))
        theta=np.zeros(np.shape(x))
        phi=np.zeros(np.shape(x))   
        for n in range(len(x[:,0])):
            for m in range(len(x[0,:])):
                r[n,m] = np.sqrt( (x[n,m])**2 + (y[n,m])**2 +(z[n,m])**2 ) 
                theta[n,m] = np.arcsin(z[n,m]/r[n,m])
                phi[n,m] = np.arctan(x[n,m]/y[n,m])
        temp_phi =  np.linspace( phi.min(), phi.max(),   np.size(px_temp)  )
        temp_theta = np.linspace(theta.min(), theta.max(),    np.size(py_temp) ) 
        mesh_phi, mesh_theta = np.meshgrid(temp_phi, temp_theta)


        # START LOOP HERE
        for j, k in enumerate (fileNumbers): 
        
            # load file 
            SpaElectronsFile = np.array(  Image.open (SpaElectronsFolder + 'tir' + str(k) + '.tiff' ) )
    
            #x, y, and values in px
            LANEXpx = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px   
        
            #we need to get rid of hot pixels
            temp = sorted(LANEXpx.flatten())[::-1]    #array -> list, then organize from highest to lowest value
            n_idx=[]
            for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
                if temp[i]-temp[i+1] >= 100:
                    n_idx=np.append(n_idx,i)
            if np.size(n_idx) > 0  :
                idx = np.where(LANEXpx >= temp[int(n_idx.max())])
                for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
                    LANEXpx[idx[0][j], idx[1][j]] = LANEXpx[idx[0][j]-1, idx[1][j]]  
            
            #downsize ie take one every dsz points
            LANEXpx = LANEXpx[0:-1:dsz, 0:-1:dsz]
        
            #with real dimensions in the LANEX plane, origin = intersection with target normal
            LANEXm = LANEXpx / (s**2)  #now in Coulomb/m²    
        
            #spherical 
            LANEXsph=np.zeros(np.shape(LANEXm))
            for n in range(len(x[:,0])):
                for m in range(len(x[0,:])):
                    LANEXsph[n,m]=LANEXm[n,m]*r[n,m]**2*np.cos(theta[n,m])  #Coulomb/sr                 
            LANEXsph2 = sc.interpolate.griddata((phi.flatten(), theta.flatten()), LANEXsph.flatten(), (mesh_phi, mesh_theta), method = 'nearest') 
            

            print(k)   #just to see the progression when you run this script
            ElecScan.append(list(sum(LANEXsph2)))
        
            #flip 90degree to get the delays on the horizontal axis, more convenient
            ElecScan_r = np.zeros ((np.shape(ElecScan)[1],np.shape(ElecScan)[0]))
            for k, osef in enumerate (np.zeros(np.shape(ElecScan)[0])):
                for j, osef_aussi in enumerate (np.zeros(np.shape(ElecScan)[1])):
                    ElecScan_r[j,k] = ElecScan[k][j]
           
            

            
            
    #average : 
    Belec = np.zeros((np.shape(ElecScan_r)[0],np.size(param_uniq)))
    for j, valuej in enumerate (param_uniq):
        Aelec = np.zeros((np.shape(ElecScan_r)[0],nb_param[j]))
        count = 0
        for i, value in enumerate (fileNumbers):
            if Params [i] == valuej :
                Aelec[:,count] = ElecScan_r[:,i]
                if count < nb_param[j] :
                    count = count +1 
        Belec[:,j]=np.sum(Aelec,1)/nb_param[j]     #average 
    
          
             

             
             #%%
    
#   initial script     
#        # load data
#        SpaElectronsFolder = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\electrons\\' 
#        SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(k) + '.tiff' ) 
#        SpaElectronsFile = np.array(SpaElectronsFile)
#
#        #x, y, and values in px
#        px_temp = np.linspace(0,int(xur-xll)-1,xur-xll)
#        py_temp = np.linspace(0,yur-yll-1,yur-yll)
#        LANEXpx = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px   
#
#        #we need to get rid of hot pixels
#        temp = sorted(LANEXpx.flatten())[::-1]    #array -> list, then organize from highest to lowest value
#        n_idx=[]
#        for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
#            if temp[i]-temp[i+1] >= 100:
#                n_idx=np.append(n_idx,i)
#        idx = np.where(LANEXpx >= temp[int(n_idx.max())])
#        for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
#            LANEXpx[idx[0][j], idx[1][j]] = LANEXpx[idx[0][j]-1, idx[1][j]]  
#        
#
#        if espectro == 0:
#        #test : downsize ie take one every dsz points
#            dsz = 2
#            LANEXpx = LANEXpx[0:-1:dsz, 0:-1:dsz]
#            px_temp = px_temp[0:-1:dsz]
#            py_temp = py_temp[0:-1:dsz]
#                    
#        #what's the convolution for?
##        dwnsmplby = 2
##        foo = np.ones((dwnsmplby,dwnsmplby), dtype=int)/dwnsmplby**2
##        LANEXpx = sc.signal.convolve2d(LANEXpx, foo,mode='valid', boundary='fill', fillvalue=0)
#
#        
#         #create coordinate mesh :
#        [px,py]=np.meshgrid(px_temp,py_temp)
#            
#                #            #plot if you want to check the ROI
#                #            fig_es = plt.figure(figsize=(6,5.6))
#                #            plt.pcolor(LANEXpx, cmap = 'inferno')
#                #            plt.colorbar()
#                #            plt.xlabel('x (pixels)') ; plt.ylabel('y (pixels)')     
#                    
#    
#        #with real dimensions in the LANEX plane, origin = intersection with target normal
#        A = a1+px*s
#        B = -(yur-yll)*s/2 + py*s
#        LANEXm = LANEXpx / (s**2)  #now in Coulomb/m²    
#        
#                #            plt.figure(figsize=(6,5.6))
#                #            plt.pcolor(a,b,LANEXm, cmap = 'viridis')
#                #            plt.xlabel('a (m)')
#                #            plt.ylabel('b (m)')    
#                
#    
#        #cartesian coordinates x,y,z : no scaling needed here for LANEXm
#        x = A*np.cos(alpha)                #along the target plane
#        y = d - A*np.sin(alpha)            #along the target normal
#        z = B                              #vertical       
#    
#    
#        #spherical 
#        r=np.zeros(np.shape(x))
#        theta=np.zeros(np.shape(x))
#        LANEXsph=np.zeros(np.shape(LANEXm))
#        phi=np.zeros(np.shape(x))
#        for n in range(len(x[:,0])):
#            for m in range(len(x[0,:])):
#                r[n,m] = np.sqrt( (x[n,m])**2 + (y[n,m])**2 +(z[n,m])**2 ) 
#                theta[n,m] = np.arcsin(z[n,m]/r[n,m])
#                phi[n,m] = np.arctan(x[n,m]/y[n,m])
#                LANEXsph[n,m]=LANEXm[n,m]*r[n,m]**2*np.cos(theta[n,m])  #Coulomb/sr
#                
#                        
#        temp_phi =  np.linspace( phi.min(), phi.max(),   np.size(px_temp)  )
#        temp_theta = np.linspace(theta.min(), theta.max(),    np.size(py_temp) ) 
#        mesh_phi, mesh_theta = np.meshgrid(temp_phi, temp_theta)
#        
#        LANEXsph2 = sc.interpolate.griddata((phi.flatten(), theta.flatten()), LANEXsph.flatten(), (mesh_phi, mesh_theta), method = 'nearest') 
#    
#    #    plt.figure(figsize=(6,5.6))
#    #    plt.pcolor(np.rad2deg(temp_phi),np.rad2deg(temp_theta),LANEXsph2, cmap = 'inferno')
#    #    plt.xlabel('Azimuth angle $\phi$ (deg)')
#    #    plt.ylabel('Polar angle' +  r' ${\Theta}$' + ' (deg)')
#    #    plt.colorbar()
#    
#    
#    #    plt.figure()
#    #    plt.plot(np.rad2deg(temp_phi), sum( LANEXsph2))
#    
#        print(k)
#        ElecScan.append(list(sum(LANEXsph2)))
#        
#            #flip 90degree to get the delays on the horizontal axis, more convenient
#    ElecScan_r = np.zeros ((np.shape(ElecScan)[1],np.shape(ElecScan)[0]))
#    for k, osef in enumerate (np.zeros(np.shape(ElecScan)[0])):
#        for j, osef_aussi in enumerate (np.zeros(np.shape(ElecScan)[1])):
#            ElecScan_r[j,k] = ElecScan[k][j]
#    
    
    #
    #    fig_e = plt.figure(figsize=(9,3))
    #    fig_e.suptitle ( 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=10 )
    #    plt.pcolor(np.append(ScanDelays,ScanDelays[-1]+1), np.rad2deg(temp_phi), ElecScan_r/np.max(ElecScan_r), cmap = 'inferno')
    #    plt.ylabel('$\phi_{elec}$ (deg)')
    #    plt.xlabel('Prepulse delay (ps)')
    #
    #
    
    
    
    
    
    
    
    
    
#%%



    







#%% ions 
ions_scan = []
if ions==1:
    for j, k in enumerate (fileNumbers):
        ionsfilepath = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\ions\\tir' + str(k) + '_spectrum_new.txt' #_new.txt' 
        if Path(ionsfilepath).is_file()== True :  
            E = np.loadtxt(ionsfilepath, skiprows = (0), usecols=(0,) )*1e3
            dNdE = np.loadtxt(ionsfilepath, skiprows = (0), usecols=(1,) )
            ions_scan.append(dNdE)
        else :
            if ions_scan:
                ions_scan.append ( np.zeros(np.size(ions_scan[0])) )
            else :
                ions_scan.append ( np.zeros(300)) 

        
            #flip 90degree to get the delays on the horizontal axis, more convenient
    IonsScan_r = np.zeros ((np.shape(ions_scan)[1],np.shape(ions_scan)[0]))
    for k, osef in enumerate (np.zeros(np.shape(ions_scan)[0])):
        for j, osef_aussi in enumerate (np.zeros(np.shape(ions_scan)[1])):
            IonsScan_r[j,k] = ions_scan[k][j]


        
    #average : 
    Bions = np.zeros((np.shape(IonsScan_r)[0],np.size(param_uniq)))
    for j, valuej in enumerate (param_uniq):
        Aions = np.zeros((np.shape(IonsScan_r)[0],nb_param[j]))
        count = 0
        for i, value in enumerate (fileNumbers):
            if Params [i] == valuej :
                Aions[:,count] = IonsScan_r[:,i]
                if count < nb_param[j] :
                    count = count +1 
        Bions[:,j]=np.sum(Aions,1)/nb_param[j]     #average 













#%%%smooth e spectro

#textsize = 14
#
#
#
#if espectro==1:
#    Ecutoff=[]
#    Belec_new = np.zeros((303,np.size(param_uniq)))
#    for v, value in enumerate(Belec[0,:]):
#        from scipy.interpolate import make_interp_spline
#        xnew = np.linspace(E_interp.min(), E_interp.max(), 303) 
#        spl = make_interp_spline(E_interp[::-1], Belec[:,v][::-1]/np.max(Belec), k=3)
#        y_smooth = spl(xnew)
##        
##        plt.plot(E_interp[::-1],  Belec[:,v][::-1]/np.max(Belec))
##        plt.plot(xnew, y_smooth, color='red')
##        plt.ylim(0,0.5)
##        plt.axhline(0.01, color='green')
#        
#        kkk = np.min( np.where(np.round(y_smooth,3) <= 0.03))
##        plt.axvline(xnew[kkk])
#        
#        Ecutoff = np.append(Ecutoff,xnew[kkk] )
#        Belec_new[:,v]=y_smooth
#    figel=plt.figure(figsize=(8,2))           
#    ax1e = figel.add_subplot(111)
#    
#    
#    pu_int=np.zeros(np.size(param_uniq))
#    lax_p = param_uniq[-1] + (param_uniq[-1] - param_uniq[-2])
#    param_uniq_n = np.append(param_uniq, lax_p)
#    for u, v in enumerate(param_uniq):
#        if u<np.size(param_uniq_n):
#            pu_int[u]=param_uniq_n[u] + (param_uniq_n[u+1]-param_uniq_n[u])/2
#           # pu_int[u]=v + (param_uniq_n[u+1]-param_uniq_n[u])/2
#    #pu_int=np.append(pu_int[1::], param_uniq[-1] + (param_uniq[-1]-param_uniq[-2])/2)
#    
#    #ecm=ax1e.pcolor(param_uniq,xnew[::-1],Belec_new[::-1], cmap=custom_cmap)
#    #ecm=ax1e.pcolor(pu_int,xnew[::-1],Belec_new[::-1], cmap=custom_cmap)
#    #
#    ecm=ax1e.pcolor(param_uniq_n,xnew[::-1],Belec_new[::-1], cmap=custom_cmap,vmin=0,vmax=1)
#    
#    
#    
#    
#    #ax1e.set_xticks([0,105,205,305,405])
#    #ax1e.set_xticklabels( ['0','100','200','300','400'])
#    #ecm=ax1e.pcolor(pu_int,xnew[::-1],Belec_new[::-1], cmap=custom_cmap)
#    ax1e.set_xticks([-2.75,-1.75,-0.75,0.25,1.25,2.25,3.25])
#    ax1e.set_xticklabels( ['-3','-2','-1','0','1','2','3'])
#    #ax1e.set_xlim(0,15)
#    
#    
#    ax2e = ax1e.twinx()
#    #ax2e.set_xticks([-2.75,-1.75,-0.75,0.25,1.25,2.25,3.25], ['-3','-2','-1','0','1','2','3'])
#    
#    pu_int_bis=np.zeros(np.size(pu_int))
#    lax_p = pu_int[-1] + (pu_int[-1] - pu_int[-2])/2
#    param_uniq_nn = np.append(pu_int, lax_p)
#    for u, v in enumerate(param_uniq_nn):
#        if u<np.size(pu_int):
#            pu_int_bis[u]=v + (param_uniq_nn[u+1]-param_uniq_nn[u])/2
#    
#    
#    
#    #scat = ax2e.scatter(param_uniq,Ecutoff, color='blue', linewidth=0.2)
#    
#    scat = ax2e.scatter(pu_int,Ecutoff, color='blue', linewidth=0.2)
#    ax2e.tick_params(axis='y', labelright=False)
#    ax2e.set_ylim(np.min(E_interp),2.1)
#    #ax1e.set_xlim(0,455)
#    ax1e.set_ylim(np.min(E_interp),2.1)
#    ax1e.set_yticks((0.5,1,1.5,2))
#
#    ax1e.set_ylabel('E (MeV)', fontsize=textsize+1)
#    ax1e.set_xlabel('relative CEP (rad)', fontsize=textsize+1)
#    ax1e.tick_params(labelsize=textsize)
#    
#    
#    cax = ax1e.inset_axes([1.02, 0.1, 0.03, 0.8], transform=ax1e.transAxes)
#    cbar= figel.colorbar(ecm,cax=cax)
#    cbar.set_label('dN/dE (a.u.)', rotation=270+180)
#    cbar.ax.set_yticks([0.01, 0.5,0.99]) # ('0','0.5','1') )
#    figel.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'_ELECSPEC.png', dpi=500, bbox_inches = 'tight')

#%%  electron spectra 

#if espectro==1:
#
#    plt.figure(figsize=(4,3))
#
#    #plt.plot(E_interp, Belec[:,36]/np.max(Belec))
#    plt.plot(xnew[::-1],Belec_new[::-1][:,i1], color='navy', label=str(param_uniq[i1]) + 'rad')
#    plt.plot(xnew[::-1],Belec_new[::-1][:,i2], color='orange', label=str(param_uniq[i2]) + 'rad')
#    
#   # plt.plot(xnew[::-1],Belec_new[::-1][:,8], color='green', label='+1 rad')
#    
#    plt.axvline(Ecutoff[i1], linestyle='--', color='navy') #, label='cutoff')
#    plt.axvline(Ecutoff[i2], linestyle='--', color='orange') #, label='cutoff')
#    
#    plt.axhline(0.03,linestyle='--', color='red') #, label='threshold')
#    plt.legend(fontsize=textsize-2, loc='lower left')
#    #plt.title(r'$\tau$=190ps, $E_{max}$=1.6MeV', fontsize=textsize)
#    plt.xlim(xnew.min(),1.6)
#    plt.ylim(ymin=0)
#    plt.xlabel('E(MeV)', fontsize=textsize)
#    plt.ylabel('dN/dE (a.u)', fontsize=textsize)
#    
#    #plt.xticks((0.4,1,1.6,2.2))
#    plt.xticks((0.5,1,1.5,2), ('0.5','1','1.5','2'))
#    plt.tick_params(labelsize=textsize)
#    plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'_elec_spec_2values.png', dpi=500, bbox_inches = 'tight')
#    
#    #PB colorbar a cause du plt scatter
    
#%% corresponding HHG plot
##
#    plt.figure(figsize=(4,1.5))
#    plt.plot(Ev[::-1], Bharm[:,i1]/1, color='navy', label=str(param_uniq[i1]) + 'rad')
#    plt.plot(Ev[::-1], Bharm[:,i2]/1, color='orange', label=str(param_uniq[i2]) + 'rad')
#    plt.xlim(10.5,32)
#    plt.ylim(ymin=0)
#    plt.xlabel('Photon energy (eV)', fontsize=textsize+1)
#    plt.tick_params(labelsize=textsize)
#    plt.yticks((0,0.5e7,1e7))
#    plt.ylabel('Intensity (a.u.)', fontsize=textsize+1)
#    plt.legend(fontsize=textsize-1, loc='upper right')
#    plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'_HHG_spec_2values.png', dpi=500, bbox_inches = 'tight')
#
###%%
##
##param_uniq = param_uniq+396

#%%    Final figure
fig = plt.figure(figsize=(7,7)) #(4,7))   #10,10  or#14,10 pour les allongés
textsize = 15
fig.text (0.5,1, 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=textsize -6)


param_uniqb=np.zeros(np.size(param_uniq)+1)
bigp=np.append(param_uniq, param_uniq[-1]+(param_uniq[-1]-param_uniq[-2]))
bigp=np.append(param_uniq[0] - (param_uniq[1]-param_uniq[0]), bigp)
for i, value in enumerate (bigp) :
    if i <= (np.size(param_uniq)):
        param_uniqb[i]= value + (bigp[i+1]-value)/2
param_uniqb=np.append(param_uniqb, param_uniq[-1]+(param_uniq[-1]-param_uniq[-2])/2)


#plt.pcolor(param_uniqb[:], Ev[::-1], Bharm[:,:],cmap=custom_cmap, vmin=0,vmax=np.max(Bharm[1:-1]))
#plt.ylim(11,25)
#plt.xlim(0,20)

if scan_param == 'Gradient':
##n=31
    n=np.where(param_uniq==10)[0][0] +1
    xmaxlim=bigp[n]
              
              
          
#plt.pcolor(bigp[0:n+1], Ev[::-1], Bharm[:,0:n],cmap=custom_cmap, vmin=0,vmax=np.max(Bharm[1:-1]))
#plt.pcolor(param_uniqb[0:n+1], Ev[::-1], Bharm[:,0:n],cmap=custom_cmap, vmin=0,vmax=np.max(Bharm[1:-1]))
#plt.pcolor(param_uniqb[:], Ev[::-1], Bharm[:,:],cmap=custom_cmap, vmin=0,vmax=np.max(Bharm[1:-1]))
##plt.xlim(bigp[0], bigp[5])
#plt.xlim(0,bigp[n])

#xmaxlim=np.max((np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])))
xmaxlim=20


##xmaxlim=np.max(param_uniq)
#
#
##
#im = plt.imshow(Bharm, aspect=1/20)
#
#ax = plt.gca();
##
### Major ticks
##ax.set_xticks(np.arange(0, 12, 1))
##ax.set_yticks(np.arange(11, 28, 1))


#xticklab=np.arange(0,14,2)
#idxs=[np.where(param_uniq==xticklab[i]) for i in range(np.size(xticklab))]
#tickpos=[ param_uniq[idxs[x][0][0]]+(param_uniq[idxs[x][0][0]+1]-param_uniq[idxs[x][0][0]])/2 for x in range(np.size(idxs)) ]
#
#

if harm==1:
    #harmonics
    gs = gridspec.GridSpec(3, 1, left=0.06, top=0.88, height_ratios=[3,1,1],hspace=0.3, figure=fig)
    fig.add_subplot(gs[0, 0])
    #plt.pcolor( (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])), Ev[::-1], Bharm/np.max(Bharm), cmap = 'jet',vmin=0,vmax=1)
    #plt.pcolor( (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])), Ev[::-1], Bharm/1,  cmap='inferno',vmin=0,vmax=np.max(Bharm[1:-1]))
    #plt.pcolor(param_uniq,   Ev[::-1], Bharm/1,  cmap='inferno',vmin=0,vmax=np.max(Bharm[1:-1]))
#    plt.pcolor(param_uniqb ,   Ev[::-1], Bharm/1,  cmap='magma_r',vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1]))
  # )
    
    #plt.pcolor( param_uniq, Ev[::-1], Bharm/1, cmap = custom_cmap,vmin=0,vmax=np.max(Bharm[1:-1]))
#    plt.xlim(0,12)

    if scan_param == 'CEP' :
        #plt.pcolor(    (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])) ,   Ev[::-1], Bharm/1,  cmap='magma_r',vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1]))
        plt.pcolor(    (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])) ,   Ev[::-1], np.transpose( [Bharm[:,i]/np.max(Bharm[:,i]) for i, j in enumerate (Bharm[0,:]) ]),  cmap='magma_r', vmin=0, vmax=1)#vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1]))
        plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
        plt.xlabel('relative CEP (rad)', fontsize=textsize )

        #plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
        #plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
        plt.xlim(-3,3.5)
        plt.tick_params(axis='both', labelsize=textsize)
    elif scan_param == 'Gradient' :
        # add plasma gradient scale length ticks on top
        plt.pcolor(param_uniqb[0:-1] ,   Ev[::-1], Bharm/1,  cmap=custom_cmap,vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1])) #custom_cmap
        #plt.pcolor((np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])) ,   Ev[::-1], Bharm/1,  cmap='magma_r',vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1]))  
        plt.tick_params(axis='both', labelsize=textsize)   
        #plt.xticks([elt +(param_uniq[idx+1]-elt)/2 for idx, elt in enumerate(param_uniq[0:-2]) if np.round(elt,1)%2==0],[int(elt) for idx, elt in enumerate(param_uniq[0:-2]) if np.round(elt,1)%2==0])
     #   plt.xticks([elt +(param_uniq[0]-param_uniq[1])/-2 for idx, elt in enumerate(param_uniq) if np.round(elt,1)%2==0],[int(elt) for idx, elt in enumerate(param_uniq) if np.round(elt,1)%2==0])
       # plt.xticks(tickpos,xticklab )

     #plt.xticks(np.arange(0.5,14.5,2), np.arange(0,14,2))
#        plt.xticks([elt +(param_uniq[idx+1]-elt)/2 for idx, elt in enumerate(param_uniq[0:-2]) if np.round(elt,1)%100==0],[int(elt) for idx, elt in enumerate(param_uniq[0:-2]) if np.round(elt,1)%100==0])
        plt.xlim(0,xmaxlim)
         
#      #  plt.xticks([elt+0.5 for idx, elt in enumerate(param_uniq) if np.round(elt,1) % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if np.round(elt,1) % 2 == 0])
#        #plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
##        plt.xlim(0,xmaxlim)
#        plt.twiny()
#        #plt.xticks([elt  for idx  , elt in enumerate(param_uniq) if idx%2 == 0])
#        #plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(np.round(int(elt)*plasma_speed)) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
#        plt.xticks(  np.linspace(0,xmaxlim,5)+(param_uniq[1]-param_uniq[0])/2   , ( x for x in (np.round( np.linspace(0,xmaxlim,5)*plasma_speed/800,3) )))# , ( [float(str(te[x])[5:7])*plasma_speed for x,k in enumerate (te)] ))
#      #  plt.xticks([(elt+((param_uniq[-1]-param_uniq[-2])/2))*plasma_speed/800 for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt)*plasma_speed/800 for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
#        #plt.xlim(0,12)
#      #  plt.xticks(  np.linspace(0,xmaxlim,5)   , ( x for x in (np.round( np.linspace(0,xmaxlim,5)*plasma_speed/800,3) )))# , ( [float(str(te[x])[5:7])*plasma_speed for x,k in enumerate (te)] ))
##        plt.xlim(0,12)
#        plt.xlim(0,xmaxlim)
#        plt.xlabel('$(L_g-L_0)/ \lambda$', labelpad=15, fontsize=textsize)
    elif scan_param == 'Parabola_pos':
        plt.xticks([elt for idx, elt in enumerate(param_uniq) if idx % 2 == 1], [int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 1])
        plt.pcolor(param_uniqb[0:-1] ,   Ev[::-1], Bharm/1,  cmap=custom_cmap,vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1])) #custom_cmap
    plt.ylabel('Photon energy (eV)', fontsize=textsize )
    plt.axhline(y=30, color='red', ls ='--',linewidth='1', label='$E_{p,max}$')
    plt.tick_params(axis='both', labelsize=textsize)
    plt.ylim(ymin=10.8, ymax=35) #35 max
    #♠plt.yticks((12,15,18,21,24))
            
            
    #plt.legend(fontsize=textsize-2, loc='upper right')
    if elec ==0 and ions==0 :
        if scan_param == 'CEP' :
            plt.xlabel('CEP (rad)',  fontsize=textsize)
        elif scan_param == 'Gradient' :
            plt.xlabel('Prepulse lead (ps)', fontsize=textsize)
        elif scan_param == 'Parabola_pos':
            plt.xlabel('OAP distance (µm)',  fontsize=textsize)         
    def ter(lg):
        return lg/plasma_speed
        
    if scan_param == 'Gradient':
        ax2=plt.twiny()
        ax2.set_xlim(0,xmaxlim)
 #       ax2.set_xticks((ter(0),ter(0.02*800),ter(0.04*800), ter(0.06*800),ter(0.08*800),ter(0.1*800)))#(which='major', [ter(i) for i in [10,9,8,7,6,5,4,3,2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1] ]) # , (te(ter(1), ter(0.5), ter(0.3))) )
#        ax2.set_xticklabels(('0', '0.02','0.04','0.06', '0.08','0.1'), fontsize=textsize)            
        grad_list = np.round(np.arange(0,np.round(plasma_speed*xmaxlim/800,1),0.1),1) #[0,0.1,0.2,0.3] #,0.4,0.5] #,0.5,0.6]
        #grad_list = np.round(np.arange(0,np.round(plasma_speed*xmaxlim/800,1)+0.1,1),1) #[0,0.1,0.2,0.3] #,0.4,0.5] #,0.5,0.6]    
        ax2.set_xticks([ter(x*800) for i, x in enumerate(grad_list)]) #(ter(0),ter(0.04*800),ter(0.08*800), ter(0.12*800),ter(0.16*800),ter(0.2*800),ter(0.24*800),ter(0.28*800)))#(which='major', [ter(i) for i in [10,9,8,7,6,5,4,3,2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1] ]) # , (te(ter(1), ter(0.5), ter(0.3))) )    
        #ax2.set_xticklabels([str(int(x)) for i, x in enumerate(grad_list)], fontsize=textsize) 
        ax2.set_xticklabels([str(x) for i, x in enumerate(grad_list)], fontsize=textsize) 

     #   ax2.set_xlim(0,plasma_speed*xmaxlim/800 +0.1)

        ax2.set_xlabel(r'$(L_g - L_0)/\lambda$', fontsize=textsize)               

  #  cbaxesh = fig.add_axes([0.915, 0.472, 0.03, 0.408 ]) #left, bottom, width, height 
    cbaxesh = fig.add_axes([0.95, 0.55, 0.02, 0.26 ]) #left, bottom, width, height 
    cbaxesh = plt.colorbar(cax=cbaxesh,orientation='vertical')
    cbaxesh.formatter.set_powerlimits((0, 0))
    cbaxesh.ax.yaxis.get_offset_text().set_fontsize(textsize-4)
    cbaxesh.ax.tick_params(labelsize=textsize-2)
    #cbaxesh.set_ticks((0,0.2,0.4,0.6,0.8,1), (('0','0.2','0.4','0.6','0.8','1')))
    cbaxesh = plt.ylabel('Spectral Intensity (a.u)',  fontsize=textsize-2 )
    cbaxesh =  plt.tick_params(labelsize=textsize-2)
    #cbaxesh =  plt.yticks
    
#    cbaxesh.formatter.set_powerlimits((0, 0))
#
#    cbaxesh.ax.yaxis.get_offset_text().set_fontsize(35)
    
    
    if elec==1:
        #electrons
        fig.add_subplot(gs[1, 0])  
        if espectro == 0 :
            if scan_param=='Gradient' or 'Parabola_pos':
            #plt.pcolor((np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])), np.rad2deg(temp_phi),  Belec, cmap = 'viridis') #, vmin=0, vmax=1)
                plt.pcolor(param_uniqb[0:-1], np.rad2deg(temp_phi),  Belec/1, cmap = 'viridis', vmin=0, vmax=np.max(Belec))
           # (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2]))
           
            elif scan_param=='CEP':
                plt.pcolor((np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])), np.rad2deg(temp_phi),  Belec, cmap = 'viridis') #, vmin=0, vmax=1)

           
            #plt.axhline(39, linestyle='--', color='black')
            plt.ylabel('$\phi_{az}$ (deg)',  fontsize=textsize)
#            plt.xlim(0,12)
#            plt.ylim(30,48.1)
            plt.yticks(np.arange(np.round(np.min(np.rad2deg(temp_phi)),0)+4, np.max(np.rad2deg(temp_phi)),5))
          #  plt.yticks((35,40,45,50))
            #plt.yticks((35,40,45,50))
        elif espectro == 1 :
            #plt.pcolor((np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])), E_interp,  Belec/np.max(Belec), cmap = custom_cmap, vmin=0, vmax=0.8)
            plt.pcolor((np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])),xnew[::-1],Belec_new[::-1], cmap=custom_cmap, vmin=0, vmax=np.max(Belec_new))
            plt.ylabel('$E_{electrons}$ (MeV)',  fontsize=textsize)
            #plt.ylim(ymin=E_interp[-1], ymax=1.8)
            plt.ylim(ymin=0.5, ymax=1.8)
        if scan_param == 'CEP' :
            #plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
            plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])

        
        elif scan_param == 'Gradient':
            #plt.xticks([elt +(param_uniq[idx+1]-elt)/2 for idx, elt in enumerate(param_uniq[0:-2]) if np.round(elt,1)%100==0],[int(elt) for idx, elt in enumerate(param_uniq[0:-2]) if np.round(elt,1)%100==0])
            plt.xlim(0,xmaxlim)
          #  plt.xticks([elt +(param_uniq[0]-param_uniq[1])/-2 for idx, elt in enumerate(param_uniq) if np.round(elt,1)%2==0],[int(elt) for idx, elt in enumerate(param_uniq) if np.round(elt,1)%2==0])
       #     plt.ylim(31,48)
        
            #plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
            #plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
           # plt.xticks([elt + (param_uniq[idx+1]-elt)/2 for idx, elt in enumerate(param_uniq) if np.round(elt,1) % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if np.round(elt,1) % 2 == 0])
            #plt.xlim(0,30)

        
        elif scan_param == 'Parabola_pos':
            plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0], [int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])

     
            #plt.xticks(ScanDelays+(ScanDelays[1]-ScanDelays[0])/2  , [round(ScanDelays[i],1) for i in range(len(ScanDelays))] )
            
        plt.tick_params(axis='both', labelsize=textsize)
        if ions==0:
            if scan_param == 'CEP':
                plt.xlabel('relative CEP (rad)'  ,  fontsize=textsize)
            elif scan_param == 'Gradient':
                plt.xlabel('Prepulse lead (ps)',  fontsize=textsize)
            elif scan_param == 'Parabola_pos':
                plt.xlabel('OAP distance (µm)',  fontsize=textsize)
##                    
#        ax2=plt.twiny()
#        ax2.set_xlim(0,xmaxlim)
#        ax2.set_xticks([ter(x*800) for i, x in enumerate(grad_list)]) #(ter(0),ter(0.04*800),ter(0.08*800), ter(0.12*800),ter(0.16*800),ter(0.2*800),ter(0.24*800),ter(0.28*800)))#(which='major', [ter(i) for i in [10,9,8,7,6,5,4,3,2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1] ]) # , (te(ter(1), ter(0.5), ter(0.3))) )    
#        ax2.set_xticklabels([str(int(x)) for i, x in enumerate(grad_list)], fontsize=textsize) 
##        ax2.set_xticks((ter(0),ter(0.5*800),ter(1*800), ter(1.5*800),ter(2*800),ter(2.5*800)))#(which='major', [ter(i) for i in [10,9,8,7,6,5,4,3,2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1] ]) # , (te(ter(1), ter(0.5), ter(0.3))) )
##        ax2.set_xticklabels(('0', '0.5','1','1.5', '2','2.5'),  fontsize=textsize-2 )            
#        ax2.set_xlabel(r'$(L_g - L_0)/\lambda$', fontsize=textsize)               
##                        
##                        
##                
#            
            
                
    
        #cbaxese = fig.add_axes([0.915, 0.292, 0.03, 0.135 ]) #left, bottom, width, height 
        cbaxese = fig.add_axes([0.95, 0.315, 0.02, 0.1 ]) #left, bottom, width, height 
        cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
        cbaxese.formatter.set_powerlimits((0, 0))
        cbaxese.ax.yaxis.get_offset_text().set_fontsize(textsize-4)
        cbaxese = plt.tick_params(labelsize=textsize-2)
        cbaxese = plt.ylabel('int. charge (a.u)',  fontsize=textsize-2 )

    
    if ions==1:
        #protons
        if elec==1:
            fig.add_subplot(gs[2,0])
        elif elec==0:
            fig.add_subplot(gs[1,0])
        #plt.pcolor(  np.append(ScanDelays,ScanDelays[-1]+(ScanDelays[-1] - ScanDelays[-2])  )   ,   E/1000.,   Bions/ np.max(Bions),  cmap='jet')
        plt.pcolor(param_uniqb[0:-1], E/1000.,   Bions/1,  cmap=custom_cmap, vmin=0, vmax= np.max(Bions))

#        plt.pcolor( (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])), E/1000.,   Bions/1,  cmap=custom_cmap, vmin=0, vmax= np.max(Bions))
#        plt.xlim(0,10)
        if scan_param == 'Parabola_pos' :
            plt.xlabel('Parabola position (µm)',  fontsize=textsize)
        elif scan_param == 'Gradient':
            plt.xlabel('Prepulse lead (ps)',  fontsize=textsize)
       #     plt.xticks([elt +(param_uniq[idx+1]-elt)/2 for idx, elt in enumerate(param_uniq[0:-2]) if np.round(elt,1)%2==0],[int(elt) for idx, elt in enumerate(param_uniq) if np.round(elt,1)%2==0])
            plt.xlim(0,xmaxlim)
          #  plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])

        elif scan_param == 'CEP':
            plt.xlabel('relative CEP (rad)',  fontsize=textsize)
            plt.xticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])

        plt.ylabel('$E_{protons}$ (MeV)'  , fontsize=textsize)
        plt.tick_params(axis='both', labelsize=textsize)
        plt.ylim(np.min(E)/1000,0.25)
        plt.yticks((0.1,0.15,0.2,0.25), ('0.1','','0.2',''))
        #plt.xticks(ScanDelays+(ScanDelays[1]-ScanDelays[0])/2  , [round(ScanDelays[i],1) for i in range(len(ScanDelays))] )
        if elec == 1:
            cbaxesi = fig.add_axes([0.95, 0.12, 0.02, 0.11 ]) #left, bottom, width, height 
        elif elec == 0 :
            cbaxesi = fig.add_axes([0.915, 0.292, 0.02, 0.135 ])
        cbaxesi = plt.colorbar(cax=cbaxesi,orientation='vertical')
        cbaxesi.formatter.set_powerlimits((0, 0))
        cbaxesi.ax.yaxis.get_offset_text().set_fontsize(textsize-4)
        cbaxesi = plt.ylabel('dN/dE (a.u)',  fontsize=textsize-2 )
        cbaxesi = plt.tick_params(labelsize=textsize-2)



        


if save == 1 :
    fig = plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'.png', dpi=500, bbox_inches = 'tight')

#fig = plt.savefig(r'F:\THESE_TRAVAIL\manuscrit_these_images\SHHG_data'+'\\'+scan_param+'Scan_' + str(dateCode) + '_full_shots' + str(file1) + '_' + str(file_last) +'.png', dpi=500, bbox_inches = 'tight')



#%% figure plus allongée papier CEP

figCEP=plt.figure(figsize=(5.3,2.6))
gs = gridspec.GridSpec(3, 3, left=0.06, top=0.88, height_ratios=[1.8,0.3,3],width_ratios=[3,0.2,0.7], hspace=0, wspace=0, figure=fig)

figCEP.add_subplot(gs[2, :])
plt.pcolor(    Ev[::-1],  (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])) ,  [Bharm[:,i]/np.max(Bharm[:,:]) for i, j in enumerate (Bharm[0,:]) ],  cmap='inferno', vmin=0, vmax=1)#vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1]))
#plt.pcolor(    Ev[::-1],  (np.append(param_uniq,param_uniq[-1]+param_uniq[-1]-param_uniq[-2])) ,  [Bharm[:,i]/np.max(Bharm[:,i]) for i, j in enumerate (Bharm[0,:]) ],  cmap='inferno', vmin=0, vmax=1)#vmin=Bharm[-1,-1],vmax=np.max(Bharm[1:-1]))
plt.yticks([elt+((param_uniq[-1]-param_uniq[-2])/2) for idx, elt in enumerate(param_uniq) if idx % 2 == 0],[int(elt) for idx, elt in enumerate(param_uniq) if idx % 2 == 0])
plt.xlim(10.9,24.5)
plt.ylabel('CEP (rad)')
plt.xlabel('photon energy (eV)')
plt.colorbar()

#figCEP.add_subplot(gs[2, 2])
#plt.scatter(  np.sum(Bharm,0)/np.max(np.sum(Bharm,0)),param_uniq, color='black')
#plt.plot(  np.sum(Bharm,0)/np.max(np.sum(Bharm,0)),param_uniq, color='black')
#plt.xlim(0.4,1.1)
#plt.ylim(-3.25,3.25)
#plt.xticks((0.4,0.7,1), ('0.4','0.7','1'))
#plt.yticks( (-3,-2,-1,0,1,2,3) )
#plt.tick_params(labelleft=False, labelbottom=False, labeltop=True, bottom=False, top=True)
#plt.xlabel('harm yield (a.u.)')

figCEP.add_subplot(gs[0, 0])
###"-1ps
#plt.plot(Ev[::-1], Bharm[:,4]/np.max(Bharm[:,4]), color='red', linewidth='2',alpha=1, label='$\Phi$=-1 rad')
#plt.plot(Ev[::-1]-0, Bharm[:,10]/np.max(Bharm[:,10])*1, color='navy', linewidth='2', alpha=1, label='$\Phi$=2 rad')
###-2ps
#plt.plot(Ev[::-1], Bharm[:,6]/np.max(Bharm[:,6]), color='red', linewidth='2',alpha=1, label='$\Phi$=0 rad')
#plt.plot(Ev[::-1]-0, Bharm[:,12]/np.max(Bharm[:,12])*1, color='navy', linewidth='2', alpha=1, label='$\Phi$=3 rad')
#-3ps
#☺plt.plot(Ev[::-1], Bharm[:,2]/np.max(Bharm[:,2]), color='red', linewidth='2',alpha=1, label='$\Phi$=-2 rad')
#plt.plot(Ev[::-1]-0, Bharm[:,8]/np.max(Bharm[:,8])*1, color='navy', linewidth='2', alpha=1, label='$\Phi$=1 rad')

plt.plot(Ev[::-1], Bharm[:,1]/np.max(Bharm[:,1]), color='red', linewidth='2',alpha=1, label='$\Phi$=-2.5 rad')
plt.plot(Ev[::-1]-0, Bharm[:,7]/np.max(Bharm[:,7])*1, color='navy', linewidth='2', alpha=1, label='$\Phi$=0.5 rad')


plt.xlim(10.95,24.5)
plt.yticks((0,0.2,0.4,0.6,0.8,1), ('0','0.2','0.4','0.6','0.8','1'))
plt.legend(loc='upper right', fontsize=10)
plt.ylabel('I (a.u)')
plt.ylim(0,1)
plt.tick_params(labelbottom=False)
plt.savefig(outdir+'\\'+scan_param+'CEPscan2ps_linear_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'.png', dpi=500, bbox_inches = 'tight')



#figCEP.add_subplot(gs[2, 2])



#
#save in a txt file
an_array= [Bharm[:,i]/np.max(Bharm[:,:]) for i, j in enumerate (Bharm[0,:]) ]
an_array2 = Ev[::-1]


#a_file = open("E_calib.txt", "w")
a_file = open("-1ps.txt", "w")
#np.savetxt('2ps prepulse delay (slice every 640 lines to obtain a 13 x 640 array)')
for row in an_array:
    np.savetxt(a_file, row)
#for i, j in enumerate(an_array2):
#    a_file.write(str(j)+ '\n')
#
a_file.close()
#
#








#%% 
textsize=18
plt.figure(figsize=(6.5,3))
plt.plot(Ev[::-1], Bharm[:,2]/np.max(Bharm[:,2]), color='red', linewidth='3',alpha=1, label='-2rad')
plt.plot(Ev[::-1]-0, Bharm[:,8]/np.max(Bharm[:,8])*1, color='navy', linewidth='3', alpha=1, label='1rad')
plt.xlim(10.95,22)
plt.legend(loc='upper right', fontsize=textsize)
plt.ylabel('I (a.u)', fontsize=textsize)
plt.xlabel('E (eV)', fontsize=textsize)
plt.tick_params(labelsize=textsize)

plt.ylim(0,1.1)
plt.yticks((0,0.5,1), ('0','0.5','1'))
#plt.title('prepulse lead time = 2ps')


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth




##################smooth:
#smoothed=smooth(Bharm[:,12]/np.max(Bharm[:,12]),15)
#plt.plot(Ev[::-1][88:-1], smoothed[88:-1], linestyle='--', linewidth='2', color='navy', alpha=1)
#
#
#smoothed2=smooth(Bharm[:,6]/np.max(Bharm[:,6]),15)
#plt.plot(Ev[::-1], smoothed2, linestyle='--', linewidth='2', color='red', alpha=1)
#


#############find extrema:
#import scipy
#from scipy import signal
#pp= scipy.signal.find_peaks(smoothed, height=0.24, threshold=None, distance=5, prominence=None, width=None, wlen=None, rel_height=0.5)
#pp2= scipy.signal.find_peaks(smoothed2, height=0.24, threshold=None, distance=5, prominence=None, width=None, wlen=None, rel_height=0.5)
#
#for i, j in enumerate (pp[0][0:8]):
#    if j:
#        if i==0:
#            print(j)
#            plt.scatter(Ev[::-1][85], Bharm[85,12]/np.max(Bharm[:,12]), marker='x',s=150, c='navy' )
#            plt.text(Ev[::-1][85]-0.1, Bharm[85,12]/np.max(Bharm[:,12])+0.12, str(round(Ev[::-1][82],1)) , color='navy' , fontsize=textsize-2)
#        
#        else:
#            plt.scatter(Ev[::-1][pp[0][0:8][i]], smoothed[pp[0][0:8][i]], marker='x',s=150, c='navy' )
#            plt.text(Ev[::-1][pp[0][0:8][i]]-0.1, smoothed[pp[0][0:8][i]]+0.06, str(round(Ev[::-1][pp[0][0:8][i]],1)) , color='navy' , fontsize=textsize-2)
#
#for i, j in enumerate (pp2[0][0:8]):
#    if j:
#        plt.scatter(Ev[::-1][pp2[0][0:8][i]], smoothed2[pp2[0][0:8][i]],marker='o', s=60, c='red' )
#        plt.text(Ev[::-1][pp2[0][0:8][i]]-0.1, smoothed2[pp2[0][0:8][i]]+0.06, str(round(Ev[::-1][pp2[0][0:8][i]],1)) , color='red' , fontsize=textsize-2)
#
#
#





#for j, k in enumerate (smoothed) :
#    if j > 0 and j<np.size(signs)-1:
#        if signs [j+1] != signs[j]:
#            sign_change.append(j)
#for n, p in enumerate (sign_change):
#    if n==0:
#        extrema_value.append(np.max(np.abs(smoothed[0:p]))   *   signs[sign_change[0]]  )       
#    else :
#        extrema_value.append(np.max(np.abs(smoothed[sign_change[n-1]:p]))  *   signs[sign_change[n]]   )
#extrema_idx = [np.where(smoothed == i) for i in extrema_value]
#extrema_delay = [  delay[ extrema_idx[i][0][0] ] for i, j in enumerate (extrema_idx) ]
#plt.scatter(extrema_delay, extrema_value, color='red',linewidth=s-10, label='extrema')
#
#




plt.tick_params(labelsize=textsize)
plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots___2ps__' + str(file1) + '_' + str(file_last) +'.png', dpi=500, bbox_inches = 'tight')






















#%%% TEST modulation depth
#from matplotlib.pyplot import cm
##
##plt.figure(figsize=(5,2))
#n=np.size(Bharm[0,0:5])
#colors=iter(cm.magma(np.linspace(0,1,n+1)))
#
##for i in range(n):
##    plt.plot(Ev[::-1],(Bharm[:,i]-Bharm[-1,i])/np.max(Bharm[:,0:n]),  c=c)
##    c=next(colors)
##plt.xlabel('Photon energy (eV)', fontsize= textsize)
##plt.ylabel('norm. I (a.u)', fontsize= textsize)
##plt.xlim(10.5,27)
##plt.ylim(0,1)
##plt.tick_params(labelsize= textsize-2)
##fig = plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'_overlapped_harm.png', dpi=500, bbox_inches = 'tight')
##
##
#
##c=colors
##for j in enumerate(Bharm[:,0:7]):
#for i in range(n):
##        if i<n:
#    c=next(colors)
#    #plt.plot(Ev[::-1],Bharm[:,i]/np.max(Bharm[:,i]),  c=c, label=str(i))  
#    plt.plot(Ev[::-1],Bharm[:,i],  c=c, label=str(param_uniq[i])+'ps')
#    #c=next(colors)
# #☻    c=next(colors)
#plt.xlabel('Photon energy (eV)', fontsize=textsize+4)
#plt.xlim(10,30)
#plt.ylim(ymin=0)
#plt.legend(fontsize=textsize+4)
#plt.tick_params(labelsize= textsize+4)
#
#fig = plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'_0to2_other_ps.png', dpi=500, bbox_inches = 'tight')

#%%
#plt.figure(figsize=(6,3))
#plt.semilogy( Ev[::-1],(Bharm[:,11]-Bharm[-1,11])/np.max(Bharm[:,11]), color='green', label=str(int(param_uniq[11]))+'µm')
#plt.semilogy( Ev[::-1],(Bharm[:,4]-Bharm[-1,4])/np.max(Bharm[:,4]), color='navy',label=str(int(param_uniq[4]))+'µm')
#plt.tick_params(labelsize= textsize)
#plt.legend(fontsize=textsize)
#plt.xlabel('E(eV)', fontsize=textsize)
#plt.ylabel('normalized intensity', fontsize=textsize)
#plt.title('Vertically integrated harmonic spectra', fontsize=textsize)
#plt.xlim(10.5,35)
##plt.ylim(0,1)
#plt.ylim(1.6e-2,1)
#         
#plt.figure(figsize=(5,2.1))
##plt.title('Vertically integrated harmonic spectra', fontsize=textsize)
#plt.plot( Ev[::-1],(Bharm[:,11]-Bharm[-1,11])/np.max(Bharm[:,11]), color='green', label=str(int(param_uniq[11]))+'µm')
#plt.plot( Ev[::-1],(Bharm[:,4]-Bharm[-1,4])/np.max(Bharm[:,4]), color='navy',label=str(int(param_uniq[4]))+'µm')
#plt.xlim(10.5,32)
#plt.ylim(0,1)
#plt.xlabel('E(eV)', fontsize=textsize)
#plt.ylabel('normalized intensity', fontsize=textsize)
#plt.tick_params(labelsize= textsize)
#plt.legend(fontsize=textsize)
#         
#if save == 1 :
#    fig = plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'_no_ions_12ps___least_mostmod.png', dpi=500, bbox_inches = 'tight')
#





#%%   integrated signal Z scan
param_uniq = param_uniq-8
#%%
textsize=13


fig_z, ax = plt.subplots(figsize=(6,3))
fig_z.subplots_adjust(right=0.75)

twin1 = ax.twinx()
twin2 = ax.twinx()
twin2.spines['right'].set_position(("axes", 1.25))


#fig_z=plt.figure(figsize=(6,3))
fig_z.text (0.2,0.94, 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=textsize -6)
ax.scatter(param_uniq, [ np.sum(Bharm[:,ee]) for ee in enumerate ( np.arange(0,np.shape(Bharm)[1])) ] , color='navy',label='HHG')
ax.set_ylabel('Harmonic signal (a.u.)', color='navy', fontsize=textsize)
ax.set_xlabel('OAP distance (µm)', color='black',  fontsize=textsize)
ax.tick_params(labelsize=textsize)
#ax.legend(fontsize=textsize, loc='upper left')
ax.tick_params(axis='y', labelcolor='navy',labelsize=textsize)
#ax.set_xticks((-10,-8,-6,-4,-2,0,2,4,6,8))
ax.set_xticks((-15,-10,-5,0,5,10))

twin1.scatter(param_uniq, [ np.sum(Belec[:,ee]) for ee in enumerate ( np.arange(0,np.shape(Belec)[1])) ] , marker='+', s=95, color='red', label='electrons')
twin1.set_ylabel('Electron signal (a.u.)', color='red',  fontsize=textsize)
twin1.tick_params(labelsize=textsize)
twin1.tick_params(axis='y', labelcolor='red',labelsize=textsize)
#






w0 =1.3e-6  #waist du faisceau à la foc (m)
l = 0.00005 #distance (m)
lambd = 800e-9
z0 = np.pi * w0**2 / lambd
                
#propagation des faisceaux gaussiens
def w(z):
    return w0*np.sqrt(1+(z/z0)**2)

def zr (la,om):
    return np.pi * om**2/la

r_list=np.linspace(-5*w0,5*w0,301)
zs = np.linspace(-2.9*zr(lambd,w0), 2.5*zr(lambd,w0), 301)
twin1.set_xlim(-19,14)
#twin1.legend(fontsize=textsize, loc='lower center')

twin2.plot(zs*1e6, 2*1*2.5e-3/(np.pi * 4e-15 * (w(zs)*1e2)**2), color='green')
twin2.set_ylabel('Laser peak intensity (W/cm²)', color='green',  fontsize=textsize)
twin2.tick_params(axis='y', labelcolor='green',labelsize=textsize)
twin2.get_yaxis().get_offset_text().set_position((1.3,0))


figz = plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'_no_ions_12ps____z.png', dpi=500, bbox_inches = 'tight')


















































#%%  figures supp; with 2 different CEP values, tirA and tirB
#num1 = 341   ; idx1 = np.where(fileNumbers == num1) ; idx1=int(idx1[0])
#num2 = 336   ; idx2 = np.where(fileNumbers == num2) ; idx2=int(idx2[0])
#CEP1 = +2.5
#CEP2 = -0.5
#
num1 = 262   ; idx1 = np.where(fileNumbers == num1) ; idx1=int(idx1[0])
num2 = 263  ; idx2 = np.where(fileNumbers == num2) ; idx2=int(idx2[0])
CEP1 = +3
CEP2 = 0




fig1 = plt.figure(figsize=(12,6))
gs = gridspec.GridSpec(2, 3, left=0.06, top=0.88, hspace=0.4, wspace = 0.3 , figure=fig1)    #y, x

fig1.add_subplot(gs[0, 0]) 
HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
H1 = np.loadtxt (HarmonicsFolder + 'tir' + str(num1) + '.dat' )
plt.pcolor(  H1[:,110:385] , cmap='inferno'  )
plt.colorbar()
plt.title('raw image, CEP=' + str(CEP1) + 'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
plt.ylim(80,570)


fig1.add_subplot(gs[0, 1]) 
H2 = np.loadtxt (HarmonicsFolder + 'tir' + str(num2) + '.dat' )
plt.pcolor(  H2[:,110:385] , cmap='inferno'  )
plt.ylim(80,570)
plt.colorbar()
plt.title('raw image, CEP=' + str(CEP2) + 'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')



fig1.add_subplot(gs[0, 2]) 
Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_tir601_VerticallyIntegrated.txt', delimiter=',', dtype='float')
Ev = Ev[0:np.size(Ev)-1]
rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
rescale = np.append(rescale, rescale[-1]) #1element longer
rescale=rescale[::-1]
H1int = (list (np.sum(H1[:,110:385],1)) * rescale)  #integrated and put in an array
H2int = (list (np.sum(H2[:,110:385],1)) * rescale)  #integrated and put in an array
plt.plot(Ev[::-1],H1int/max(H1int), color='purple', label = str(CEP1) + 'rad')
plt.plot(Ev[::-1],H2int/max(H2int), color='teal',  label = str(CEP2) + 'rad')
plt.xlim(11,35)
plt.legend()
plt.ylabel('normalized intensity (a.u.)')
plt.xlabel('photon energy (eV)')




fig1.add_subplot(gs[1, 0]) 
SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(num1) + '.tiff' ) 
SpaElectronsFile = np.array(SpaElectronsFile)
LANEXpx1 = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px  
#we need to get rid of hot pixels
temp = sorted(LANEXpx1.flatten())[::-1]    #array -> list, then organize from highest to lowest value
n_idx=[]
for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
    if temp[i]-temp[i+1] >= 100:
        n_idx=np.append(n_idx,i)
idx = np.where(LANEXpx1 >= temp[int(n_idx.max())])
for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
    LANEXpx1[idx[0][j], idx[1][j]] = LANEXpx1[idx[0][j]-1, idx[1][j]]       
plt.pcolor(LANEXpx1)
plt.title('raw image, CEP='+str(CEP1)+'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
plt.colorbar( )            # the mystery step ???????????
plt.xlabel('Azimuth angle $\phi$')
plt.ylabel('vertical axis')

                          
                          
fig1.add_subplot(gs[1, 1]) 
SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(num2) + '.tiff' ) 
SpaElectronsFile = np.array(SpaElectronsFile)
LANEXpx2 = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px 
#we need to get rid of hot pixels
temp = sorted(LANEXpx2.flatten())[::-1]    #array -> list, then organize from highest to lowest value
n_idx=[]
for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
    if temp[i]-temp[i+1] >= 100:
        n_idx=np.append(n_idx,i)
idx = np.where(LANEXpx2 >= temp[int(n_idx.max())])
for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
    LANEXpx2[idx[0][j], idx[1][j]] = LANEXpx2[idx[0][j]-1, idx[1][j]]      
plt.pcolor(LANEXpx2)  
plt.title('raw image, CEP='+str(CEP2)+'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
plt.xlabel('Azimuth angle $\phi$')
plt.colorbar()
plt.ylabel('vertical axis')



fig1.add_subplot(gs[1, 2]) 
plt.plot(np.rad2deg(temp_phi[22:-8]) ,ElecScan_r[22:-8,idx1]/max(ElecScan_r[22:-8,idx1]), color='purple', label = str(CEP1) + 'rad')
plt.plot(np.rad2deg(temp_phi[22:-8]) ,ElecScan_r[22:-8,idx2]/max(ElecScan_r[22:-8,idx2]), color='teal', label = str(CEP2) + '0rad')
plt.xlim(30,49)
#plt.ylim(0.65,1)
plt.legend()
plt.ylabel('electrons charge density (a.u.)')
plt.xlabel('Azimuth angle $\phi$ (deg)')






fig1 = plt.savefig(outdir+'\\CEPfigure_' + str(dateCode) + '_shots' + str(num1) + '_' + str(num2) +'.png', dpi=1000, bbox_inches = 'tight')






#



#%%
#%%  figures supp; with 4 different CEP values, tirA and tirB
num1 = 715   ; idx1 = np.where(fileNumbers == num1) ; idx1=int(idx1[0])
num2 = 737   ; idx2 = np.where(fileNumbers == num2) ; idx2=int(idx2[0])
num3 = 695   ; idx3 = np.where(fileNumbers == num3) ; idx3=int(idx3[0])
num4 = 669   ; idx4 = np.where(fileNumbers == num4) ; idx4=int(idx4[0])

CEP1 = -1
CEP2 = 0.5
CEP3 = 2
CEP4 = -2.5
           

#                               
#num1 = 262   ; idx1 = np.where(fileNumbers == num1) ; idx1=int(idx1[0])
#num2 = 263   ; idx2 = np.where(fileNumbers == num2) ; idx2=int(idx2[0])
#num3 = 262   ; idx3 = np.where(fileNumbers == num3) ; idx3=int(idx3[0])
#num4 = 263   ; idx4 = np.where(fileNumbers == num4) ; idx4=int(idx4[0])
#
#CEP1 = -1
#CEP2 = 0.5
#CEP3 = 2
#CEP4 = -2.5
#           


fig1 = plt.figure(figsize=(12,6))
gs = gridspec.GridSpec(2, 4, left=0.06, top=0.88, hspace=0.3, wspace = 0.1 , width_ratios=[1,1,1,1.2],figure=fig1)    #y, x
s=12
                      
                      
                      
                      
                      
HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
H1 = np.loadtxt (HarmonicsFolder + 'tir' + str(num1) + '.dat' )
H2 = np.loadtxt (HarmonicsFolder + 'tir' + str(num2) + '.dat' )
H3 = np.loadtxt (HarmonicsFolder + 'tir' + str(num3) + '.dat' )
H4 = np.loadtxt (HarmonicsFolder + 'tir' + str(num4) + '.dat' )

vmax =int( np.max([H3[:,110:385],H2[:,110:385],H4[:,110:385] ,H1[:,110:385]]))
vmin =int( np.min([H3[:,110:385],H2[:,110:385],H4[:,110:385] ,H1[:,110:385]]))

                      
                      
fig1.add_subplot(gs[0, 0]) 
plt.pcolor(  H1[:,110:385] , cmap='inferno', vmin=vmin, vmax=vmax  )
#plt.colorbar()
plt.title('raw image, CEP=' + str(CEP1) + 'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
plt.ylim(80,570)


fig1.add_subplot(gs[0, 1]) 
plt.pcolor(  H2[:,110:385] , cmap='inferno', vmin=vmin, vmax=vmax  )
plt.ylim(80,570)
#plt.colorbar()
plt.title('raw image, CEP=' + str(CEP2) + 'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')


fig1.add_subplot(gs[0, 2]) 
plt.pcolor(  H3[:,110:385] , cmap='inferno' , vmin=vmin, vmax=vmax )
plt.ylim(80,570)
#plt.colorbar()
plt.title('raw image, CEP=' + str(CEP3) + 'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')

fig1.add_subplot(gs[0, 3]) 
plt.pcolor(  H4[:,110:385] , cmap='inferno', vmin=vmin, vmax=vmax  )
plt.ylim(80,570)
plt.colorbar()
plt.title('raw image, CEP=' + str(CEP4) + 'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')







#
#fig1.add_subplot(gs[0, 2]) 
#Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_tir601_VerticallyIntegrated.txt', delimiter=',', dtype='float')
#Ev = Ev[0:np.size(Ev)-1]
#rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
#rescale = np.append(rescale, rescale[-1]) #1element longer
#rescale=rescale[::-1]
#H1int = (list (np.sum(H1[:,110:385],1)) * rescale)  #integrated and put in an array
#H2int = (list (np.sum(H2[:,110:385],1)) * rescale)  #integrated and put in an array
#plt.plot(Ev[::-1],H1int/max(H1int), color='purple', label = str(CEP1) + 'rad')
#plt.plot(Ev[::-1],H2int/max(H2int), color='teal',  label = str(CEP2) + 'rad')
#plt.xlim(11,35)
#plt.legend()
#plt.ylabel('normalized intensity (a.u.)')
#plt.xlabel('photon energy (eV)')



SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(num1) + '.tiff' ) 
SpaElectronsFile = np.array(SpaElectronsFile)
LANEXpx1 = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px  
#we need to get rid of hot pixels
temp = sorted(LANEXpx1.flatten())[::-1]    #array -> list, then organize from highest to lowest value
n_idx=[]
for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
    if temp[i]-temp[i+1] >= 100:
        n_idx=np.append(n_idx,i)
idx = np.where(LANEXpx1 >= temp[int(n_idx.max())])
for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
    LANEXpx1[idx[0][j], idx[1][j]] = LANEXpx1[idx[0][j]-1, idx[1][j]]  
    

SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(num2) + '.tiff' ) 
SpaElectronsFile = np.array(SpaElectronsFile)
LANEXpx2 = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px 
#we need to get rid of hot pixels
temp = sorted(LANEXpx2.flatten())[::-1]    #array -> list, then organize from highest to lowest value
n_idx=[]
for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
    if temp[i]-temp[i+1] >= 100:
        n_idx=np.append(n_idx,i)
idx = np.where(LANEXpx2 >= temp[int(n_idx.max())])
for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
    LANEXpx2[idx[0][j], idx[1][j]] = LANEXpx2[idx[0][j]-1, idx[1][j]]  
    
    
SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(num3) + '.tiff' ) 
SpaElectronsFile = np.array(SpaElectronsFile)
LANEXpx3 = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px 
#we need to get rid of hot pixels
temp = sorted(LANEXpx3.flatten())[::-1]    #array -> list, then organize from highest to lowest value
n_idx=[]
for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
    if temp[i]-temp[i+1] >= 100:
        n_idx=np.append(n_idx,i)
idx = np.where(LANEXpx3 >= temp[int(n_idx.max())])
for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
    LANEXpx3[idx[0][j], idx[1][j]] = LANEXpx3[idx[0][j]-1, idx[1][j]]  


SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(num4) + '.tiff' ) 
SpaElectronsFile = np.array(SpaElectronsFile)
LANEXpx4 = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px 
#we need to get rid of hot pixels
temp = sorted(LANEXpx4.flatten())[::-1]    #array -> list, then organize from highest to lowest value
n_idx=[]
for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
    if temp[i]-temp[i+1] >= 100:
        n_idx=np.append(n_idx,i)
idx = np.where(LANEXpx4 >= temp[int(n_idx.max())])
for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
    LANEXpx4[idx[0][j], idx[1][j]] = LANEXpx4[idx[0][j]-1, idx[1][j]] 
    
    
    
    
    
    
    
    
vmaxe = int(np.max( [LANEXpx1, LANEXpx2, LANEXpx3, LANEXpx4]  ))
vmine = int(np.min( [LANEXpx1, LANEXpx2, LANEXpx3, LANEXpx4]  ))    
    
fig1.add_subplot(gs[1, 0])      
plt.pcolor(LANEXpx1, vmin=vmine, vmax=vmaxe)
plt.title('raw image, CEP='+str(CEP1)+'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')            
plt.xlabel('Azimuth angle $\phi$', fontsize=s)
plt.ylabel('vertical axis', fontsize=s)
                         
                          
fig1.add_subplot(gs[1, 1])     
plt.pcolor(LANEXpx2, vmin=vmine, vmax=vmaxe)  
plt.title('raw image, CEP='+str(CEP2)+'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')


fig1.add_subplot(gs[1, 2])     
plt.pcolor(LANEXpx3, vmin=vmine, vmax=vmaxe)  
plt.title('raw image, CEP='+str(CEP3)+'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')


fig1.add_subplot(gs[1, 3])      
plt.pcolor(LANEXpx4, vmin=vmine, vmax=vmaxe)  
plt.title('raw image, CEP='+str(CEP4)+'rad')
plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
plt.colorbar()



#
#
#fig1.add_subplot(gs[1, 2]) 
#plt.plot(np.rad2deg(temp_phi[22:-8]) ,ElecScan_r[22:-8,idx1]/max(ElecScan_r[22:-8,idx1]), color='purple', label = str(CEP1) + 'rad')
#plt.plot(np.rad2deg(temp_phi[22:-8]) ,ElecScan_r[22:-8,idx2]/max(ElecScan_r[22:-8,idx2]), color='teal', label = str(CEP2) + '0rad')
#plt.xlim(30,49)
##plt.ylim(0.65,1)
#plt.legend()
#plt.ylabel('electrons charge density (a.u.)')
#plt.xlabel('Azimuth angle $\phi$ (deg)')
#





fig1 = plt.savefig(outdir+'\\CEPfigure_' + str(dateCode) + '_shots' + str(num1) + '_' + str(num2) +'_'+str(num3)+'_'+str(num4)+'.png', dpi=1000, bbox_inches = 'tight')

























#%% figure article CEP, espectro + harmonics pour 2 valeurs
#fig = plt.figure(figsize=(10,10))
#textsize = 16
#fig.text (0.2,0.98, 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=textsize -4)
#
#
#gs = gridspec.GridSpec(2, 2, left=0.06, top=0.88, height_ratios=[1,1], figure=fig)
#
#fig.add_subplot(gs[0, 0])
#plt.plot(Ev[::-1], Bharm[:,0]/np.max(Bharm))  #harm, CEP=-3, modulated
#plt.xlim(10,36)
#plt.ylim(0,1)
#plt.xlabel('Photon energy (eV)')
#
#
#fig.add_subplot(gs[0, 1])
#plt.plot(E_interp,  Belec[:,0]/np.max(Belec))
#
#
#
#
#
#fig.add_subplot(gs[1, 0])
#plt.plot(Ev[::-1], Bharm[:,6]/np.max(Bharm))   #harm, CEP=0, continuum
#plt.xlim(10,36)
#plt.ylim(0,1)
#plt.xlabel('Photon energy (eV)')
#
#
#
#fig.add_subplot(gs[1, 1])
#plt.plot(E_interp,  Belec[:,6]/np.max(Belec))
#plt.xlim(0.3,2)
#







################################################### same but together
#fig = plt.figure(figsize=(6,6))
#textsize = 16
#fig.text (0.2,0.95, 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=textsize -4)
#
#n1=11  #mod
#n2=5   #cont
#
#gs = gridspec.GridSpec(2, 1, left=0.06, top=0.88, height_ratios=[1,1], hspace=0.3, figure=fig)
#
#fig.add_subplot(gs[0, 0])
#plt.plot(Ev[::-1], Bharm[:,n1]/np.max(Bharm), label=str(param_uniq[n1])+'rad')
#plt.plot(Ev[::-1], Bharm[:,n2]/np.max(Bharm), label=str(param_uniq[n2])+'rad')  
#plt.xlim(10,36)
#plt.ylim(0,1)
#plt.xlabel('Photon energy (eV)')
#plt.legend()
#
#
#fig.add_subplot(gs[1, 0])
#plt.plot(E_interp,  Belec[:,n1]/np.max(Belec), label=str(param_uniq[n1])+'rad')
#plt.plot(E_interp,  Belec[:,n2]/np.max(Belec), label=str(param_uniq[n2])+'rad')
#plt.xlim(E_interp[-1],2)
#plt.ylim(ymin=0)
#plt.xlabel('Electron energy (MeV)')
#plt.legend()
#
#
#plt.savefig(r'C:\Users\ouille\Desktop\article_CEP_effects' + '\\' + str(dateCode)+ '_' + str(file1) + '_' + str(file_last) + '_cep-' + str(n1) +'_'+str(n2)+'_espectro' + '.png')











#%% 20210429   Figure CLEO spectres d'électrons  VS CEP
#
#fig_elec = plt.figure(figsize=(np.size(param_uniq)*1.5,3.5))    #5 for big
#gs = gridspec.GridSpec(1,np.size(param_uniq), left=0.06, top=0.88, wspace=0.04, hspace=0.02, figure=fig_elec)
#
#for s, value in enumerate(param_uniq):
#    fig_elec.add_subplot(gs[ 0,s])  
#    plt.plot(Belec[:,s]/Belec.max(), E_interp, color='Navy');plt.ylim(E_interp[-1],1.5);plt.xlim(1, -0.02)
#    if s!=0:
#        plt.tick_params( top='off', labeltop='off', left='off', labelleft='off', labelbottom='off')
#    else:
#        plt.ylabel('Electron energy (MeV)', fontsize=textsize-3)
#        plt.xlabel('Normalized intensity', fontsize=textsize-3)
#        plt.tick_params(labelsize=textsize-3)
#    
#
#plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python' + '\\' + str(dateCode)+ '_' + str(file1) + '_' + str(file_last) + '_cep_espectro_CLEO' + '.png', bbox_inches='tight')
