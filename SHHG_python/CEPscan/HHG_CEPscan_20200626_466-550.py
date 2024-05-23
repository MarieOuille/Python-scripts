# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 17:13:37 2019

@author: ouille


For now, this script is to plot HH integrated spectra VS prepulse delays

goal = 
add electrons saptial and spectral, protons and so on, optional



THINGS TO DO : 
check that it gets organized nicely for random scans
"""



#%%   #import libraries
import numpy as np
from PIL import Image
import sys
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
if not r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' in sys.path : 
    sys.path.insert(0,r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' )
from phys_constants import ct
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['#BF360C','#F4511E', '#f3c1bb', 'white' ,'#d9e1f9' , '#a0b4f0','#6d8ce8','royalblue', '#2250d9' ,'#0000e6','darkblue', 'midnightblue','black' ])



#%%  ##indicate file paths and load data.... things to change :
remoteDrive = r'Z:'
dataPath = 'SHHG\Jaismeen_Marie' 
figSavePath = r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files'


#data to work on :
dateCode = '20200626'#'20190605' #'20190829'
file1 = 466 #531 #66 
file_last = 550 #552 #67 


##load calibration file:
Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_tir601_VerticallyIntegrated.txt', delimiter=',', dtype='float')




#%%
fileNumbers = np.arange(file1, file_last + 1,1)



#load DataLog
DataLog = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\\DataLog_with_CEP.txt' 
                                  
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
ScanDelays = -PrepulseDelay[fileNumbers-offset_shot]
CEP = np.loadtxt(DataLog, dtype= str, skiprows = (2), usecols=(10,) ) 
CEPs = CEP[fileNumbers-offset_shot]
CEPs = CEPs.astype(np.float)            
                   
                   
#%%                   
#load calibration for harmonics energies :
Ev = Ev[0:np.size(Ev)-1]
rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
rescale = np.append(rescale, rescale[-1]) #1element longer


HarmonicsScan = []
for j, k in enumerate (fileNumbers):
    HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
    HarmonicsFile = np.loadtxt (HarmonicsFolder + 'tir' + str(k) + '.dat' ) 
    HarmonicsScan.append(list (np.sum(HarmonicsFile,1)) * rescale)  #integrated and put in an array
#calibrate Y axis : énergies en eV à gauche et à droite wavelength ou Harmonics number au choix




#%% PLOT HARMONICS scans


##original plot, with delays on the y axis
##plt.pcolor(np.arange(0,np.shape(HarmonicsScan)[1],1), ScanDelays, HarmonicsScan, cmap = 'jet')
#plt.pcolor(Ev[::-1], ScanDelays, HarmonicsScan, cmap = custom_cmap)
#plt.title ('Harmonics Gradient Scan')
#plt.ylabel('Prepulse delay (ps)')
#plt.xlabel('Photon energy (eV)')
#plt.xlim(12,37)









#flip 90degree to get the delays on the horizontal axis, more convenient
HarmScan_r = np.zeros ((np.shape(HarmonicsScan)[1],np.shape(HarmonicsScan)[0]))
for k, osef in enumerate (np.zeros(np.shape(HarmonicsScan)[0])):
    for j, osef_aussi in enumerate (np.zeros(np.shape(HarmonicsScan)[1])):
        HarmScan_r[j,k] = HarmonicsScan[k][j]
        
        
        
    
B = np.zeros((640,13))
ar_CEP = np.arange(-3,3.5,0.5)
for j, valuej in enumerate (ar_CEP):
    A = np.zeros((640,5))
    count = 0
    for i, value in enumerate (fileNumbers):
        if CEPs [i] == valuej :
            A[:,count] = HarmScan_r[:,i]
            if count < 4 :
                count = count +1 
    B[:,j]=(A[:,0]+A[:,1]+A[:,2]+A[:,3]+A[:,4])/5

        
     
     
        

#plot the figure
fig_h = plt.figure(figsize=(9,6))
fig_h.suptitle ( 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=10 )
plt.pcolor(np.arange(-3,4,0.5), Ev[::-1], B/np.max(B), cmap = 'inferno')
#plt.title ('Harmonics Gradient Scan', fontsize=11)
plt.xlabel('relative CEP (rad)')
plt.ylabel('Photon energy (eV)')
plt.ylim(ymax=39)

#line for CWE cutoff at max palsma frequency for a given wavelenght
wp = np.sqrt (ct.e**2*ct.ne / (ct.me*ct.e0)) 
lp = 2*np.pi*ct.c/wp
Ev_max = ct.h*ct.c/lp/ct.e                
plt.axhline(y=Ev_max, color='red', label='$E_{p,max}$')
plt.legend()


#colorbar
cbaxes = fig_h.add_axes([0.915, 0.125, 0.03, 0.76]) #left, bottom, width, height
cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
cbaxes = plt.ylabel('Spectral Intensity (a.u)')


#savefigure
fig_h = plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\Harmonics_CEPscan_from_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last)  +'.png'  ,dpi=400,  bbox_inches='tight'    )                
                          
                          






plt.figure(figsize=(6,5))
plt.plot(Ev[::-1], B[:,0]/np.max(B), label='CEP =' + str(ar_CEP[0]) + 'rad' )
plt.plot(Ev[::-1], B[:,6]/np.max(B),  label='CEP =' + str(ar_CEP[6]) + 'rad')
plt.legend() 
plt.xlabel('photon energy (eV)')
plt.ylabel('spectral intensity (a.u)')
plt.xlim(10,39)
plt.ylim(0,1)
plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'CEP'+str(ar_CEP[0])+'_and_'+str(ar_CEP[6]) +'.png'  ,dpi=400,  bbox_inches='tight'    )                



plt.figure(figsize=(6,5))
plt.plot(Ev[::-1], B[:,1]/np.max(B), label='CEP =' + str(ar_CEP[1]) + 'rad' )
plt.plot(Ev[::-1], B[:,7]/np.max(B),  label='CEP =' + str(ar_CEP[7]) + 'rad')
plt.legend() 
plt.xlabel('photon energy (eV)')
plt.ylabel('spectral intensity (a.u)')
plt.xlim(10,39)
plt.ylim(0,1)
plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'CEP'+str(ar_CEP[1])+'_and_'+str(ar_CEP[7]) +'.png'  ,dpi=400,  bbox_inches='tight'    )                



plt.figure(figsize=(6,5))
plt.plot(Ev[::-1], B[:,2]/np.max(B), label='CEP =' + str(ar_CEP[2]) + 'rad' )
plt.plot(Ev[::-1], B[:,8]/np.max(B),  label='CEP =' + str(ar_CEP[8]) + 'rad')
plt.legend() 
plt.xlabel('photon energy (eV)')
plt.ylabel('spectral intensity (a.u)')
plt.xlim(10,39)
plt.ylim(0,1)
plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'CEP'+str(ar_CEP[2])+'_and_'+str(ar_CEP[8]) +'.png'  ,dpi=400,  bbox_inches='tight'    )                



plt.figure(figsize=(6,5))
plt.plot(Ev[::-1], B[:,3]/np.max(B), label='CEP =' + str(ar_CEP[3]) + 'rad' )
plt.plot(Ev[::-1], B[:,9]/np.max(B),  label='CEP =' + str(ar_CEP[9]) + 'rad')
plt.legend() 
plt.xlabel('photon energy (eV)')
plt.ylabel('spectral intensity (a.u)')
plt.xlim(10,39)
plt.ylim(0,1)
plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'CEP'+str(ar_CEP[3])+'_and_'+str(ar_CEP[9]) +'.png'  ,dpi=400,  bbox_inches='tight'    )                



plt.figure(figsize=(6,5))
plt.plot(Ev[::-1], B[:,4]/np.max(B), label='CEP =' + str(ar_CEP[4]) + 'rad' )
plt.plot(Ev[::-1], B[:,10]/np.max(B),  label='CEP =' + str(ar_CEP[10]) + 'rad')
plt.legend() 
plt.xlabel('photon energy (eV)')
plt.ylabel('spectral intensity (a.u)')
plt.xlim(10,39)
plt.ylim(0,1)
plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'CEP'+str(ar_CEP[4])+'_and_'+str(ar_CEP[10]) +'.png'  ,dpi=400,  bbox_inches='tight'    )                



plt.figure(figsize=(6,5))
plt.plot(Ev[::-1], B[:,5]/np.max(B), label='CEP =' + str(ar_CEP[5]) + 'rad' )
plt.plot(Ev[::-1], B[:,11]/np.max(B),  label='CEP =' + str(ar_CEP[11]) + 'rad')
plt.legend() 
plt.xlabel('photon energy (eV)')
plt.ylabel('spectral intensity (a.u)')
plt.xlim(10,39)
plt.ylim(0,1)
plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'CEP'+str(ar_CEP[5])+'_and_'+str(ar_CEP[11]) +'.png'  ,dpi=400,  bbox_inches='tight'    )                



plt.figure(figsize=(6,5))
plt.plot(Ev[::-1], B[:,6]/np.max(B), label='CEP =' + str(ar_CEP[6]) + 'rad' )
plt.plot(Ev[::-1], B[:,12]/np.max(B),  label='CEP =' + str(ar_CEP[12]) + 'rad')
plt.legend() 
plt.xlabel('photon energy (eV)')
plt.ylabel('spectral intensity (a.u)')
plt.xlim(10,39)
plt.ylim(0,1)
plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'CEP'+str(ar_CEP[6])+'_and_'+str(ar_CEP[12]) +'.png'  ,dpi=400,  bbox_inches='tight'    )                



##single plot
#plt.figure()
#plt.plot(Ev[::-1], HarmonicsScan[10] ) 
#plt.xlim(10,37)
#plt.xlabel('Photon energy (eV)')
#plt.axvline(x=30.17, color='red', label='$w_{p,max}$')
#plt.legend()



#plt.figure(figsize=(6,5))
#plt.plot(Ev[::-1], B[:,1]/np.max(B), label='CEP =' + str(ar_CEP[1]) + 'rad' )
#plt.plot(Ev[::-1], B[:,4]/np.max(B),  label='CEP =' + str(ar_CEP[4]) + 'rad')
#plt.plot(Ev[::-1], B[:,10]/np.max(B),  label='CEP =' + str(ar_CEP[10]) + 'rad')
#plt.plot(Ev[::-1], B[:,7]/np.max(B),  label='CEP =' + str(ar_CEP[7]) + 'rad')
#plt.legend() 
#plt.xlabel('photon energy (eV)')
#plt.ylabel('spectral intensity (a.u)')
#plt.xlim(10,39)
#plt.ylim(0,1)
#plt.savefig(  r'Z:\SHHG\Analyzed Data\20200626\\int_harm_' + str(dateCode) + '_shots' + str(file1) + '_to_' + str(file_last) + 'multiCEP' +'.png'  ,dpi=400,  bbox_inches='tight'    )                
#
#






#%% NEXT : electrons ?
##for j, k in enumerate (fileNumbers):
#k = 55
#SpaElectronsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\electrons\\' 
#SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(k) + '.tiff' ) 
#SpaElectronsFile = np.array(SpaElectronsFile)
#fig_es = plt.figure(figsize=(6,5.6))
#plt.pcolor(SpaElectronsFile[200:800,450:1010], cmap = 'viridis')
#plt.xlabel('Azimuth angle $\phi$ (deg)')
#plt.ylabel('Polar angle' +  r' ${\Theta}$' + ' (deg)')
#
#
##colorbar
#cbaxes = fig_es.add_axes([0.915, 0.125, 0.03, 0.76]) #left, bottom, width, height
#cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
#cbaxes = plt.ylabel('Radiant Charge Density (nC/sr)')
#
#
