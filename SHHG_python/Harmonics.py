# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 15:19:06 2020

@author: ouille
"""

"""
Created on Wed Jul 31 17:13:37 2019

@author: ouille
This is a (almost) a copy of the matlab scripts for SHHG analysis

THINGS TO DO : 
- check that it gets organized nicely for random scans
- add protons when there is not specturm.txt file
- comment on fait quand il y a deux tirs au meme delai? moyenner!
- add electrons spectrometer 
- ajouter plasma speed/gradient length
- make it faster for the electrons analysis
- créer dossier date pour sauvegarde des données s'il n'existe pas deja dans outdir
"""



#%%   #import libraries
import numpy as np
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




#%%  ##indicate file paths and load data.... things to change :
remoteDrive = r'Z:'        #file path    Z: or C: drive 
#dataPath = 'SHHG\Jaismeen_Marie'      #dataPath = e.g : 'SHHG\Stefan_and_Marie\\2nd-campaign-May2019'  
dataPath = 'SHHG\Stefan_and_Marie\\2nd-campaign-May2019'  
outdir = r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python'

save = 1    # save=0 : don't save, save = 1 : save figures
harm = 1     ;  elec = 0    ; ions = 0  ; espectro = 0

##data to work on : #image exemple
#dateCode ='20190605' #exemple: '20200629',  chirped: '20190605' #'20190829'
#file = 531 #531 #66   #exemple: 98

#scan iris
#dateCode = '20200619'
#file=470
#    y_low = 110
#    y_high = 385
#
##scan, GDD+500fs2
dateCode = '20190605'
#file=189  #chirped
file=31  #unchurped
y_low = 39
y_high = 359

  
#example image
#dateCode='20200629'
#file=66
#y_low = 108
#y_high = 388



##load calibration file for harmonics :
Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_VIntegratedharmonics_tir31.txt', delimiter=',', dtype='float')   #for the iris scan
#Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_tir38_VerticallyIntegrated.txt', delimiter=',', dtype='float')   #for the GDD scan





#load DataLog
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


ScanDelays = -PrepulseDelay[file-offset_shot]
scan_param = 'Gradient'
if scan_param == 'CEP':
    CEP = np.loadtxt(DataLog, dtype= str, skiprows = (2), usecols=(9,) ) 
    CEP = CEP[file - offset_shot]
    CEP = [float(CEP[i]) for i in range(len(CEP))]




 
                   
#%%                   
if harm == 1 : 

#    
   #load calibration for harmonics energies :
    Ev = Ev[0:np.size(Ev)]
    rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
    rescale = np.append(rescale, rescale[-1]) #1element longer
    #rescale=rescale[::-1]
    
    HarmonicsScan = []
    k = file
    HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
    HarmonicsFile = np.loadtxt (HarmonicsFolder + 'tir' + str(k) + '.dat' ) 
    
    
#%%
#calib divergecne? iris scan
    plt.figure()
    plt.imshow(HarmonicsFile, cmap='inferno', vmin=0, vmax=1000)
    
    #iris scan
#    plt.axvline(40, color='orange')
#    plt.axvline(424,color='orange')
#    plt.axvline(110, color='orange')
#    plt.axvline(382,color='orange')
    
    #GDD scan
    plt.axvline(39, color='orange')
    plt.axvline(359,color='orange')
    plt.axvline(0, color='orange')
    plt.axvline(415,color='orange')
    
    
#    #example scan
#    plt.axvline(36, color='red')
#    plt.axvline(430,color='red')
#    plt.axvline(108, color='orange')
#    plt.axvline(388,color='orange')
    
#%%


    
    
    
    HarmonicsScan.append(list (np.sum(HarmonicsFile[:,y_low:y_high],1)) * rescale)  #integrated and put in an array
    #calibrate Y axis : énergies en eV à gauche et à droite wavelength ou Harmonics number au choix
#    
#    
    #flip 90degree to get the delays on the horizontal axis, more convenient
    HarmFile_r = np.zeros ((np.shape(HarmonicsFile)[1],np.shape(HarmonicsFile)[0]))
    for i, osef in enumerate (np.zeros(np.shape(HarmonicsFile)[0])):
        for j, osef_aussi in enumerate (np.zeros(np.shape(HarmonicsFile)[1])):
            HarmFile_r[j,i] = HarmonicsFile[i][j]
    
    
    #line for CWE cutoff at max palsma frequency for a given wavelenght
    wp = np.sqrt (ct.e**2*ct.ne / (ct.me*ct.e0)) 
    lp = 2*np.pi*ct.c/wp
    Ev_max = ct.h*ct.c/lp/ct.e                
    
    
    
    
    
    
    Ev=Ev[::-1]
    Ev_round = [round(Ev[i],1) for i in range(len(Ev))]
    idx_wmax = np.where(Ev_round==round(Ev_max,1))
    
    
    
    
    xlist = np.arange(0, np.shape(HarmFile_r[y_low:y_high,70:570])[0])
    ### divergence : ((xlist-138)/(138/43))
    
    
    
    
    
    s=18 #16
    
    
 #   fig = plt.figure(figsize=(6.5,5))  #example
        
    fig = plt.figure(figsize=(5,5))
#    fig.suptitle('Image brute', fontsize=s+2)
    gs = gridspec.GridSpec(2, 2, left=0.06, top=0.88, height_ratios=[1,3], width_ratios=[7.3,1], hspace=0, figure=fig)
    
    fig.add_subplot(gs[1, :])
    #plt.pcolor(HarmFile_r[y_low-20:y_high+20,70:570]/  (HarmFile_r[y_low:y_high,70:570]).max(), cmap='inferno')
    plt.pcolor(HarmFile_r[y_low-20:y_high+20,70:570], cmap='inferno')

    plt.ylabel('divergence (mrad)', fontsize=s)
    #plt.yticks([])
    plt.xlabel('$E_{photon}$ (eV)', fontsize=s)
   # plt.yticks([20,69+20, 138+20,138+69+20,2*138+20] , [-32,-16,0,16,32], fontsize=s)   #EXAMPLE
    plt.yticks([20,110,160+20,270,2*160+20] , ['-35','','0','','35'], fontsize=s)  #GDD
   # plt.yticks([10,138,266] , [-40,0,40], fontsize=s)
    plt.xticks( (0,100,200,300,400,500)  , (round(Ev[0+70],1),  round(Ev[100+70],1),  round(Ev[200+70],1) ,  round(Ev[300+70],1) , round(Ev[400+70],1) ,  round(Ev[500+70],1)  ), fontsize=s  )
    plt.axvline(x=idx_wmax[0][0]-70, linestyle= '--', color='orange', label='$\omega_p$')
#    plt.legend(loc='upper right', fontsize=s)
    #plt.colorbar()
    #plt.xlim(xmin=25)
    #test colorbar
    cbaxese = fig.add_axes([0.97, 0.125, 0.03, 0.48]) #left, bottom, width, height 
    cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
    cbaxese.formatter.set_powerlimits((0, 0))
    cbaxese.ax.yaxis.get_offset_text().set_fontsize(s-4)
    cbaxese = plt.tick_params(labelsize=s-4)
    cbaxese = plt.ylabel('counts (a.u)',  fontsize=s-2 )
        
    
    
#    plt2 = plt.twiny()
    list_Ev = (round(Ev[0+70],1),  round(Ev[100+70],1),  round(Ev[200+70],1) ,  round(Ev[300+70],1) , round(Ev[400+70],1) ,  round(Ev[500+70],1)  )
    h = 6.62607004 * 1e-34 #m2 kg / s
    c = 299792458 #m/s 
    e = 1.602176565*1e-19 #Joules dans 1eV
    list_l = [ int(round(h*c/(list_Ev[k]*e)*1e9,0) ) for k in range(len(list_Ev))]
#    plt.xticks((0,100,200,300,400,500)  ,  list_l*2 , fontsize=s )
#    plt.xlabel('$\lambda$ (nm)', fontsize=s)
    plt.tick_params(labelsize=s)
    
    
    #integrated
    ax=fig.add_subplot(gs[0,:])
    ax.yaxis.tick_right()
    plt.plot(sum(HarmFile_r[y_low-10:y_high+10,70:570])/np.max(sum(HarmFile_r[y_low-10:y_high+10,70:570])), color='purple')  
    plt.axvline(x=idx_wmax[0][0]-70, linestyle= '--', color='orange', label='$\omega_p$')
    plt.legend(loc='upper right', fontsize=s)
    plt.yticks([0,0.5,1], fontsize=s)
    plt.ylabel('I (a.u)', fontsize=s)
    plt.xlim(0, len( sum(HarmFile_r[y_low-10:y_high+10,70:570])))
    plt.xticks([])
    plt2 = plt.twiny()
    plt.xticks((0,100,200,300,400,500)  ,  [x for x in list_l], fontsize=s )
    plt.xlabel('$\lambda$ (nm)', fontsize=s)


    
  #  fig = plt.savefig(r'F:\Abtracts_confs_posters\mini-soutenance\images_mini_soutenance\\raw_harm.png',bbox_inches='tight', dpi=800)
    fig = plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python\\raw_harm_tir' +str(file) +'_' + str(dateCode) + '.png',bbox_inches='tight', dpi=600)

#
##%%tes div mesure... mauvaise idéee
#
#plt.figure()
#plt.pcolor(Ev[70:570],np.arange(-137,138,1)*0.1953/834*1000,HarmFile_r[y_low:y_high,70:570])
#plt.axvline(13.4)
##
#px_posh=int(np.mean(np.where(np.round(Ev,1)==13.4)))
#print(px_posh)
#
#px_posr=int(np.mean(np.where(np.round(Ev,1)==35)))
#print(px_posr)
#
#
#plt.figure()
#plt.plot(np.arange(-137,138,1)*0.1953/834*1000, HarmFile_r[y_low:y_high,px_posh] - HarmFile_r[y_low:y_high,px_posr] )
#
#maxi = np.max(HarmFile_r[y_low:y_high,px_posh] - HarmFile_r[y_low:y_high,px_posr] )
#x_m = np.where( HarmFile_r[y_low:y_high,px_posh] - HarmFile_r[y_low:y_high,px_posr]== maxi)
#             
#             
#plt.axhline(maxi/2)
#plt.axvline((-137+x_m[0][0])*0.1953/834*1000)
#plt.axvline(12)
#
##    fig2 = plt.figure(figsize=(8,4))
##    fig2.suptitle('image brute - zoom à droite')
##    plt.pcolor(HarmFile_r[y_low-10:y_high+10,400:570]/ (HarmFile_r[y_low-10:y_high+10,70:570]).max(), cmap='inferno')
##    plt.ylabel('vertical axis')
##    plt.yticks([])
##    plt.xticks((0,20,40,60,80,100,120,140,160),(   round(Ev[0+400],1),  round(Ev[20+400],1),  round(Ev[40+400],1) ,  round(Ev[60+400],1) , round(Ev[80+400],1) ,  round(Ev[100+400],1), round(Ev[120+400],1)  , round(Ev[140+400],1) ,   round(Ev[160+400],1)        )    )
##    plt.xlabel('photon energy (eV)')
##    plt.axvline(x=idx_wmax[0][0]-400, linestyle= '--', color='white', label='$\omega_p$')
##    plt.legend(loc='upper right')
##    plt.colorbar()
###    