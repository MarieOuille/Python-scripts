# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 10:47:18 2020

@author: ouille

this is also for the big lanex plots
"""


#%%   #import libraries
import numpy as np
import re
from PIL import Image
import sys
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
if not r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' in sys.path : 
    sys.path.insert(0,r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' )
from phys_constants import ct
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['#BF360C','#F4511E', '#f3c1bb', 'white' ,'#d9e1f9' , '#a0b4f0','#6d8ce8','royalblue', '#2250d9' ,'#0000e6','darkblue', 'midnightblue','black' ])
import scipy as sc 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from matplotlib import cm



#%%  ##indicate file paths and load data.... things to change :
remoteDrive = r'Z:'        #file path    Z: or C: drive 
dataPath = 'SHHG\Jaismeen_Marie'      #dataPath = e.g : 'SHHG\Stefan_and_Marie\\2nd-campaign-May2019'  
outdir = r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python'

save = 1        # save=0 : don't save, save = 1 : save figures
harm = 1        #1 to look at harmonics, otherwise 0
elec = 1        #1 to look at spatial electrons, otherwise 0

#data to work on :
    
    ##  BIG LANEX CEP scans
#dateCode='20201005' ; fileNumbers=np.arange(75,102)
#dateCode='20201005' ; fileNumbers=np.arange(102,129)
#dateCode='20201005' ; fileNumbers=np.arange(132,158)


#dateCode='20201001' ; fileNumbers=np.arange(205,232)    #2ps
#dateCode='20201001' ; fileNumbers=np.arange(261,290)   #3ps
#dateCode='20201001' ; fileNumbers=np.arange(291,319)   #-1ps
#dateCode='20201001' ; fileNumbers=np.arange(321,347)   #4ps

   # LANEX full screen scans
#dateCode='20200923' ; fileNumbers=np.arange(154,180+1)
#dateCode = '20200824' ; fileNumbers = np.arange(364 , 389+1)   #2ps
#dateCode = '20200824' ; fileNumbers = np.arange(326 , 351+1)   #3ps
#dateCode = '20200824' ; fileNumbers = np.append ( np.arange(457,462+1), np.arange(465,484+1) )    #4ps
#dateCode = '20200824' ; fileNumbers = np.arange(411 , 438+1)   #1ps                                            

dateCode = '20200902' ; fileNumbers = np.arange(325, 350+1)   #2ps   confs shots  THE FIGURE!!!!

                                                
                                                
                                                
                                                
#dateCode = '20200902' ; fileNumbers = np.arange(356, 382+1)   #3ps                               
                                            
#dateCode = '20200923' ; fileNumbers = np.arange(288, 313+1)   #2ps                                             
                                                                   
#dateCode = '20200904' ; fileNumbers = np.arange(289, 315+1)   #2ps                                                 
                                               
#dateCode='20200904' ; fileNumbers=np.array([297,301])
#dateCode='20200824'  ; fileNumbers=np.array([329, 341, 351, 336])
#dateCode='20200626' ; fileNumbers = np.array([715,737,695,669])
#dateCode = '20200824' ; fileNumbers = np.array([507,523,513,520])
#dateCode = '20200902' ; fileNumbers = np.array([332,346,334,325])


#dateCode = '20200902' ; fileNumbers = np.array([346,325])
#dateCode = '20200824' ; fileNumbers = np.array([341, 336])
#dateCode = '20200626' ; fileNumbers = np.array([737,669])

#dateCode='20200629'; fileNumbers = np.array([262,263])

#dateCode = '20200626' ; fileNumbers = np.array([715,737,695,669])

#dateCode = '20200824' ; fileNumbers = np.array ([460, 465])



##load calibration file for harmonics :
Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_VIntegratedharmonics_tir31.txt', delimiter=',', dtype='float')















#%%   #load DataLog

DataLog = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\\DataLog.txt' 
                                  
ShotNo = np.loadtxt(DataLog, skiprows = (2), usecols=(0,) )
if ShotNo[0] >0:   #offset in shot num.....
    offset_shot=int(ShotNo[0])
#PulseEnergy = np.loadtxt(DataLog, skiprows = (2), usecols=(1,) )  #mJ
#PulseDuration = np.loadtxt(DataLog, skiprows = (2), usecols=(2,) )  #fs
#Chirp = np.loadtxt(DataLog, skiprows = (2), usecols=(3,) )  #fs²
#WedgeInsert = np.loadtxt(DataLog, skiprows = (2), usecols=(4,) )  #mm
#ParabolaPos = np.loadtxt(DataLog, skiprows = (2), usecols=(5,) )  #µm
#PrepulseDelay = np.loadtxt(DataLog, skiprows = (2), usecols=(6,) ) #ps
#SeqLen = np.loadtxt(DataLog, skiprows = (2), usecols=(7,) )   #ms
#CEP = np.loadtxt(DataLog, dtype= str, skiprows = (2), usecols=(9,) ) 
#coms=np.loadtxt(DataLog, dtype= str, comments='CEP=',delimiter='\t',skiprows = (2), usecols=(11,) )
coms=np.loadtxt(DataLog, dtype= str,delimiter='\t',skiprows = (2), usecols=(11,) )




#ScanDelays = -PrepulseDelay[fileNumbers-offset_shot]
#scan_param = 'CEP' #'CEP'
#CEPs = CEP[fileNumbers-offset_shot]
#CEPs = [float(CEPs[i]) for i in range(len(CEPs))]
#


CEPs=np.zeros(np.size(ShotNo))
for r, t in enumerate (np.arange(0,np.size(ShotNo)))  :
   # a = re.findall("[+-]\d+", coms[r])
    
   
   #a=re.findall("[+-]?\d+\.\d+", coms[r])
    a=re.findall("[+-]?\d\.\d", coms[r]) 
    #a=re.findall("[+-]?\d\.\d", coms[r][58:-1]) 

    if np.shape(a)[0] == 0 :
        #a_temp = (re.findall("[+-]?\d", coms[r]))
        a_temp = re.findall("[+-]?\d", coms[r])
        #a_temp = re.findall("[+-]?\d", coms[r][58:-1])
        if a_temp :
            a.append(a_temp[0])
            
    if np.shape(a)[0] == 0 :
        CEPs[r] = 0
    else :
        CEPs[r] = float(a[0])
        

cep_uniq = np.unique(CEPs[fileNumbers-offset_shot])









 
                   
#%%                   
if harm == 1 : 
   #load calibration for harmonics energies :
    Ev = Ev[0:np.size(Ev)]
    rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
    rescale = np.append(rescale, rescale[-1]) #1element longer
    rescale=rescale[::-1]
    if dateCode == '20200626' :
        y_low = 110 ; y_high = 385      ##checked
    if dateCode == '20200629' :
        y_low = 110 ; y_high = 385      
    elif dateCode >= '20200902' and dateCode <='20200904' :
        y_low =99 ; y_high = 376
    elif dateCode=='20200824':
        y_low =98 ; y_high = 376
    Harmpics=[]
    HarmonicsScan = []
    for j, k in enumerate (fileNumbers):
        HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
        HarmonicsFile = np.loadtxt (HarmonicsFolder + 'tir' + str(k) + '.dat' ) 
        #Harmpics.append(HarmonicsFile[:,95:377])
        #Harmpics.append(HarmonicsFile[:,100:375])
        Harmpics.append(HarmonicsFile[:,y_low:y_high])        
        #Harmpics.append(HarmonicsFile[:,108:386])
        #HarmonicsScan.append(list (np.sum(HarmonicsFile[:,95:377],1)) * rescale)  #integrated and put in an array
        HarmonicsScan.append(list (np.sum(HarmonicsFile[:,y_low:y_high],1)) * rescale)  #integrated and put in an array
        #HarmonicsScan.append(list (np.sum(HarmonicsFile[:,108:386],1)) * rescale)  #integrated and put in an array
    #calibrate Y axis : énergies en eV à gauche et à droite wavelength ou Harmonics number au choix
    
    
    
    #flip 90degree to get the delays on the horizontal axis, more convenient
    HarmScan_r = np.zeros ((np.shape(HarmonicsScan)[1],np.shape(HarmonicsScan)[0]))
    for k, osef in enumerate (np.zeros(np.shape(HarmonicsScan)[0])):
        for j, osef_aussi in enumerate (np.zeros(np.shape(HarmonicsScan)[1])):
            HarmScan_r[j,k] = HarmonicsScan[k][j]
                    
                
    
    #line for CWE cutoff at max palsma frequency for a given wavelenght
    wp = np.sqrt (ct.e**2*ct.ne / (ct.me*ct.e0)) 
    lp = 2*np.pi*ct.c/wp
    Ev_max = ct.h*ct.c/lp/ct.e                











#%% figure with harm linear
#
#plt.figure()
#plt.pcolor(Harmpics[1][80:565,:], cmap='jet')
#
#plt.figure()
#plt.plot(1/Ev[80:565][::-1], list (np.sum(Harmpics[1][80:565,:],1)))
#
#
#plt.figure()
#plt.plot(1/Ev[:][::-1], list (np.sum(Harmpics[1][:,:],1))*rescale)

## if you want to rescale the itnensity
#for k, t in enumerate(fileNumbers):
#    for e in np.arange(0, np.shape(Harmpics[0])[1] ) :
#        Harmpics[k][:,e] = Harmpics[k][:,e]*rescale 

#plt.figure()
#plt.plot(1/Ev[80:565][::-1], HarmScan_r[80:565,0])
#plt.plot(1/Ev[80:565][::-1], HarmScan_r[80:565,2])

#plt.figure()
#plt.semilogy(1/Ev[::-1], HarmScan_r[:,1]);plt.ylim(2e5,2e7)
#plt.semilogy(1/Ev[::-1], HarmScan_r[:,3])

#plt.figure()
#plt.pcolor(Harmpics[1][80:565,:], cmap='jet')
#
#
#plt.figure()
#plt.plot(1/Ev[80:565][::-1], list (np.sum(Harmpics[1][80:565,:],1)))
#

#xlist = np.arange(0,282)
xlist = np.arange(0,np.shape(Harmpics)[2])
color_map= 'inferno'
s=12
#
#fig = plt.figure(figsize=(12,4))
#gs = gridspec.GridSpec(2, 4, left=0.06, top=0.88, height_ratios=[3,1], figure=fig)
#
#
#fig.add_subplot(gs[0, 0])  
#plt.pcolor(Ev[80:565][::-1], xlist, np.rot90(Harmpics[0][80:565,:], k=1, axes=(0, 1)), cmap=color_map)
#plt.xlim(Ev[80:565].min() ,35)
#
#
#fig.add_subplot(gs[1, 0]) 
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,0])
#plt.xlim(Ev[80:565].min() ,35)
#
#
#
#fig.add_subplot(gs[0, 1])  
#plt.pcolor(Ev[80:565][::-1], xlist, np.rot90(Harmpics[1][80:565,:], k=1, axes=(0, 1)), cmap=color_map)
#plt.xlim(Ev[80:565].min() ,35)
#
#fig.add_subplot(gs[1, 1]) 
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,1])
#plt.xlim(Ev[80:565].min() ,35)
#
#
#
#fig.add_subplot(gs[0, 2])  
#plt.pcolor(Ev[80:565][::-1], xlist, np.rot90(Harmpics[2][80:565,:], k=1, axes=(0, 1)), cmap=color_map)
#plt.xlim(Ev[80:565].min() ,35)
#
#fig.add_subplot(gs[1, 2]) 
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,2])
#plt.xlim(Ev[80:565].min() ,35)
#
#
#
#fig.add_subplot(gs[0, 3])  
#plt.pcolor(Ev[80:565][::-1], xlist, np.rot90(Harmpics[3][80:565,:], k=1, axes=(0, 1)), cmap=color_map)
#plt.xlim(Ev[80:565].min() ,35)
#
#fig.add_subplot(gs[1, 3]) 
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,3])
#plt.xlim(Ev[80:565].min() ,35)
#
#
#
#
#cbaxese = fig.add_axes([0.915, 0.35, 0.02, 0.55 ]) #left, bottom, width, height 
#cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
#cbaxese = plt.ylabel('MCP intensity (a.u)',  fontsize=s)
#
#
#




##%%harmonics, square fig with 4 CEP values + 2 figs with vert integ
#
#xlist = np.arange(0,282)
#color_map= 'inferno'
#s=12
#
#fig = plt.figure(figsize=(6,6))
#gs = gridspec.GridSpec(3, 2, left=0.06, top=0.88, height_ratios=[1,1,1], wspace=0.04, hspace=0.02, figure=fig)
#
#cmax = np.amax(Harmpics)
#cmax2 = np.amax(HarmScan_r)
#
#
##divergence calib : 110.16mm MCP height and we have 281 pixels, distance from target = 919mm
#
#fig.add_subplot(gs[0, 0])  
#plt.pcolor(Ev[80:565][::-1], np.arctan((xlist-141)*0.39/919.)*1000, np.rot90(Harmpics[0][80:565,:], k=1, axes=(0, 1))/cmax, cmap=color_map, vmin=0, vmax=1)
#plt.xlim(Ev[80:565].min() ,35)
#plt.tick_params(bottom='off', labelbottom='off')
#plt.text(25,40,r'CEP $\approx 0$', color='white', fontsize=s)
#plt.ylabel(r'$\alpha$ (mrad)', fontsize=s)
#plt.yticks([-40,0,40])
#
#fig.add_subplot(gs[1, 0])  
#plt.pcolor(Ev[80:565][::-1], np.arctan((xlist-141)*0.39/919.)*1000, np.rot90(Harmpics[2][80:565,:], k=1, axes=(0, 1))/cmax, cmap=color_map, vmin=0, vmax=1)
#plt.xlim(Ev[80:565].min() ,35)
#plt.tick_params(bottom='off', labelbottom='off')
#plt.text(25,40,r'CEP $\approx \pi$', color='white', fontsize=s)
#plt.ylabel(r'$\alpha$ (mrad)', fontsize=s)
#plt.yticks([-40,0,40])
#
#fig.add_subplot(gs[2, 0]) 
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,0]/cmax2, color='purple', label=r'CEP $\approx 0$')
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,2]/cmax2, color='teal', label=r'CEP $\approx \pi$')
#plt.xlim(Ev[80:565].min() ,35)
#plt.xlabel('Photon energy (eV)',  fontsize=s)
#plt.legend(fontsize=s-2)
##plt.tick_params( left='off',labelleft='off', top='off')
#plt.ylim(0,1)
#plt.ylabel('I (a.u.)', fontsize=s)
#
#
#
#
#fig.add_subplot(gs[0, 1])  
#plt.pcolor(Ev[80:565][::-1], np.arctan((xlist-141)*0.39/919.)*1000, np.rot90(Harmpics[1][80:565,:], k=1, axes=(0, 1))/cmax, cmap=color_map, vmin=0, vmax=1)
#plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
#plt.xlim(Ev[80:565].min() ,35)
#plt.text(25,40,r'CEP $\approx \pi/2$', color='white',  fontsize=s)
#
#fig.add_subplot(gs[1, 1])  
#plt.pcolor(Ev[80:565][::-1], np.arctan((xlist-141)*0.39/919.)*1000, np.rot90(Harmpics[3][80:565,:], k=1, axes=(0, 1))/cmax, cmap=color_map, vmin=0, vmax=1)
#plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
#plt.xlim(Ev[80:565].min() ,35)
#plt.text(25,40,r'CEP $\approx 3\pi/2$', color='white',  fontsize=s)
#
#fig.add_subplot(gs[2, 1]) 
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,1]/cmax2, color='purple',label=r'CEP $\approx \pi/2$')
#plt.plot(Ev[80:565][::-1], HarmScan_r[80:565,3]/cmax2, color='teal', label=r'CEP $\approx 3\pi/2$')
#plt.tick_params( left='off',labelleft='off', top='off', labeltop='off')
#plt.xlim(Ev[80:565].min() ,35)
#plt.xlabel('Photon energy (eV)',  fontsize=s)
#plt.ylim(0,1)
#plt.legend(fontsize=s-2)
#
#
#
#
#cbaxese = fig.add_axes([0.915, 0.38, 0.03, 0.5 ]) #left, bottom, width, height 
#cmin = np.amin(Harmpics)
#cmax = np.amax(Harmpics)
#cbax = plt.colorbar(cax=cbaxese,orientation='vertical', ticks=[0,0.2,0.4,0.6,0.8,1]) #, ticks=[cmin, cmin + 0.2*(cmax-cmin), cmin + 0.4*(cmax-cmin), cmin + 0.6*(cmax-cmin) , cmin + 0.8*(cmax-cmin) ,cmax])
#cbax.ax.set_ylabel('MCP intensity (a.u)',  fontsize=s)
##cbax.ax.set_yticklabels( labels=[0,0.2,0.4,0.6,0.8,1])
#
#
#
#fig.text(0,0.9, dateCode + '_shots' + str(fileNumbers))
#
#if save==1:
#    fig = plt.savefig(outdir+'\\4CEPfigure_' + str(dateCode) + '_shots' + str(fileNumbers[0]) + '_' + str(fileNumbers[1]) +'_'+str(fileNumbers[2])+'_'+str(fileNumbers[3])+'.png', dpi=1000, bbox_inches = 'tight')
##)
#
#    
#    
#    
#    
#    
#    
#    
#    
#    
#
#
#
#
#





#%% harmonics, for all CEP values
if harm==1:
    
    xlist = np.arange(0, np.shape(Harmpics)[2])#278)  #282   275)
    color_map= 'magma'
    size=15
    
    
    vmax = np.max(Harmpics)
    vmin = np.min(Harmpics)
#    
    fig = plt.figure(figsize=(3, np.size(fileNumbers)*3))
    gs = gridspec.GridSpec( np.size(fileNumbers), 2,left=0.06, top=0.88, width_ratios=[1,3], hspace=0.02, wspace=0.02, figure=fig)
    
    
    for s, value in enumerate(cep_uniq[::-1]):
        print(str(value))
        for w, q in enumerate (fileNumbers) :
            if CEPs[q - offset_shot] == value :
    
                        
                fig.add_subplot(gs[s, 0]) 
                plt.plot( -HarmScan_r[80:565,w]/np.max(HarmScan_r[80:565,w]), Ev[80:565][::-1], color='black')
                plt.ylim(Ev[80:565].min() ,28)
               # plt.tick_params(axis='y', which='both', right=True, left=False, labelright=True, labelleft = False)
                plt.ylabel('Photon energy (eV)', fontsize=size)
                
#                plt.tick_params( top=False, labeltop=False, labelsize=size)
                if s == np.size(cep_uniq) -1:
                    plt.xlabel('I (a.u)', fontsize=size)
                    plt.xlim(-1,0)
                    plt.xticks([-1,-0.5,0])
                    plt.tick_params( top=False, labeltop=False, bottom=True, labelbottom=False, labelsize=size)
                    #plt.tick_params(axis='both', labelsize=size)
                
                else:
                    plt.tick_params( top=True, labeltop=False, bottom=False, labelbottom=False, labelsize=size)
                    plt.xticks([-1,-0.5,0])
                    plt.xlim(-1,0)

                fig.add_subplot(gs[s, 1])  
                #np.arctan((xlist-141)*0.39/919.)*1000
                plt.pcolor(((xlist-138)/(138/35)),   Ev[80:565][::-1], Harmpics[w][80:565,:]/vmax, cmap=color_map)
                plt.xticks([-30,-15,0,15,30], ['-30','', '0','','30'])
                plt.ylim(Ev[80:565].min() ,28)
                plt.text(-31,26,'CEP=' + str(value), fontsize=size, color= 'white')
    
                if s != np.size(cep_uniq) -1:
                    plt.tick_params(  left=False, labelleft=False,top=False, labeltop=False, bottom=False, labelbottom=False)
                elif s == np.size(cep_uniq) -1:
                    plt.tick_params(  left=False, labelleft=False,top=False, labeltop=False)
                    plt.xlabel('Divergence (mrad)', fontsize=size) 
                    plt.tick_params(axis='both', labelsize=size)
                break
    
                
    
    cbaxese = fig.add_axes([0.915, 0.125, 0.03, 0.755 ]) #left, bottom, width, height 
    cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
    cbaxese = plt.ylabel('MCP intensity (a.u)',  fontsize=size)
    cbaxese = plt.tick_params(labelsize=size)
    
    
    fig.text(0,0.9, dateCode + '_shots' + str(fileNumbers), fontsize=size-3)
    
    
    if save ==1:
        fig =plt.savefig(outdir+'\\Harm_CEPfigure_' + str(dateCode) + '_shots' + str(fileNumbers[0]) + '_'+str(fileNumbers[-1])+'.png', dpi=200, bbox_inches = 'tight')
    
    
    
    
    
    
#%% CLEO Harm fig

    size = 20
    
    fig = plt.figure(figsize=(np.size(cep_uniq)*1.5,3.5))
    gs = gridspec.GridSpec(1, np.size(cep_uniq),left=0.06, top=0.88, wspace=0.02, figure=fig)
    
    
    for s, value in enumerate(cep_uniq):
        print(str(value))
        for w, q in enumerate (fileNumbers) :
            if CEPs[q - offset_shot] == value :
    
                        
                hhg1=fig.add_subplot(gs[ 0, s])  
                #c= hhg1.pcolor(np.arctan((xlist-141)*0.39/919.)*1000,   Ev[80:565][::-1], Harmpics[w][80:565,:]/vmax, cmap=color_map)
                c= hhg1.pcolor(((xlist-138)/(138/32)),   Ev[80:565][::-1], Harmpics[w][80:565,:]/vmax, cmap= 'magma_r')#custom_cmap)#'viridis') #'magma_r') #color_map)   # +- 43 mrad acceptance angle of the XUV gratings
                hhg1.set_ylim(Ev[80:565].min() ,33)
                hhg1.set_xlim(-30,30)
                hhg1.set_yticks((15,20,25,30))
               # hhg1.set_xticks( [ np.mean( np.arctan((xlist-141)*0.39/919.)*1000 ) ] , minor=False )
                hhg1.set_xticks( [ np.mean( ((xlist-138)/(138/32)) ) ] , minor=False )
                #plt.text(-55,34,'CEP=' + str(value), fontsize=size, color= 'white')
                hhg1.tick_params(top=False, labeltop=False)   #cep en bas
                if value - int(value)==0:                    
                    hhg1.set_xticklabels ( [str(int(value))]  ) 
                else :
                    hhg1.set_xticklabels ( [str(value)]  ) 
                hhg1.tick_params(axis='both', labelsize=size)
                
                hhg2 = hhg1.twiny()    #divergence en haut
                hhg2.set_xlim(hhg1.get_xlim())
                hhg2.tick_params(labelsize=size)
                hhg2.set_xticks( ([ -30,0,30 ]) , minor=False )
                hhg2.set_xticks([-20,-10,10,20], minor=True )
                hhg2.xaxis.labelpad=10
    
                if s !=0:
                    #plt.tick_params(  left='off', labelleft='off',top='off', labeltop='off', labelbottom='off')
                    hhg2.tick_params(  left=False, labelleft=False, labeltop=False)
                    hhg1.tick_params( left=False, labelleft=False )
                elif s == 0:
                    hhg2.set_xlabel('Divergence (mrad)', fontsize=size) 
                    hhg1.set_ylabel('Photon energy (eV)', fontsize=size)
                    hhg1.tick_params(axis='both', labelsize=size)
                break
    
                
    
#    cbaxese = fig.add_axes([0.915, 0.125, 0.01, 0.755 ]) #left, bottom, width, height 
#    cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
#    cbaxese = plt.ylabel('MCP intensity (a.u)',  fontsize=size)
    
    
    cbaxese = fig.add_axes([0.915, 0.125, 0.01, 0.755 ]) #left, bottom, width, height 
    cbax = plt.colorbar(c,cax=cbaxese,orientation='vertical', ticks=[0,0.2,0.4,0.6,0.8,1])
    cbax.ax.set_ylabel('MCP intensity (a.u)',  fontsize=size)
    cbax.ax.set_yticklabels( labels=[0,0.2,0.4,0.6,0.8,1], fontsize=size-2)
    
    #fig.text(0,0.9, dateCode + '_shots' + str(fileNumbers), fontsize=size-3)
    fig.text(0,1.2, 'HHG spectrum', fontsize=size+1)
    fig.text(0.42,-0.1, 'relative CEP (rad)', fontsize=size+1)
    
    if save ==1:
        fig =plt.savefig(outdir+'\\Harm_CEPfigure_CLEO' + str(dateCode) + '_shots' + str(fileNumbers[0]) + '_'+str(fileNumbers[-1])+'.png', dpi=200, bbox_inches = 'tight')
    
    
    
    #%% test big integrated figure
    
    

    fig = plt.figure(figsize=(np.size(cep_uniq)*1.5,5))
    
    HarmFigCEP=[]
    
    for s, value in enumerate(cep_uniq):
        print(str(value))
        for w, q in enumerate (fileNumbers) :
            if CEPs[q - offset_shot] == value :
                
    
                        
                HarmFigCEP.append(HarmScan_r[200:565,w]/vmax)
                break
#                
#    plt.pcolor(cep_uniq,Ev[80:565][::-1],   np.rot90(HarmFigCEP,-1)[:,::-1], cmap=color_map)
#    plt.ylim(Ev[80:565].min() ,35)
#    
    
    
    




#
#
################### TEST PLOT 3D num 1
#    #from mayavi import mlab
#    #E_XUV = Ev[80:565][::-1]
#    E_XUV = Ev[200:565][::-1]
#    CEP = cep_uniq
#    
#    E_mesh, CEP_mesh = np.meshgrid(E_XUV, CEP)
#    Int =  HarmFigCEP[:][::-1]
#    Int=(Int-np.min(Int))/np.max(Int)
#    
#    
#    # figure    
#    #fig = mlab.figure(bgcolor=(1,1,1))
#    fig = plt.figure(figsize=(8,8))
#    ax=fig.gca(projection='3d')
#
#    
##    #proportions des axes
##    x_scale=2 #E
##    y_scale=2.5 #cep
##    z_scale=1  #int
##    
##    scale=np.diag([x_scale, y_scale, z_scale, 1.0])
##    scale=scale*(1.0/scale.max())
##    scale[3,3]=1.0
##    
##    def short_proj():
##      return np.dot(Axes3D.get_proj(ax), scale)
##    ax.get_proj=short_proj
##
#
#
#
#    ax.set_box_aspect((8,12,4))
#
#    #test add lightsource
#        # Create light source object.
#    ls = LightSource(azdeg=0, altdeg=70, hsv_min_val=0, hsv_max_val=1, hsv_min_sat=1, hsv_max_sat=0)
#        # Shade data, creating an rgb array.
#    rgb = ls.shade(Int, plt.cm.jet)
#
#
#
#    # Plot the data
#    #ax.scatter(E_mesh, CEP_mesh,np.array(Int), color='black',marker='.', linewidth=0.1)
#    #ax.add_collection3d(np.array(Int), CEP)
#    ax.plot_surface(E_mesh, CEP_mesh,np.array(Int),rcount=10000000, ccount=10000000, alpha=1, shade=True, facecolors=rgb)
#    #ax.contour3D(E_mesh, CEP_mesh, Int,100, cmap='inferno')
#    ax.scatter(E_mesh, CEP_mesh,np.array(Int), color='black',marker='.', alpha=0.1,linewidth=0.05)
#    #ax.plot_wireframe(E_mesh, CEP_mesh,np.array(Int), cmap=cm.jet,rstride=1, cstride=1)
#    
#
#    #labels ....
#    size=12
#    ax.set_xlabel('Photon energy (eV)', fontsize=size)
#    ax.set_ylabel('relative CEP (rad)', fontsize=size)
#    ax.set_zlabel('Intensity (a.u.)', fontsize=size)
#    ax.xaxis.labelpad=10
#    ax.yaxis.labelpad=15
#    ax.zaxis.labelpad=5
#    ax.tick_params(labelsize=size, size=size)
#    labels=([0,0.2,0.4,0.6,0.8,1,1.2])
##    ax.set_zticklabels(labels,rotation=0,verticalalignment='center',horizontalalignment='left')
##    
#    ax.set_xlim(Ev[200:565].min(),Ev[200:565].max() )
#    ax.set_zlim(np.min(Int),np.max(Int))
#    ax.set_ylim(-3,3) 
##    ax.view_init(30, -70) #enovyé par mail à rod
#    #ax.view_init(30, -70)   #angle de vue pour le plot
#    ax.view_init(40,40)
#   
#    
#    fig = plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python\\3dtest.png',transparent=True, dpi=400,bbox_inches='tight')
##########################################
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#    
#    
#    
#    plt.text(-55,34,'CEP=' + str(value), fontsize=size, color= 'white')
#                
#                
#                plt.tick_params( top='off', labeltop='off')
#    
#                if s !=0:
#                    plt.tick_params(  left='off', labelleft='off',top='off', labeltop='off')
#                elif s == 0:
#                    plt.xlabel('Divergence (mrad)', fontsize=size) 
#                    plt.ylabel('Photon energy (eV)', fontsize=size)
#                    plt.tick_params(axis='both', labelsize=size)
#                break
#    
#                
#    
#    cbaxese = fig.add_axes([0.915, 0.125, 0.01, 0.755 ]) #left, bottom, width, height 
#    cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
#    cbaxese = plt.ylabel('MCP intensity (a.u)',  fontsize=size)
#    
#    
#    fig.text(0,0.9, dateCode + '_shots' + str(fileNumbers), fontsize=size-3)
#    
#    
#
#
#
#











#    
    if save ==1:
        fig =plt.savefig(outdir+'\\Harm_CEPfigure_CLEO_int' + str(dateCode) + '_shots' + str(fileNumbers[0]) + '_'+str(fileNumbers[-1])+'.png', dpi=200, bbox_inches = 'tight')
    




    
    
    
    
    #%% electrons
LANEXsph2 = []
if elec == 1:
 
#    j=0
#    k=fileNumbers[0]

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
    elif int(dateCode) >= 20200902 and int(dateCode) <= 20200922  :                                         #validated
        px_ll_lanex =  [350,467]#[350,450] #lower left corner of the LANEX : y,x
        px_ur_lanex = [925,1000] #upper right corner of the LANEX : y,x
        a1 = 0.127 ; d = 0.254 #m
        alpha = np.deg2rad(36.54) #rad
    elif int(dateCode) == 20200923:                                         #validated
        px_ll_lanex = [340,260] #lower left corner of the LANEX : y,x
        px_ur_lanex = [950,850] #upper right corner of the LANEX : y,x
        a1 = 0.0275 ; d =  0.112#m
        alpha = np.deg2rad(36.6) #rad
    elif int(dateCode) <= 20201005 and int(dateCode) >= 20201001:
        px_ll_lanex = [345,0] #lower left corner of the LANEX : y,x
        px_ur_lanex = [825,1391] #upper right corner of the LANEX : y,x
        a1 = 0.0626 ; d = 0.2015 #m
        alpha = np.deg2rad(35) #rad
    else:
        print ('NO ELEC CALIB !')
                          
                          
                          
    LANEX_width = 0.065 #0.177 #0.065 # 0.177 # 0.065 #m for sthe small one or 177mm for the big one.
    
    
    
    
    yll =  px_ll_lanex[0] ; xll =  px_ll_lanex[1]  #ROI
    yur = px_ur_lanex [0] ; xur = px_ur_lanex [1]  #ROI
    s = round(LANEX_width/(xur-xll),7)  #m/pixel   scaling
        
    # load data
    SpaElectronsFolder = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\electrons\\' 
    
    dsz=4                                             
                                                 
    
    for j, k in enumerate (fileNumbers):
#    j=0
#    k=fileNumbers[j]
#    
        print(k) 
        
                       

        SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(k) + '.tiff' ) 
        SpaElectronsFile = np.array(SpaElectronsFile)
      
        #x, y, and values in px
        px_temp = np.linspace(0,int(xur-xll)-1,xur-xll)[0:-1:dsz]
        py_temp = np.linspace(0,yur-yll-1,yur-yll)[0:-1:dsz]
        LANEXpx = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px            
                      
        #we need to get rid of hot pixels
        temp = sorted(LANEXpx.flatten())[::-1]    #array -> list, then organize from highest to lowest value
        n_idx=[]
        for c, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
            if temp[c]-temp[c+1] >= 100:
                n_idx=np.append(n_idx,c)
        if np.shape(n_idx)[0] > 0 :  
            idx = np.where(LANEXpx >= temp[int(n_idx.max())])
            for g,f in enumerate (np.arange(0,int(np.size(idx)/2))):
                LANEXpx[idx[0][f], idx[1][f]] = LANEXpx[idx[0][f]-1, idx[1][f]]  
        
        LANEXpx = LANEXpx[0:-1:dsz, 0:-1:dsz]  
        
        
        #create coordinate mesh :
        [px,py]=np.meshgrid(px_temp,py_temp)
                
                #plot if you want to check the ROI  ... replace yll yur, xll and xur with values, or check for LANEXpx or..
        if j==0:
            fig_es = plt.figure(figsize=(5,5))
            plt.pcolor(SpaElectronsFile[yll:yur,xll:xur], cmap = 'viridis')
            plt.colorbar()
            plt.xlabel('x (pixels)') ; plt.ylabel('y (pixels)')     
#    
    
    #
    #%    #%%with real dimensions in the LANEX plane, origin = intersection with target normal
        A = a1+px*s
        B = -(yur-yll)*s/2 + py*s
        LANEXm = LANEXpx / (s**2)  #now in Coulomb/m²        
    #    plt.figure(figsize=(6,5.6))
    #    plt.pcolor(a,b,LANEXm, cmap = 'viridis')
    #    plt.xlabel('a (m)')
    #    plt.ylabel('b (m)')    
    
    
    #%  cartesian coordinates x,y,z : no scaling needed here for LANEXm
        x = A*np.cos(alpha)                #along the target plane
        y = d - A*np.sin(alpha)            #along the target normal
        z = B                              #vertical       
    
    
    #%     spherical 
        r=np.zeros(np.shape(x))
        theta=np.zeros(np.shape(x))
        LANEXsph=np.zeros(np.shape(LANEXm))
        phi=np.zeros(np.shape(x))
        for n in range(len(x[:,0])):
            for m in range(len(x[0,:])):
                r[n,m] = np.sqrt( (x[n,m])**2 + (y[n,m])**2 +(z[n,m])**2 ) 
                theta[n,m] = np.arcsin(z[n,m]/r[n,m])
                phi[n,m] = np.arctan(x[n,m]/y[n,m])
                LANEXsph[n,m]=LANEXm[n,m]*r[n,m]**2*np.cos(theta[n,m])  #Coulomb/sr
                
                        
        temp_phi =  np.linspace( phi.min(), phi.max(),   np.size(px_temp)  )
        temp_theta = np.linspace(theta.min(), theta.max(),    np.size(py_temp) ) 
        mesh_phi, mesh_theta = np.meshgrid(temp_phi, temp_theta)
        
        LANEXsph2.append(sc.interpolate.griddata((phi.flatten(), theta.flatten()), LANEXsph.flatten(), (mesh_phi, mesh_theta), method = 'nearest') )
    
    
    
    #%%   plot en vertical 
    r = None
    theta = None
    LANEXsph = None
    phi = None
    SpaElectronsFile = None
    LANEXpx = None
    px_temp = None
    py_temp= None
    px = None
    py = None
    A = None
    B=None
    LANEXm = None
    x = None
    y= None
    z=None
    mesh_phi = None
    mesh_theta = None
    
    
    
    
    j=0;k=0;
    size = 15
    vmax = np.max(LANEXsph2)
    vmin = np.min(LANEXsph2)

   
    #uncommented from here today
    
    fig_elec = plt.figure(figsize=(2.7, np.size(fileNumbers)*3))#(figsize=(3.5, np.size(cep_uniq)*1.6))  #4.5 pour les figures dans ma thèse    #2 for small lanex, 5 for big
    gs = gridspec.GridSpec(np.size(cep_uniq), 1, left=0.06, top=0.88, wspace=0.04, hspace=0.02, figure=fig_elec)
    


#    for s, value in enumerate(cep_uniq[::-1]):
#        print(str(value))
#        for j,k in enumerate (fileNumbers):
#            if CEPs[k - offset_shot] == value :
#                fig_elec.add_subplot(gs[s, 0])  
#                plt.pcolor(np.rad2deg(temp_phi),-np.rad2deg(temp_theta),LANEXsph2[j], cmap = 'ocean', vmin=vmin, vmax=vmax)
#                
#                
##                ##Vertical line for spec
##                plt.axvline(55, linestyle = '-', lw=0.6, c='purple')
#
#                plt.text(31,7.5,'CEP=' + str(CEPs[k-offset_shot]), fontsize=size, color='white')   #small lanex
#               # plt.text(18,0,'CEP=' + str(CEPs[k-offset_shot]), fontsize=size-1, color='white')    # big lanex
#                if s != np.size(cep_uniq) -1:
#                    plt.xticks([30,35,40,45,50])
#                    plt.tick_params( top=False, labeltop=False, bottom=False, labelbottom=False, labelleft=False)    #left='off',labelleft='off',
#                    #plt.ylabel('Polar angle' +  r' ${\Theta}$' + ' (deg)', fontsize=size)
#                    plt.tick_params(axis='both', labelsize=size)
#                elif s == np.size(cep_uniq) -1:
#                    plt.xlabel('Azimuth angle $\phi$ (deg)', fontsize=size)
#                    plt.xticks([30,35,40,45,50])
#                    plt.ylabel('Polar angle' +  r' ${\Theta}$' + ' (deg)', fontsize=size) 
#                    plt.tick_params(axis='both', labelsize=size)
#                break
    
    cbaxese = fig_elec.add_axes([0.915, 0.125, 0.05, 0.755 ]) #left, bottom, width, height 
    cbax = plt.colorbar(cax=cbaxese,orientation='vertical', ticks=[vmin, vmin + 0.2*(vmax-vmin),vmin + 0.4*(vmax-vmin), vmin + 0.6*(vmax-vmin) , vmin + 0.8*(vmax-vmin) ,vmax])
    cbax.ax.set_ylabel('electron charge (a.u)',  fontsize=size)
    cbax.ax.set_yticklabels( labels=[0,0.2,0.4,0.6,0.8,1], fontsize=size)
    
    
    fig_elec.text(0,0.9, dateCode + '_shots' + str(fileNumbers), fontsize=size-3)
    
    
    if save==1:
    
        fig_elec = plt.savefig(outdir+'\\Elec_CEPfigure_' + str(dateCode) + '_shots' + str(fileNumbers[0]) + '_'+str(fileNumbers[-1])+'.png', dpi=200, bbox_inches = 'tight')
    
    
    
   # %%  new fig for CLEO, cmap=hsv_r, images totales pour un scan entier de cep de -3 à +3
    
    size=20
    
    fig_elec = plt.figure(figsize=(np.size(cep_uniq)*1.5,1.8))    #5 for big
    gs = gridspec.GridSpec(1,np.size(cep_uniq), left=0.06, top=0.88, wspace=0.04, hspace=0.02, figure=fig_elec)
    INT_charge=[]

#    vmax = np.max(LANEXsph2)*1
    for s, value in enumerate(cep_uniq):
        print(str(value))
        for j,k in enumerate (fileNumbers):
            if CEPs[k - offset_shot] == value :
                fig_elec1=plt.subplot(gs[ 0,s]) 
                print('plot number=' + str(j))
                c = fig_elec1.pcolor(np.rad2deg(temp_phi),-np.rad2deg(temp_theta),LANEXsph2[j], cmap = 'hsv_r', vmin=vmin, vmax=vmax)
                fig_elec1.set_xticks( [ np.mean(np.rad2deg(temp_phi)) ] , minor=False )
                INT_charge = np.append( INT_charge, np.sum(LANEXsph2[j]))
                if value - int(value)==0:                    
                    fig_elec1.set_xticklabels ( [str(int(value))]  ) 
                    
                elif value==-2.5:
                    print(str(j) + 'c est çaaaaaaaaaaa')
                else :
                    fig_elec1.set_xticklabels ( [str(value)]  ) 
                fig_elec1.tick_params(axis='both', labelsize=size)
                fig_elec2 = fig_elec1.twiny()
                fig_elec2.set_xlim(fig_elec1.get_xlim())
                fig_elec2.tick_params(labelsize=size)
#                fig_elec2.set_xticks( ([ 30,  40 ]) , minor=False )
#                fig_elec2.set_xticks([32,34,36,38,42,44,46], minor=True )
#                
#                ##LANEX closer:
                fig_elec2.set_xticks( ([ 20,  40 ]) , minor=False )
                fig_elec2.set_xticks([16,18,20,22,24,26,28,30,32,34,36,38,42,44,46,48,50], minor=True )
                


                
                ##Vertical line for spec :
                #plt.axvline(55, linestyle = '-', lw=0.4, c='purple')
                
                #write CEP value :
                #plt.text(33,7,'CEP=' + str(CEPs[k-offset_shot]), fontsize=size, color='white')   #small lanex
                #plt.text(20,15,'CEP=' + str(CEPs[k-offset_shot]), fontsize=size-1, color='white')    # big lanex
                if s != 0:
                    fig_elec1.tick_params(labelright=False, right=False, left=False, labelleft=False, labeltop=False)
                    fig_elec2.tick_params(labeltop=False)

                elif s ==0:
                    fig_elec1.set_ylabel('${\Theta}_{pol}$' + ' (deg)', fontsize=size) 
                    fig_elec2.set_xlabel('$\phi_{az}$ (deg)', fontsize=size)
                    
                    
                break
#    
    cbaxese = fig_elec.add_axes([0.915, 0.125, 0.01, 0.755 ]) #left, bottom, width, height 
    cbax = plt.colorbar(c,cax=cbaxese,orientation='vertical', ticks=[vmin, vmin + 0.2*(vmax-vmin),vmin + 0.4*(vmax-vmin), vmin + 0.6*(vmax-vmin) , vmin + 0.8*(vmax-vmin) ,vmax])
    cbax.ax.set_ylabel('charge (a.u)',  fontsize=size)
    cbax.ax.set_yticklabels( labels=[0,0.2,0.4,0.6,0.8,1], fontsize=size-2)
    
    fig_elec.text(0.42  , -0.27, 'relative CEP (rad)', fontsize=size )
    #fig_elec.text(0,0.9, dateCode + '_shots' + str(fileNumbers), fontsize=size-3)
    fig_elec.text(0,1.35, 'Electron spatial distribution', fontsize=size+1)
    
    if save==1:
    
        fig_elec = plt.savefig(outdir+'\\Elec_CEPfigure_CLEO_' + str(dateCode) + '_shots' + str(fileNumbers[0]) + '_'+str(fileNumbers[-1])+'.png', dpi=200, bbox_inches = 'tight')
    
    
   
    
#%% TESTS 3D electrons !!!!!!!!!

import mpl_toolkits.mplot3d.axes3d as axes3d


s=14
#fig,ax = plt.subplots(subplot_kw=dict(projection='3d'),figsize = (12,8))
#ax.plot_surface(
X=np.rad2deg(temp_phi[6::])
Y=-np.rad2deg(temp_theta[5::])
X,Y = np.meshgrid(X,Y)
Z=LANEXsph2[0] # ) #rstride = 1,cstride = 1,\

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, Z[5::,6::], cmap='hsv_r', vmin=vmin, vmax=vmax,linewidth=0, antialiased=False)

ax.view_init(elev=40., azim=230)                
ax.set_zlim(2e10,4.2e10)
ax.set_xlabel(r'$\phi_{az}$ (deg)', fontsize=s)
ax.set_ylabel(r'$\theta_{pol}$ (deg)', fontsize=s)
ax.set_title(r'Electron charge', fontsize=s)
ax.tick_params(labelsize=s) 
Z2=LANEXsph2[20] # ) #rstride = 1,cstride = 1,\             
plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\0-5rad.png', bbox_inches='tight')  
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, Z2[5::,6::], cmap='hsv_r', vmin=vmin, vmax=vmax,linewidth=0, antialiased=False)

ax.view_init(elev=40., azim=230)                
ax.set_zlim(2e10,4.2e10)
ax.set_ylabel(r'$\theta_{pol}$ (deg)', fontsize=s)
ax.set_xlabel(r'$\phi_{az}$ (deg)', fontsize=s)
ax.set_title(r'Electron charge', fontsize=s)  
ax.tick_params(labelsize=s) 
plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\-2-5rad.png', bbox_inches='tight')    
                
                
                
                
                
                

#%%
s=20

plt.figure(figsize=(19.5,2)); plt.scatter(cep_uniq,INT_charge/np.max(INT_charge), color='red'); plt.xlabel('relative CEP (rad)', fontsize=s); plt.ylabel('int. charge (a.u.)', fontsize=s);plt.tick_params(labelsize=s); plt.yticks((0.8,0.85,0.9,0.95,1))
plt.xticks(np.arange(-3,3.5,0.5), ('-3','-2.5','-2','-1.5','-1', '-0.5', '0', '0.5', '1', '1.5','2','2.5','3'))
plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files' + '\\' + str(dateCode) + '_' + str(fileNumbers[0]) + 'int_elec_charge.png', bbox_inches='tight' )


#%% TEST unique 20200923 data 3ps   :   tirs 20 et 18


### test tir 20200902 : tirs 0 et 21 ((325 et 346 en vrai)) # +0.5rad et -2.5rad
#s1=0
#s2=21

#####%%%
s=17
###20200923      # +1.5 et -1.5rad
#s1=18
#s2=20

#####les 4 cep de 20200626     -2.5 et 0.5radrad
#s1=3
#s2=1


#dateCode = '20200824' ; fileNumbers = np.array([341, 336])
s1=1
s2=0

fig_elec_test = plt.figure(figsize=(6.5,7)  )  #5 for big
gs = gridspec.GridSpec(2,5, wspace=0, width_ratios=[1.6,1,0.7,1.6,1], hspace=0.3)


vmax1e=np.max(LANEXsph2[s1])  #20   ou  0
vmax2e=np.max(LANEXsph2[s2])  #18   ou 21
vmaxe=np.max([vmax1e, vmax2e])

vmin1e=np.min(LANEXsph2[s1])
vmin2e=np.min(LANEXsph2[s2])
vmine=np.min([vmin1e, vmin2e])


fig_elec_test1=plt.subplot(gs[0,0:2]) 
c = plt.pcolor(np.rad2deg(temp_phi),-np.rad2deg(temp_theta),LANEXsph2[s1], cmap = 'jet', vmin=vmine, vmax=vmaxe)
#plt.title('CEP=+0.5rad',fontsize=s)  #-1.5
plt.title('CEP=-0.5rad',fontsize=s)  #-1.5
plt.tick_params(labelsize=s-1)
#plt.yticks([-10,-5,0,5,10])
#plt.xticks([30,35,40,45])
plt.xticks([35,40,45])
plt.yticks([-5,0,5])
#plt.yticks([-20,-10,0,10,20])
#plt.xticks([20,30,40,50])
plt.ylabel('${\Theta}_{pol}$' + ' (deg)', fontsize=s) 
plt.xlabel('$\phi_{az}$ (deg)', fontsize=s)

fig_elec_test2=plt.subplot(gs[0,3:5]) 
plt.pcolor(np.rad2deg(temp_phi),-np.rad2deg(temp_theta),LANEXsph2[s2], cmap = 'jet', vmin=vmine, vmax=vmaxe)
#plt.xticks([30,35,40,45])
#plt.yticks([-10,-5,0,5,10])
plt.xticks([35,40,45])
plt.yticks([-5,0,5])
#plt.yticks([-20,-10,0,10,20])
#plt.xticks([20,30,40,50])
plt.title('CEP=+2.5rad', fontsize=s) #+1.5  #-2.5
plt.tick_params(labelsize=s-1)
plt.xlabel('$\phi_{az}$ (deg)', fontsize=s)


cbaxese = fig_elec_test.add_axes([0.95, 0.565, 0.02, 0.318 ]) #left, bottom, width, height 
cbax = plt.colorbar(c,cax=cbaxese,orientation='vertical') #, ticks=[vmine, vmine + 0.2*(vmaxe-vmine),vmine + 0.4*(vmaxe-vmine), vmine + 0.6*(vmaxe-vmine) , vmine + 0.8*(vmaxe-vmine) ,vmaxe])
cbax.ax.set_ylabel('charge (a.u)',  fontsize=s-2)
#cbax.ax.set_yticklabels( labels=[0,0.2,0.4,0.6,0.8,1])






vmax1h=np.max( Harmpics[s1][80:565,:])   #20
vmax2h=np.max( Harmpics[s2][80:565,:])  #18
vmaxh=np.max([vmax1h, vmax2h])

vmin1h=np.min( Harmpics[s1][80:565,:])
vmin2h=np.min( Harmpics[s2][80:565,:])
vminh=np.min([vmin1h, vmin2h])


fig_elec_test1=plt.subplot(gs[1,0]) 
v=plt.pcolor(((xlist-138)/(138/35)),   Ev[80:565][::-1], Harmpics[s1][80:565,:], cmap=color_map, vmin=vminh, vmax=vmaxh) 
plt.ylabel('Photon energy (ev)', fontsize=s)
plt.xlabel('div (mrad)', fontsize=s)
plt.xticks((-30,-15,0,15,30), ('-30','','0','','30'))
plt.tick_params(labelsize=s)
fig_elec_test1=plt.subplot(gs[1,1]) 
plt.plot(  np.sum( Harmpics[s1][80:565,:],1)*rescale[80:565]/np.max(np.sum( Harmpics[s1][80:565,:],1)*rescale[80:565]),  Ev[80:565][::-1], color='black') 
plt.xlim(0,1)
plt.xticks((0,0.5,1))
plt.ylim(np.min( Ev[80:565][::-1]), np.max( Ev[80:565][::-1]))
plt.tick_params(left=False, labelleft=False,bottom=True, labelbottom=False,labelsize=s-1)

fig_elec_test1=plt.subplot(gs[1,3]) 
plt.pcolor(((xlist-138)/(138/35)),   Ev[80:565][::-1], Harmpics[s2][80:565,:], cmap=color_map,vmin=vminh, vmax=vmaxh) 
plt.tick_params(labelsize=s)
plt.xticks((-30,-15,0,15,30), ('-30','','0','','30'))
plt.xlabel('div (mrad)', fontsize=s)
fig_elec_test1=plt.subplot(gs[1,4]) 
plt.plot(  np.sum( Harmpics[s2][80:565,:],1)*rescale[80:565]/ np.max( np.sum( Harmpics[s2][80:565,:],1)*rescale[80:565]),  Ev[80:565][::-1], color='black')    #18
plt.xlim(0,1)
plt.xticks((0,0.5,1))
plt.ylim(np.min( Ev[80:565][::-1]), np.max( Ev[80:565][::-1]))
plt.tick_params(left=False, labelleft=False,bottom=True, labelbottom=False, labelsize=s-1)



cbaxese = fig_elec_test.add_axes([0.95, 0.125, 0.02, 0.318 ]) #left, bottom, width, height 
cbax = plt.colorbar(v,cax=cbaxese,orientation='vertical')#, ticks=[vminh, vminh + 0.2*(vmaxh-vminh),vminh + 0.4*(vmaxh-vminh), vminh + 0.6*(vmaxh-vminh) , vminh + 0.8*(vmaxh-vminh) ,vmaxh])
cbax.ax.set_ylabel('spectral intensity (a.u)',  fontsize=s-2)
#cbax.ax.set_yticklabels( labels=[0,0.2,0.4,0.6,0.8,1])


#fig_elec1.set_xticks( [ np.mean(np.rad2deg(temp_phi)) ] , minor=False )
#if value - int(value)==0:                    
#    fig_elec1.set_xticklabels ( [str(int(value))]  ) 
#else :
#    fig_elec1.set_xticklabels ( [str(value)]  ) 
#fig_elec1.tick_params(axis='both', labelsize=size)
#fig_elec2 = fig_elec1.twiny()
#fig_elec2.set_xlim(fig_elec1.get_xlim())
#fig_elec2.tick_params(labelsize=size)
#
#

fig_elec_test=plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python' + '\\test_these_fig_' +str(dateCode) +'2CEPvalues.png', dpi=200, bbox_inches='tight' )

 
#%%    ############################### Test 3D :
#            
#
#
#    
#    vmin = np.min(LANEXsph2)
#    vmax=np.max(LANEXsph2)
#    
#    X,Y = np.meshgrid (  np.rad2deg(temp_phi)  ,   -np.rad2deg(temp_theta)  )
#
#    fig_elec_3d = plt.figure(figsize=(8,8))   #figsize=plt.figaspect(2)*1
#    ax = fig_elec_3d.add_subplot(111, projection='3d')
#    
#    
#    pos1 = ax.get_position() # get the original position 
#    pos2 = [pos1.x0 + 0.3, pos1.y0 + 0.3,  pos1.width , pos1.height ] 
#    ax.set_position(pos2) # set a new position
#    #ax = p3.Axes3D(fig_elec_3d)
#    
##    x_scale=1 #
##    y_scale=4 #
##    z_scale=1  #
##    
##    scale=np.diag([x_scale, y_scale, z_scale, 1])
##    scale=scale*(1.0/(scale.max()))
##    scale[3,3]=1
##    
##    def short_proj():
##      return np.dot(Axes3D.get_proj(ax), scale)
##    ax.get_proj=short_proj
##
##    
##    
#    
#
#
#    ax.set_box_aspect((2,14,2))
#    
#    for s, value in enumerate(cep_uniq):
#        print(str(value)+'_'+str(s))
#        for j,k in enumerate (fileNumbers):
#            if CEPs[k - offset_shot] == value :
#                Z=(LANEXsph2[j]-vmin)/vmax
#                cc=ax.contourf(X, Z, Y, 100, zdir='y', offset=-value,cmap = 'hsv_r', vmin=0, vmax=1)
#
#    #ax.set_zlim((0.,12.5))   #CEP
#    size=12
#    ax.set_xlabel('Phi(deg)', fontsize=size)
#    ax.set_ylabel('relative CEP (rad)', fontsize=size)
#    ax.set_zlabel('Theta (deg)', fontsize=size)
#    ax.set_ylim3d(-3.01,3.01)
#    ax.set_xlim3d( np.min(np.rad2deg(temp_phi))  ,  np.max( np.rad2deg(temp_phi)))
#    ax.set_zlim3d(np.min( -np.rad2deg(temp_theta)  )  , np.max( -np.rad2deg(temp_theta)  )  )
#    
#    ax.set_xticks( ([ 30,  40 ]) , minor=False )
#    ax.set_xticks([32,34,36,38,42,44,46], minor=True )
#    
#    
#    ax.set_zticks( ([ -10,0,10 ]) , minor=False )
#    ax.set_zticks([-5,5], minor=True )
#    
#    
#    ax.set_yticks( ([-3,-2,-1,0,1,2,3]) , minor=False )
#    ax.set_yticks([-2.5,-1.5,-0.5,0.5,1.5,2.5], minor=True )
#    
#    
#    ax.tick_params(which='major', length=7)
#    ax.tick_params(which='minor', length=2, color='gray')
#    
#    ax.xaxis.labelpad=10
#    ax.yaxis.labelpad=35
#    ax.zaxis.labelpad=10
#    ax.tick_params(labelsize=size, size=size)
#    #ax.view_init(40, 40)   #angle de vue pour le plot sympa!
#    ax.view_init(30,45)
#    #ax.view_init(-30,70)
#    
#    #ax.view_init(30,-70)
#    #ax.invert_xaxis()
#    #ax.invert_yaxis()
#    #ax.set_ylim3d(3,-3)
#    
#    m = plt.cm.ScalarMappable(cmap=cm.hsv_r)
#    m.set_array(Y)
#    m.set_clim(0., 1.)
#    cbaxese = fig_elec_3d.add_axes([1.2, 0.65, 0.015, 0.09 ]) #left, bottom, width, height 
#    #cbax = plt.colorbar(cc,cax=cbaxese,orientation='vertical', ticks=[vmin, vmin + 0.2*(vmax-vmin),vmin + 0.4*(vmax-vmin), vmin + 0.6*(vmax-vmin) , vmin + 0.8*(vmax-vmin) ,vmax])
#    cbax = plt.colorbar(m,boundaries=np.linspace(0, 1, 100), cax=cbaxese,orientation='vertical', ticks= [0,0.2,0.4,0.6,0.8,1])
#    cbax.ax.set_ylabel('charge (a.u)',  fontsize=size)
##    cbax.ax.set_yticklabels( labels=[0,0.2,0.4,0.6,0.8,1])
#    
#    
#    
#    
#    
#    fig_elec_3d = plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python\\3dtest_elec.png', transparent=True, dpi=400,bbox_inches='tight')
#
#
#
