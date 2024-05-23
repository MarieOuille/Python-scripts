# -*- coding: utf-8 -*-
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
from pathlib import Path



#%%  ##indicate file paths and load data.... things to change :
remoteDrive = r'Z:'        #file path    Z: or C: drive 
dataPath = 'SHHG\Jaismeen_Marie'      #dataPath = e.g : 'SHHG\Stefan_and_Marie\\2nd-campaign-May2019'  
outdir = r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python'

save = 0    # save=0 : don't save, save = 1 : save figures
harm = 1     ;  elec = 0    ; ions = 0  ; espectro = 0

##data to work on :
#dateCode = '20200907'#'20200618'#'20190605' #'20190829'
#file1 = 33 #531 #66 
#file_last = 45 #552 #67 
#


#data to work on :
dateCode = '20200824'#'20200618'#'20190605' #'20190829'
file1 = 353 #531 #66 
file_last = 352 #552 #67 


##load calibration file for harmonics :
Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_VIntegratedharmonics_tir31.txt', delimiter=',', dtype='float')






#%%
fileNumbers = np.arange(file1, file_last + 1,1)



#load DataLog
DataLog = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\\DataLog.txt' 
                                  
ShotNo = np.loadtxt(DataLog, skiprows = (2), usecols=(0,) )
if ShotNo[0] >0:   #offset in shot num.....
    offset_shot=int(ShotNo[0])
else :
    offset_shot = 0
PulseEnergy = np.loadtxt(DataLog, skiprows = (2), usecols=(1,) )  #mJ
PulseDuration = np.loadtxt(DataLog, skiprows = (2), usecols=(2,) )  #fs
Chirp = np.loadtxt(DataLog, skiprows = (2), usecols=(3,) )  #fs²
WedgeInsert = np.loadtxt(DataLog, skiprows = (2), usecols=(4,) )  #mm
ParabolaPos = np.loadtxt(DataLog, skiprows = (2), usecols=(5,) )  #µm
PrepulseDelay = np.loadtxt(DataLog, skiprows = (2), usecols=(6,) ) #ps
SeqLen = np.loadtxt(DataLog, skiprows = (2), usecols=(7,) )   #ms


ScanDelays = -PrepulseDelay[fileNumbers-offset_shot]
scan_param = 'Gradient'
if scan_param == 'CEP':
    CEP = np.loadtxt(DataLog, dtype= str, skiprows = (2), usecols=(9,) ) 
    CEP = CEP[fileNumbers - offset_shot]
    CEP = [float(CEP[i]) for i in range(len(CEP))]
















 
                   
#%%    XUV HARMONICS           
if harm == 1 : 
    
    #ROI
    y_low = 100
    y_high = 375
    
   #load calibration for harmonics energies :
    Ev = Ev[0:np.size(Ev)]
    rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
    rescale = np.append(rescale, rescale[-1]) #1element longer
    rescale=rescale[::-1]
    
    HarmonicsScan = []
    for j, k in enumerate (fileNumbers):
        HarmonicsFolder =  remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\harmonics\\' 
        HarmonicsFile = np.loadtxt (HarmonicsFolder + 'tir' + str(k) + '.dat' ) 
        HarmonicsScan.append(list (np.sum(HarmonicsFile[:,y_low:y_high],1)) * rescale)  #integrated and put in an array
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
    

#%%   #Check ROI ?
plt.pcolor(HarmonicsFile, cmap='inferno')
plt.axvline(y_low)
plt.axvline(y_high)
#plt.figure()
#plt.plot(Ev[::-1], np.sum(HarmonicsFile[:,y_low:y_high],1)) ; plt.xlim(10,40) ; plt.title('uncalibrated intensity')
#plt.figure()
#plt.plot(Ev[::-1], list (np.sum(HarmonicsFile[:,y_low:y_high],1)) *rescale ); plt.xlim(10,40) ; plt.title('calibrated intensity')
#
#
#








#%% electrons 
ElecScan = []
if elec == 1:
    for j, k in enumerate (fileNumbers):

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
            
                              
                              
        LANEX_width = 0.065 #m
        
        # ROI and remove hot pixels                
        yll =  px_ll_lanex[0]
        xll =  px_ll_lanex[1]
        yur = px_ur_lanex [0]
        xur = px_ur_lanex [1]
        s = round(LANEX_width/(xur-xll),7)  #m/pixel
        
        # load data
        SpaElectronsFolder = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\electrons\\' 
        SpaElectronsFile = Image.open (SpaElectronsFolder + 'tir' + str(k) + '.tiff' ) 
        SpaElectronsFile = np.array(SpaElectronsFile)
        
        #x, y, and values in px
        px_temp = np.linspace(0,int(xur-xll)-1,xur-xll)
        py_temp = np.linspace(0,yur-yll-1,yur-yll)
        LANEXpx = SpaElectronsFile[yll:yur,xll:xur]   #only take the ROI, in Coulomb/px            
                      
        #we need to get rid of hot pixels
        temp = sorted(LANEXpx.flatten())[::-1]    #array -> list, then organize from highest to lowest value
        n_idx=[]
        for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
            if temp[i]-temp[i+1] >= 100:
                n_idx=np.append(n_idx,i)
        idx = np.where(LANEXpx >= temp[int(n_idx.max())])
        for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
            LANEXpx[idx[0][j], idx[1][j]] = (LANEXpx[idx[0][j]+1, idx[1][j]] + LANEXpx[idx[0][j]-1, idx[1][j]]) /2   #interpolate 
        
        #create coordinate mesh :
        [px,py]=np.meshgrid(px_temp,py_temp)
        
#        #plot if you want to check the ROI
#        fig_es = plt.figure(figsize=(6,5.6))
#        plt.pcolor(LANEXpx, cmap = 'inferno')
#        plt.colorbar()
#        plt.xlabel('x (pixels)') ; plt.ylabel('y (pixels)')     

    
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
        
        LANEXsph2 = sc.interpolate.griddata((phi.flatten(), theta.flatten()), LANEXsph.flatten(), (mesh_phi, mesh_theta), method = 'nearest') 
    
    #    plt.figure(figsize=(6,5.6))
    #    plt.pcolor(np.rad2deg(temp_phi),np.rad2deg(temp_theta),LANEXsph2, cmap = 'inferno')
    #    plt.xlabel('Azimuth angle $\phi$ (deg)')
    #    plt.ylabel('Polar angle' +  r' ${\Theta}$' + ' (deg)')
    #    plt.colorbar()
    
    
    #    plt.figure()
    #    plt.plot(np.rad2deg(temp_phi), sum( LANEXsph2))
    
        print(k)
        ElecScan.append(list(sum(LANEXsph2)))
        
            #flip 90degree to get the delays on the horizontal axis, more convenient
    ElecScan_r = np.zeros ((np.shape(ElecScan)[1],np.shape(ElecScan)[0]))
    for k, osef in enumerate (np.zeros(np.shape(ElecScan)[0])):
        for j, osef_aussi in enumerate (np.zeros(np.shape(ElecScan)[1])):
            ElecScan_r[j,k] = ElecScan[k][j]


#
#    fig_e = plt.figure(figsize=(9,3))
#    fig_e.suptitle ( 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=10 )
#    plt.pcolor(np.append(ScanDelays,ScanDelays[-1]+1), np.rad2deg(temp_phi), ElecScan_r/np.max(ElecScan_r), cmap = 'inferno')
#    plt.ylabel('$\phi_{elec}$ (deg)')
#    plt.xlabel('Prepulse delay (ps)')
#
#




















#%% ions 
ions_scan = []
if ions==1:

#    Cutoff_Ep = []
    for j, k in enumerate (fileNumbers):
#    k = 28
        ionsfilepath = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\ions\\tir' + str(k) + '_spectrum.txt'  
        if Path(ionsfilepath).is_file()== True : 
            E = np.loadtxt(ionsfilepath, skiprows = (0), usecols=(0,) )
            dNdE = np.loadtxt(ionsfilepath, skiprows = (0), usecols=(1,) )
            ions_scan.append(dNdE)
        
        
        
#        
#        #remarque ; il faudrait smoother les courbes d'abord !
#            max_i = np.where(dNdE == dNdE.max() )
#            for n, m in enumerate (E):
#                if n > max_i[0][0] + 10 :
#                    if dNdE[n] <= dNdE.min() + 0.1* (  dNdE.max() - dNdE.min() ):
#                        Cutoff_Ep.append(m) 
#                        print(str(k)+ ' m 0.02 = ' + str(m))
#                        break
##                    elif dNdE[n] <= dNdE.min() + 0.5* (  dNdE.max() - dNdE.min() ):
##                        Cutoff_Ep.append(m) 
##                        print(str(k)+ ' m 0.5 = ' + str(m))
##                        break
#
        else :
            ions_scan.append ( np.zeros(np.size(ions_scan[0])) )
            #Cutoff_Ep.append(0)
            #print(str(k), 'done 0')
        
            #flip 90degree to get the delays on the horizontal axis, more convenient
    IonsScan_r = np.zeros ((np.shape(ions_scan)[1],np.shape(ions_scan)[0]))
    for k, osef in enumerate (np.zeros(np.shape(ions_scan)[0])):
        for j, osef_aussi in enumerate (np.zeros(np.shape(ions_scan)[1])):
            IonsScan_r[j,k] = ions_scan[k][j]


    






























#%%    Final figure
fig = plt.figure(figsize=(10,10))
textsize = 13
fig.text (0.2,0.98, 'date : ' + str(dateCode) + ' , shots : #' + str(file1) + ' to #' + str(file_last), fontsize=textsize -1)



if harm==1:
    #harmonics
    gs = gridspec.GridSpec(3, 1, left=0.06, top=0.88, height_ratios=[3,1,1], figure=fig)
    fig.add_subplot(gs[0, 0])
    plt.pcolor(np.append(ScanDelays,ScanDelays[-1]+(ScanDelays[-1] - ScanDelays[-2])  ), Ev[::-1], HarmScan_r, cmap = 'jet')
    plt.ylabel('Photon energy (eV)', fontsize=textsize )
    plt.ylim(ymin=10, ymax=35)
    plt.xticks(ScanDelays+(ScanDelays[1]-ScanDelays[0])/2  , [round(ScanDelays[i],1) for i in range(len(ScanDelays))] )
    plt.axhline(y=Ev_max, color='red', label='$E_{p,max}$')
    plt.legend()
    if elec ==0 and ions==0 :
        plt.xlabel('Prepulse lead (ps)')
    cbaxesh = fig.add_axes([0.915, 0.48, 0.03, 0.4 ]) #left, bottom, width, height 
    cbaxesh = plt.colorbar(cax=cbaxesh,orientation='vertical')
    cbaxesh = plt.ylabel('Spectral Intensity (a.u)',  fontsize=textsize )

    
    if elec==1:
        #electrons
        fig.add_subplot(gs[1, 0])    
        plt.pcolor(np.append(ScanDelays,ScanDelays[-1]+(ScanDelays[-1] - ScanDelays[-2])  ), np.rad2deg(temp_phi), ElecScan_r, cmap = 'viridis')
        plt.ylabel('$\phi_{elec}$ (deg)',  fontsize=textsize)
        if ions==0:
            plt.xlabel('Prepulse lead (ps)',  fontsize=textsize)
        plt.xticks(ScanDelays+(ScanDelays[1]-ScanDelays[0])/2  , [round(ScanDelays[i],1) for i in range(len(ScanDelays))] )
        cbaxese = fig.add_axes([0.915, 0.3, 0.03, 0.135 ]) #left, bottom, width, height 
        cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
        cbaxese = plt.ylabel('int. charge (a.u)',  fontsize=textsize )
    
    
    if ions==1:
        #protons
        if elec==1:
            fig.add_subplot(gs[2,0])
        elif elec==0:
            fig.add_subplot(gs[1,0])
        #plt.scatter( ScanDelays+(ScanDelays[1]-ScanDelays[0])/2, Cutoff_Ep, zorder=1) 
#        plt.scatter( ScanDelays+(ScanDelays[1]-ScanDelays[0])/2,  [ np.sum( IonsScan_r[0:np.where(IonsScan_r[:,h]==IonsScan_r[:,h].min())[0][0] , h]  )  ] - IonsScan_r[:,h].min() )    ,  for h in enumerate ( np.arange(0,np.size(fileNumbers))  ) ], zorder=1 )
        
        
        plt.pcolor(  np.append(ScanDelays,ScanDelays[-1]+(ScanDelays[-1] - ScanDelays[-2])  )   ,   E,   IonsScan_r,  cmap='inferno', zorder=-1)       
        plt.xlabel('Prepulse lead (ps)',  fontsize=textsize)
        plt.ylabel('$E_{protons}$ (keV)'  , fontsize=textsize)
        plt.ylim(E.min(),250)
        plt.xticks(ScanDelays+(ScanDelays[1]-ScanDelays[0])/2  , [round(ScanDelays[i],1) for i in range(len(ScanDelays))] )
        if elec == 1:
            cbaxesi = fig.add_axes([0.915, 0.125, 0.03, 0.130 ]) #left, bottom, width, height 
        elif elec == 0:
            cbaxesi = fig.add_axes([0.915, 0.304, 0.03, 0.130 ]) #left, bottom, width, height 
        cbaxesi = plt.colorbar(cax=cbaxesi,orientation='vertical')
        cbaxesi = plt.ylabel('spectrum (a.u)',  fontsize=textsize )
        cbaxesi.formatter.set_powerlimits((0, 0))
        cbaxesi.ax.yaxis.get_offset_text().set_fontsize(90)


    



if save == 1 :
    fig = plt.savefig(outdir+'\\'+scan_param+'Scan_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'.png', bbox_inches = 'tight')
