# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:51:59 2020

@author: ouille

Attempt for the XUV spatial beam data analysis

"""


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
from scipy import interpolate
from matplotlib.colors import LinearSegmentedColormap


#%%  ##indicate file paths and load data.... things to change :
remoteDrive = r'Z:'        #file path    Z: or C: drive 
dataPath = 'SHHG\Jaismeen_Marie'      #dataPath = e.g : 'SHHG\Stefan_and_Marie\\2nd-campaign-May2019'  
outdir = r'C:\Users\ouille\Desktop\GitHub_users\analyzed\Python'


#data to work on :
dateCode = '20200928'#
file1 = 178 #137 #158 #137
file_last= 192 #154 #170#154
#dateCode = '20200928' ; file1 = 69 ; file_last = 98 #CEP scan


save = 1

fileNumbers = np.arange(file1, file_last + 1,1)
#fileNumbers=np.array([137,138,140,141,142,143,144,145,146,147,148,150,151,152,153,154,155,156,157])   #wedge insertion -3ps
#fileNumbers=np.array([158,159,160,161,162,163,164,165,166,167,168,169,170]) #wedges insertion -4ps
#fileNumbers=np.array([ 69,71,73,75,77,79,81,83,86,88,92,94,98 ])  #CEP scan
                    
                    
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
CEP = np.loadtxt(DataLog, dtype= str, skiprows = (2), usecols=(9,) ) 



ScanDelays = -PrepulseDelay[fileNumbers-offset_shot]
scan_param = 'Gradient'



PrepulseDelay = - PrepulseDelay



#%%  for the gradient scans



custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'blue', 'green', 'yellow', 'orange', 'red', 'darkred'])
custom_cmap2 = LinearSegmentedColormap.from_list('mycmap2', ['white', 'yellow', 'orange', 'red', 'purple', 'black'])


s=13
fig = plt.figure(figsize=(10.5,6))  #6
gs = gridspec.GridSpec(3, 5, left=0.06, top=0.88, wspace=0, figure=fig)   #3,5
fig.text(0.4,0.07,'XUV spatial measurement, ' + str(dateCode) + ', tir' + str(file1) + '-' + str(file_last), fontsize=s-1)


minS = 100
maxS = 7400#8900





C =round(31/920,5) #mm/px camera image, assuming the MCP is 32 mm of diameter
d = 294 #mm, distance from interaction point on the target to the MCP
px0=540 #center of the MCP along X
py0=384#along y

def conv_angle(px,t):
    return np.arctan((px-t)*C/d)*1000  #mrad


for m, k in enumerate (fileNumbers):

    
    
    
#k=415
        
    SpaXUV = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\\harmonics\\tir'+ str(k)+'.tiff'
                                     
    Spa = Image.open(SpaXUV)                                 
    Spar = np.array(Spa)                                 
    
    
    #we need to get rid of hot pixels
    temp = sorted(Spar.flatten())[::-1]    #array -> list, then organize from highest to lowest value
    n_idx=[]
    for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
        if abs(temp[i]-temp[i+1]) >= 300:
            n_idx=np.append(n_idx,i)
    idx = np.where(Spar>= temp[int(n_idx.max())])
    
    
#    PrepulseDelay=-PrepulseDelay
    
    
    for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
        Spar[idx[0][j], idx[1][j]] = (Spar[idx[0][j]+1, idx[1][j]] + Spar[idx[0][j]-1, idx[1][j]]) /2
                
    cm='jet'
            


    if m < 5:
        fig.add_subplot(gs[0, m])
        plt.imshow(Spar, cmap = cm, vmin = minS, vmax = maxS)
        plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
        plt.title(str(round(PrepulseDelay[k-1],1)) + 'ps')                  
        plt.xlim(100,900)  
        
        

    elif m >= 5 :
        if m <10:
            fig.add_subplot(gs[1, m-5])
            plt.imshow(Spar, cmap = cm, vmin = minS, vmax = maxS)
            plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
            plt.title(str(round(PrepulseDelay[k-1],1)) + 'ps')
            plt.xlim(100,900)
        elif m >= 10:
            fig.add_subplot(gs[2, m-10])
            plt.imshow(Spar, cmap = cm, vmin = minS, vmax = maxS)
            plt.title(str(round(PrepulseDelay[k-1],1)) + 'ps')
            
            
            
            
            if m== 10:
                
                ticks_xpos=[191,366,px0,714,889]
                ticks_ypos=[35,210,py0,558,733]
                plt.xticks(ticks_xpos,[int(round(p)) for p in [conv_angle(x, px0) for x in ticks_xpos]])
                plt.yticks(ticks_ypos,[int(round(p)) for p in [conv_angle(x, py0) for x in ticks_ypos]])
                plt.tick_params(labelsize=s)
                plt.xlabel('$\Phi_{az}$ (mrad)', fontsize=s)
                plt.ylabel('$\Theta_{pol}$ (mrad)', fontsize=s)
                              

                
            elif m > 10:
                plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
            
            
            plt.xlim(100,900)


cbaxese = fig.add_axes([0.915, 0.13, 0.02, 0.75 ]) #left, bottom, width, height 
cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
cbaxese = plt.ylabel('counts (a.u)',  fontsize=s)




if save == 1 :
    fig = plt.savefig(outdir+'\\'+'XUVspatial_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'.png', bbox_inches = 'tight')
















#%% WEDGES insertion scan



#
#
#
#insert = np.array([WedgeInsert[p-1] for p in fileNumbers ])
#insert_sorted =  sorted((e,i) for i,e in enumerate(insert))
#insertion_sort = [insert_sorted[q][0] for q in (np.arange(0,np.size(insert)))]
#indices_sort = [insert_sorted[q][1] for q in (np.arange(0,np.size(insert)))]
#
#
#
#
#
#
#s=14
#fig = plt.figure(figsize=(10.5,8))  #6
#gs = gridspec.GridSpec(4, 5, left=0.06, top=0.88, wspace=0, figure=fig)   #3,5
#fig.text(0.3,1,'XUV spatial measurement, ' + str(dateCode) + ', tir' + str(file1) + '-' + str(file_last), fontsize=s-1)
#
#
#minS = 100
#maxS = 6800#8900
#
#
#
#
#
#C =round(31/920,5) #mm/px camera image, using the MCP diameter = 31 mm of diameter
#d = 294 #mm, distance from interaction point on the target to the MCP
#px0=540 #center of the MCP along X
#py0=384#along y
#
#def conv_angle(px,t):
#    return np.arctan((px-t)*C/d)*1000  #mrad
#
#Spar_sum = np.zeros(np.size(indices_sort))
#Spar_max = np.zeros(np.size(indices_sort))
#                             
#
#for m, k in enumerate (indices_sort):
#
#    
##k=415
#        
#    SpaXUV = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\\harmonics\\tir'+ str(fileNumbers[k])+'.tiff'
#                                     
#    Spa = Image.open(SpaXUV)                                 
#    Spar = np.array(Spa)                                 
#    if m == 11:
#        Spa256 = Spar
#    Spar_sum[m] = np.sum(Spar)
#    Spar_max[m]=np.max(Spar)
#    
#    #we need to get rid of hot pixels
#    temp = sorted(Spar.flatten())[::-1]    #array -> list, then organize from highest to lowest value
#    n_idx=[]
#    for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
#        if abs(temp[i]-temp[i+1]) >= 300:
#            n_idx=np.append(n_idx,i)
#    idx = np.where(Spar>= temp[int(n_idx.max())])
#    for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
#        Spar[idx[0][j], idx[1][j]] = (Spar[idx[0][j]+1, idx[1][j]] + Spar[idx[0][j]-1, idx[1][j]]) /2
#                
#                
#            
#
#
#    if m < 5:
#        fig.add_subplot(gs[0, m])
#        plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
#        plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
#        plt.title(str(insertion_sort[m]) + 'mm')                  
#        plt.xlim(100,900)  
#        
#        
#
#    elif m >= 5 :
#        if m <10:
#            fig.add_subplot(gs[1, m-5])
#            plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
#            plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
#            plt.title(str(insertion_sort[m]) + 'mm')
#            plt.xlim(100,900)
#            
#        elif (m <15) :
#            fig.add_subplot(gs[2, m-10])
#            plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
#            plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
#            plt.title(str(insertion_sort[m]) + 'mm')
#            plt.xlim(100,900)
#            
#        elif m >= 15:
#            fig.add_subplot(gs[3, m-15])
#            plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
#            plt.title(str(insertion_sort[m]) + 'mm')
#       
#            if m == 15:
#                ticks_xpos=[191,366,px0,714,889]
#                ticks_ypos=[35,210,py0,558,733]
#                plt.xticks(ticks_xpos,[int(round(p)) for p in [conv_angle(x, px0) for x in ticks_xpos]])
#                plt.yticks(ticks_ypos,[int(round(p)) for p in [conv_angle(x, py0) for x in ticks_ypos]])
#                plt.tick_params(labelsize=s)
#                plt.xlabel('$\Phi_{az}$ (mrad)', fontsize=s)
#                plt.ylabel('$\Theta_{pol}$ (mrad)', fontsize=s)
#
#                
#            elif m > 15:
#                plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
#                
#            
#            
#            plt.xlim(100,900)
#
#
#cbaxese = fig.add_axes([0.915, 0.13, 0.02, 0.75 ]) #left, bottom, width, height 
#cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
#cbaxese = plt.ylabel('counts (a.u)',  fontsize=s)
#
#
#
#
#if save == 1 :
#    fig = plt.savefig(outdir+'\\'+'XUVspatial_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'.png', bbox_inches = 'tight')
#
#
#
#GDD_fs = 36.16 #fs²/mm gdd, fused silica
#best_insertion= 2.51  #mm, found by dscan
#
#def pulse_elong(ins, dt):
#   return ( int( dt * np.sqrt( 1 + (4*np.log(2) * (ins-best_insertion)*GDD_fs/dt**2)**2 ) ) )
#
#
#
#s=20
#
#plt.figure()
#plt.plot(insertion_sort, Spar_sum/np.max(Spar_sum), '.', markersize=15,color='navy')
#plt.plot(insertion_sort, Spar_sum/np.max(Spar_sum), '--', markersize=2,color='navy')
#plt.ylabel('Integrated signal level (a.u.)', fontsize=s, color='navy')
#plt.xlabel('wedge insertion (mm)', fontsize=s)
#plt.tick_params(axis='x', labelcolor='black',labelsize=s)
#plt.tick_params(axis='y', labelcolor='navy',labelsize=s)
#plt.xlim(np.min(insertion_sort), np.max(insertion_sort))
#plt.ylim(np.min(Spar_sum/np.max(Spar_sum)),1)
#plt.axvline( insertion_sort [ np.where(Spar_sum == np.max(Spar_sum))[0][0] ], color='magenta', alpha=0.4 , label=str(insertion_sort [ np.where(Spar_sum == np.max(Spar_sum))[0][0] ])+' mm')
#plt.legend(fontsize=s-2)
#ax1=plt.gca()
#ax2=ax1.twiny()
#ax2.set_xticks((1.5,2,2.5,3))
#ax2.set_xlim(np.min(insertion_sort), np.max(insertion_sort))
#ax2.set_xticklabels((int(pulse_elong(1.5, 3.7)),int( pulse_elong(2,3.7)), int(pulse_elong(2.5, 3.7)), int(pulse_elong(3, 3.7)) ), size=s)
#ax2.set_xlabel('pulse duration (fs)', fontsize=s)
##plt.twiny()
##plt.xticks(insertion_sort,durations)
##plt.twinx()
##plt.plot(insertion_sort, Spar_max, 'x', markersize=15, color='#963275')
##plt.plot(insertion_sort, Spar_max, linestyle='dotted', markersize=2, color='#963275')
##plt.ylabel('Max signal count (a.u.)', fontsize=s, color='#963275')
##plt.tick_params(axis='y',labelsize=s, labelcolor='#963275')
#plt.savefig(outdir+'\\'+'XUVspatial_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) + '_integrated_plot' +'.png', bbox_inches = 'tight')
#
#
#
#
#
#fig=plt.figure(figsize=(5,4.2))
#gs = gridspec.GridSpec(2, 2, left=0.06, top=0.88, width_ratios=[3,1], height_ratios=[1,3],wspace=0, hspace=0,figure=fig)
#fig.add_subplot(gs[1, 0])
#s=16
##plt.imshow(Spa256, cmap = 'jet', vmin = minS, vmax = maxS)
#plt.pcolor(conv_angle(np.arange(0,1024),px0), conv_angle(np.arange(0,768),py0)[::-1] ,Spa256, cmap = 'jet', vmin = minS, vmax = maxS)
#plt.scatter(conv_angle(np.where(Spa256 == np.max(Spa256))[1][0],px0), -conv_angle(np.where(Spa256 == np.max(Spa256))[0][0],py0), marker= '+', linewidth=20, color='black')
##ticks_xpos=[191,366,px0,714,889]
##ticks_ypos=[35,210,py0,558,733]
#plt.axvline(conv_angle(np.where(Spa256 == np.max(Spa256))[1][0],px0), linestyle='--', color='black')
#plt.axhline(-conv_angle(np.where(Spa256 == np.max(Spa256))[0][0],py0), linestyle='--', color='black')
##plt.xticks(ticks_xpos,[int(round(p)) for p in [conv_angle(x, px0) for x in ticks_xpos]])
##plt.yticks(ticks_ypos,[-int(round(p)) for p in [conv_angle(x, py0) for x in ticks_ypos]])
#plt.tick_params(labelsize=s)
#plt.xlabel('$\Phi_{az}$ (mrad)', fontsize=s)
#plt.ylabel('$\Theta_{pol}$ (mrad)', fontsize=s)
#
#def smooth(y, box_pts):
#    box = np.ones(box_pts)/box_pts
#    y_smooth = np.convolve(y, box, mode='same')
#    return y_smooth
#
#fig.add_subplot(gs[1, 1])
#sm=smooth(Spa256[:,np.where(Spa256 == np.max(Spa256))[1][0]]-np.min(Spa256[:,np.where(Spa256 == np.max(Spa256))[1][0]]), 130)
#plt.plot(Spa256[:,np.where(Spa256 == np.max(Spa256))[1][0]] - np.min(Spa256[:,np.where(Spa256 == np.max(Spa256))[1][0]]), conv_angle(np.arange(0,768),py0)[::-1], color='0.8')   #vertical lineout
#plt.plot(sm,conv_angle(np.arange(0,768),py0)[::-1], color='black', linewidth = 3)
#plt.xlim(xmin=0)
#plt.ylim( np.min(conv_angle(np.arange(0,768),py0)[::-1]), np.max(conv_angle(np.arange(0,768),py0)[::-1]))
##plt.axvline((np.max(sm)-np.min(sm))/2 + np.min(sm), linewidth=0.08,color='red')
#plt.tick_params(labelsize=s)
#plt.tick_params(labelleft=False, left=False, bottom=False, labelbottom=False)
##plt.title('44mrad FWHM')
#plt.text(5500,-20,'44mrad FWHM', fontsize=s-2, rotation=-90)
#
#
#fig.add_subplot(gs[0, 0])
#smx=smooth(Spa256[np.where(Spa256 == np.max(Spa256))[0][0],:]-np.min(Spa256[np.where(Spa256 == np.max(Spa256))[0][0],:]), 50)
#plt.plot(conv_angle(np.arange(0,1024),px0), Spa256[np.where(Spa256 == np.max(Spa256))[0][0],:] - np.min(Spa256[np.where(Spa256 == np.max(Spa256))[0][0],:]), color='0.8')   #vertical lineout
#plt.plot(conv_angle(np.arange(0,1024),px0),smx ,color='black', linewidth = 3)
#plt.xlim(np.min(conv_angle(np.arange(0,1024),px0)), np.max(conv_angle(np.arange(0,1024),px0)))
#plt.tick_params(labelsize=s)
#plt.tick_params(labelleft=False, left=False, bottom=False, labelbottom=False)
##plt.axhline((np.max(smx)-np.min(smx))/2 + np.min(smx), linewidth=0.08, color='red')
##plt.text(0,6000,'77mrad FWHM', rotation=90)
#plt.title('77mrad FWHM', fontsize=s-2)
#
#
#
#
#
#plt.savefig(outdir+'\\'+'XUVspatial_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) + '_one_plot_with_lineouts_256mm' +'.png', bbox_inches = 'tight')
##%% CEP scan
##
##
##
##
##
##
##
##CEP_val = np.array([float(CEP[p-1]) for p in fileNumbers ])
##CEP_sorted =  sorted((e,i) for i,e in enumerate(CEP_val))
##CEP_sort = [CEP_sorted[q][0] for q in (np.arange(0,np.size(CEP_val)))]
##indices_sort = [CEP_sorted[q][1] for q in (np.arange(0,np.size(CEP_val)))]
##
##
##
##
##
##
##s=12
##fig = plt.figure(figsize=(10.5,6))  #6
##gs = gridspec.GridSpec(3, 5, left=0.06, top=0.88, wspace=0, figure=fig)   #3,5
##fig.text(0.3,1,'XUV spatial measurement, ' + str(dateCode) + ', tir' + str(file1) + '-' + str(file_last), fontsize=s-1)
##
##
##minS = 100
##maxS = 8600#8900
##
##
##
##
##
##C =round(31/920,5) #mm/px camera image, using the MCP diameter = 31 mm of diameter
##d = 294 #mm, distance from interaction point on the target to the MCP
##px0=540 #center of the MCP along X
##py0=384#along y
##
##def conv_angle(px,t):
##    return np.arctan((px-t)*C/d)*1000  #mrad
##
##Spar_sum = np.zeros(np.size(indices_sort))
##Spar_max = np.zeros(np.size(indices_sort))
##for m, k in enumerate (indices_sort):
##
##    
##
##    
###k=415
##        
##    SpaXUV = remoteDrive + '\\' + str(dataPath) + '\\' + dateCode + '\\harmonics\\tir'+ str(fileNumbers[k])+'.tiff'
##                                     
##    Spa = Image.open(SpaXUV)                                 
##    Spar = np.array(Spa)                                 
##    Spar_sum[m] = np.sum(Spar)
##    Spar_max[m]=np.max(Spar)
##    #we need to get rid of hot pixels
##    temp = sorted(Spar.flatten())[::-1]    #array -> list, then organize from highest to lowest value
##    n_idx=[]
##    for i, value in enumerate (temp[0:10]):  #search only within the 10 highest values for >= 100 counts differences
##        if abs(temp[i]-temp[i+1]) >= 300:
##            n_idx=np.append(n_idx,i)
##    idx = np.where(Spar>= temp[int(n_idx.max())])
##    for i,j in enumerate (np.arange(0,int(np.size(idx)/2))):
##        Spar[idx[0][j], idx[1][j]] = (Spar[idx[0][j]+1, idx[1][j]] + Spar[idx[0][j]-1, idx[1][j]]) /2
##
##            
##                
##
##
##    if m < 5:
##        fig.add_subplot(gs[0, m])
##        plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
##        plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
##        plt.title(str(CEP_sort[m]) + 'rad')                  
##        plt.xlim(100,900)  
##        
##        
##
##    elif m >= 5 :
##        if m <10:
##            fig.add_subplot(gs[1, m-5])
##            plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
##            plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
##            plt.title(str(CEP_sort[m]) + 'rad')
##            plt.xlim(100,900)
###            
###        elif (m <15) :
###            fig.add_subplot(gs[2, m-10])
###            plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
###            plt.tick_params(bottom='off', left='off',labelleft='off', labelbottom='off')
###            plt.title(str(CEP_sort[m]) + 'rad')
###            plt.xlim(100,900)
##            
##        elif m >= 10:
##            fig.add_subplot(gs[2, m-10])
##            plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
##            plt.title(str(CEP_sort[m]) + 'rad')
##       
##            if m == 10:
##                ticks_xpos=[191,366,px0,714,889]
##                ticks_ypos=[35,210,py0,558,733]
##                plt.xticks(ticks_xpos,[int(round(p)) for p in [conv_angle(x, px0) for x in ticks_xpos]])
##                plt.yticks(ticks_ypos,[int(round(p)) for p in [conv_angle(x, py0) for x in ticks_ypos]])
##                plt.tick_params(labelsize=s)
##                plt.xlabel('$\Phi_{az}$ (mrad)', fontsize=s)
##                plt.ylabel('$\Theta_{pol}$ (mrad)', fontsize=s)
##
##                
##            elif m > 10:
##                plt.tick_params(bottom=False, left=False,labelleft=False, labelbottom=False)
##                
##            
##            
##            plt.xlim(100,900)
##
##
##cbaxese = fig.add_axes([0.915, 0.13, 0.02, 0.75 ]) #left, bottom, width, height 
##cbaxese = plt.colorbar(cax=cbaxese,orientation='vertical')
##cbaxese = plt.ylabel('counts (a.u)',  fontsize=s)
##
##
##
##
##if save == 1 :
##    fig = plt.savefig(outdir+'\\'+'XUVspatial_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) +'.png', bbox_inches = 'tight')
##
##
##
##
##
##
##
##
###%%
##plt.figure()
##s=16
##plt.imshow(Spar, cmap = 'jet', vmin = minS, vmax = maxS)
##ticks_xpos=[191,366,px0,714,889]
##ticks_ypos=[35,210,py0,558,733]
##plt.xticks(ticks_xpos,[int(round(p)) for p in [conv_angle(x, px0) for x in ticks_xpos]])
##plt.yticks(ticks_ypos,[int(round(p)) for p in [conv_angle(x, py0) for x in ticks_ypos]])
##plt.tick_params(labelsize=s)
##plt.xlabel('$\Phi_{az}$ (mrad)', fontsize=s)
##plt.ylabel('$\Theta_{pol}$ (mrad)', fontsize=s)
##
##
##
##
##plt.figure()
##plt.plot(CEP_sort, Spar_sum, '.', markersize=15,color='navy')
##plt.plot(CEP_sort, Spar_sum, '--', markersize=2,color='navy')
##plt.ylabel('Integrated signal level (a.u.)', fontsize=s, color='navy')
##plt.xlabel('relative CEP (rad)', fontsize=s)
##plt.tick_params(axis='x', labelcolor='black',labelsize=s)
##plt.tick_params(axis='y', labelcolor='navy',labelsize=s)
##plt.twinx()
##plt.plot(CEP_sort, Spar_max, 'x', markersize=15, color='#963275')
##plt.plot(CEP_sort, Spar_max, linestyle='dotted', markersize=2, color='#963275')
##plt.ylabel('Max signal count (a.u.)', fontsize=s, color='#963275')
##plt.tick_params(axis='y',labelsize=s, labelcolor='#963275')
##plt.xlim(-3,3)
##plt.savefig(outdir+'\\'+'XUVspatial_' + str(dateCode) + '_shots' + str(file1) + '_' + str(file_last) + '_integrated_plot' +'.png', bbox_inches = 'tight')
