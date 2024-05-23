# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:10:32 2019

@author: ouille


This plots raw and calibrated harmonics spectrum (image + vertically integrated)

TO DO :
    y axis = divergence

"""

#%% # import libraries
import sys
from matplotlib.gridspec import GridSpec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from scipy import ndimage
if not r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' in sys.path : 
    sys.path.insert(0,r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\python_functions' )
from phys_constants import ct
#custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['black', 'white' ])
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['#BF360C','#F4511E', '#f3c1bb', 'white' ,'#d9e1f9' , '#a0b4f0','#6d8ce8','royalblue', '#2250d9' ,'#0000e6','darkblue', 'midnightblue','black' ])







#%% ## 5 things to modify :
#" point at the harmonics file you want to plot :
filepath = r'Z:\SHHG\Jaismeen_Marie\20200629\harmonics'
filename = 'tir248.dat'   #50, 46 ,65 super !! ! from 20191202 harmonics beyond max plasma freq

#image borders:
#l1=50;h1=550;l2=50;h2=350     #Dan 2nd campaign
#l1=0;h1=550;l2=80;h2=410      #data LILU December
l1=70;h1=550;l2=80;h2=410    #Jaismeen Marie

#indicate the calibration file :
Ev = np.genfromtxt(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\\NEW_Ev-calib_tir601_VerticallyIntegrated.txt', delimiter=',', dtype='float')

#save or not? 1 for yes, 0 for no
save = 0

#where?
outdir = r'Z:\SHHG\Analyzed Data\20200629'
#that's it! you're done, nothing left to do








#%% load data
im = np.loadtxt (filepath + '\\' + filename)
Ev = Ev[0:np.size(Ev)-1]

#max plasma frequency in eV
wp = np.sqrt (ct.e**2*ct.ne / (ct.me*ct.e0)) 
lp = 2*np.pi*ct.c/wp
Ev_max = ct.h*ct.c/lp/ct.e



                   





#%% #RAW HARMONICS
#together
fig = plt.figure(figsize=(7,5))
fig.suptitle('Raw harmonics spectrum (as recorded on the MCP)', fontsize = 10)
gs = fig.add_subplot(212)
gs = GridSpec(2,1, height_ratios=[1,3],hspace=0, wspace=0)


#raw image recorded on the MCP
image_rec = plt.subplot(gs[1,0])
image_rec= plt.pcolor(ndimage.rotate(im[l1:h1,l2:h2],90), cmap='inferno')   
image_rec = plt.tick_params(left='off', labelleft = 'off')


#vertically integrated
int_rec = plt.subplot(gs[0,0])
im_v = np.zeros(np.size(im[:,0]))
for i, j in enumerate (im[:,0]):
    im_v[i] = sum (im[i,l2:h2])
int_rec = plt.plot(im_v/np.max(im_v), c='navy')
int_rec = plt.xlim(l1,h1)
int_rec = plt.tick_params( bottom='off', labelbottom='off', left='off', labelleft='off', labeltop='off')
int_rec = plt.title( 'shot # ' + str(filename[3:-4]) + ' from ' + str(filepath[-18:-10])  , fontsize=8)


#colorbar :
cbaxes = fig.add_axes([0.915, 0.125, 0.03, 0.56]) #left, bottom, width, height
cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
labels = (0,0.2,0.4,0.6,0.8,1)
loc = np.arange(0, im.max() + im.max()/5, im.max()/5)
cbaxes.set_ticks(loc)
cbaxes.set_ticklabels(labels)


#savefig
if save == 1:
    fig = plt.savefig(   str(outdir) + '\\Harm_raw_' + str(filename) +'.png'  ,dpi=400,  bbox_inches='tight'    )                
                          
                          
                          
                          
                          
                          
                

##independant plots
##
##plt.figure()
##plt.imshow(ndimage.rotate(im[0:550,80:410], 90), cmap=custom_cmap)
##plt.ylabel('vertical direction')
##plt.xlabel('photon energy')
##
#
#
###integrate vertically :
##plt.figure()
##im_v = np.zeros(np.size(im[:,0]))
##for i, j in enumerate (im[:,0]):
##    im_v[i] = sum (im[i])
##plt.plot(im_v)
##plt.xlim(0,550)









#%% #RESCALED HARMONICS
#figure parameters
fig2 = plt.figure(figsize=(7,5))
fig2.suptitle('  Calibrated harmonics spectrum', fontsize = 10)
gs2 = fig2.add_subplot(212)
gs2 = GridSpec(2,1, height_ratios=[1,3],hspace=0, wspace=0)


#vertically integrated and rescaled
int_rescaled = plt.subplot(gs2[0,0])
rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
rescale = np.append(rescale, rescale[-1]) #1element longer
int_rescaled = plt.plot(Ev[::-1], im_v*rescale/np.max(im_v*rescale)-min(im_v*rescale/np.max(im_v*rescale)), c='navy')
int_rescaled = plt.xlim(Ev[::-1][l1],np.max(Ev[::-1][l1:h1]))  
int_rescaled = plt.axvline(Ev_max,color = 'red', label = '$E_{p,max}$') 
int_rescaled = plt.tick_params( left='off', labelleft='off', bottom ='off', labelbottom='off')
int_rescaled = plt.title( 'shot # ' + str(filename[3:-4]) + ' from ' + str(filepath[-18:-10])  , fontsize=8)
int_rescaled = plt.legend(fontsize=12) 


#rescaled image
for r, value in enumerate (im[0,:]):
    im[:,r] = im[:,r]*rescale  #rescale the intensity on the image too
image_rec_rescaled = plt.subplot(gs2[1,0])
image_rec_rescaled = plt.pcolor(Ev[::-1][l1:h1], np.arange(l2,h2,1), ndimage.rotate(im[l1:h1,l2:h2], 90), cmap = 'inferno')
image_rec_rescaled = plt.axvline(Ev_max, color = 'red')
image_rec_rescaled = plt.tick_params( left='off', labelleft = 'off',labeltop='off',top='off')
image_rec_rescaled = plt.xlabel ('Photons energy (eV)')
image_rec_rescaled = plt.xticks( (15,20,25,30,35) , ('15','20','25', '30', '35') )                         
                     
                          
#colorbar :
cbaxes = fig2.add_axes([0.915, 0.125, 0.03, 0.56]) #left, bottom, width, height
cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
labels = (0,0.2,0.4,0.6,0.8,1)
loc = np.arange(0, im.max() + im.max()/5, im.max()/5)
cbaxes.set_ticks(loc)
cbaxes.set_ticklabels(labels)

                          
#savefig
if save==1:
    fig2 = plt.savefig(  str(outdir) + '\\Harm_resc_' + str(filename) +'.png'  ,dpi=400,  bbox_inches='tight'    )                
                          
                          
                          
                          
                          
                          
                          
                          


#indep plots
#avec l'axe des x rescalé 
plt.figure()
plt.pcolor(Ev[::-1][0:550], np.arange(80,410,1), ndimage.rotate(im[0:550,80:410], 90), cmap = custom_cmap)
plt.axvline(Ev_max, color = 'red')




#rescale with energy x axis in Ev
plt.figure()
rescale = abs(1/np.diff(Ev))   ##=dpx/dE  ,to rescale "dI/dpx"  to dI/dE
rescale = np.append(rescale, rescale[-1]) #1element longer
plt.plot(Ev[::-1], im_v*rescale)
plt.xlim(10.6922,np.max(Ev[::-1][0:550]))                           
plt.axvline(Ev_max,color = 'red') 