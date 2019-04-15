# -*- coding: utf-8 -*-
"""
Created on Mon May 14 10:40:15 2018

@author: ouille
"""
import os
from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage
from matplotlib.animation import FuncAnimation

video = 0 #0 if you don't want to save it, 1 if you do

## FILEPATH for the focal spot image
filepath = r'Z:\Laser\Beamprofiles\DM50\20181012\Chiller 19°C\seq7-later_in_the_afternoon'
#filename = 'Image0908 12-40-44.tiff'
#outdir = r'C:\Users\ouille\Desktop\test_DM50_beam_pointing'

scaling_x = 3.75 #µm per pixel
scaling_y = 3.75
s = 16 #text and labels size

#Create your own color map :
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'blue', 'green', 'yellow', 'orange', 'red', 'darkred'])



##Function to obtain the FWHM :
def FWHM (x, y):
    ymax = np.max(y)
    #print('\n', ymax*0.5 , 'max*0.5')
    for i, value in enumerate (y):
        if value == ymax:
            index_max = i
   #         print(index_max)
    for i, value in enumerate ( y[0:index_max] ):
        if value >= 0.5*y[index_max]:
            low = i
            break
    #print(y[low] , 'y[low]')
    for i, value in enumerate (y):
        if i > index_max:
            if value <= 0.5*y[index_max]:
                high = i
                break    
   # print(y[high] , 'y[high]', '\n')
    return (x[high] - x[low])


#creation of 2 lists to save the X and Y coordinates of the hottest pixel
Xlist = []
Ylist = []


#itération pour plotter toutes les images dans le dossier :
list_ims = os.listdir(filepath) ; list_im=[]
for i, value in enumerate (list_ims):
    if 'Image' in value :
        filename = value
        print (i)
        list_im.append(filename)

##filename = 'Image0908 12-40-44.tiff'
#
        im = Image.open(str(filepath) +'//' +  str(filename))
        imarray = np.array(im)
        
        
        #Find the center of mass :                        
        for k in np.arange(0,np.size(imarray[0])):
            for j in np.arange(0,np.size(imarray[:,0])):
                if imarray [j,k] < 15:
                    imarray [j,k] = 0
                else :
                    imarray[j,k] = imarray[j,k] - 15
        
        com = ndimage.measurements.center_of_mass(imarray)
        ###plt.scatter(com[1], com[0], color = 'black')
        
        
        
#        #Big figure parameters
#        fig = plt.figure(figsize=(9,7))
#        gs = fig.add_subplot(111)
#        gs = GridSpec(2,2, width_ratios=[5,1], height_ratios=[1,4], wspace =0, hspace=0)
#
##
##        
#        #cut along X axis
#        Xcut = plt.subplot(gs[0,0])
#        Xcutx = np.arange(0,np.size(imarray[0])) * scaling_x
#        Xcuty = imarray[com[0]]
#        Xcut = plt.plot(Xcutx, Xcuty,  linewidth =2, color = 'black')
#        Xcut = plt.tick_params(axis='both', left='off', bottom='off', labelleft='off', labelbottom='off', labelsize=12)
#        fx = interpolate.interp1d(Xcutx, Xcuty)    #fonction "lissée"
#         #   Xcut = plt.plot( Xcutx, fx (Xcutx), color = 'blue', alpha = 0.5, linestyle = '--')  #make sure the interpolation is ok
#        Xcut = plt.xlim(0,Xcutx[np.size(Xcutx)-1])
#               
#        
#        #cut along Y axis with Gaussian fit to estimate the dimension
#        Ycut = plt.subplot(gs[1,1])
#        Ycutx = np.arange(0,np.size(imarray[:,0])) * scaling_y
#        Ycuty = imarray[:,com[1]]
#        Ycut = plt.plot(Ycuty, Ycutx, linewidth =2, color = 'black')
#        Ycut = plt.tick_params(axis='both', left='off', bottom='off', labelleft='off', labelbottom='off', labelsize=12)
#        fy = interpolate.interp1d(Ycutx, Ycuty)
#          #  Ycut = plt.plot( fy (Ycutx), Ycutx, color = 'blue', alpha = 0.5, linestyle = '--')  #make sure the interpolation is ok
#        Ycut = plt.ylim(0,Ycutx[np.size(Ycutx)-1])
#                               
#        
#        #focal spot image
#        ax1 = plt.subplot(gs[1,0])
#        ax1 = plt.pcolor(imarray, cmap=custom_cmap)
#        ax1 = plt.axvline(np.mean(com[1]),color='black')               #to visualize the axis for the Y cut
#        ax1 = plt.axhline(np.mean(com[0]),color='black')               #to visualize the axis for the X cut    
#        
#                           
#        #estimation of the FWHM from interpolation :
#        ax1 = plt.text( 50, 50, str(np.round(FWHM (Xcutx, fx (Xcutx)),1))+ ' µm x ' + str(np.round(FWHM (Ycutx, fy (Ycutx)),1))+ ' µm' , color='black', fontsize = s+2) 
#        
#        
#        #colorbar :
#        cbaxes = fig.add_axes([0.915, 0.125, 0.03, 0.605]) #left, bottom, width, height
#        cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
#        labels = np.arange(0, 1.2, 0.2)
#        loc = np.arange(0, imarray.max() + imarray.max()/5, imarray.max()/5)
#        cbaxes.set_ticks(loc)
#        cbaxes.set_ticklabels(labels)
#        
#        
#        
#        fig.text(0.45, 0.94, str(filepath) + str(filename), ha="center", va="bottom", size=s-5,color="black")
        #save the figure
#        plt.savefig('test_DM50_beam_pointing\\' + str(filename[0:np.size(filename)-6]),dpi=100,  bbox_inches='tight')
#        fig = plt.close()

        Xlist.append(com[1]) ; Ylist.append(com[0])
             
Xrel = (Xlist - np.mean (Xlist))*scaling_x ; Yrel = (Ylist - np.mean (Ylist))*scaling_y  #écart à la moyenne en micromètres
distance = np.sqrt (Xrel**2 + Yrel**2) 
distance_std = np.sqrt ( np.mean(distance**2) - (np.mean(distance))**2 )
print (distance_std , ' µm  rms') 

plt.figure()
plt.scatter(Xrel, Yrel)
plt.title('Position of the center of mass for 100 shots ($\sigma$ = '+ str(np.round(distance_std,1)) + ' µm rms)')
plt.xlabel('x position (µm)')
plt.ylabel('y position (µm)')
plt.savefig('test_DM50_beam_pointing\rms_2.png')





##Generate and save a video showing the beam pointing stability
if video==1:
    fig1 = plt.figure(figsize=(9,7))  
    gs = fig1.add_subplot(111)
    gs = GridSpec(2,2, width_ratios=[5,1], height_ratios=[1,4], wspace =0, hspace=0)
    
    im = Image.open(str(filepath) +'//' +  str(list_im[0]))
    imarray = np.array(im)
    for k in np.arange(0,np.size(imarray[0])):
        for j in np.arange(0,np.size(imarray[:,0])):
            if imarray [j,k] < 15:
                imarray [j,k] = 0
            else :
                imarray[j,k] = imarray[j,k] - 15
    com = ndimage.measurements.center_of_mass(imarray)   
    
#    #cut along X axis
#    Xcut = plt.subplot(gs[0,0])
#    an1, = Xcut.plot(np.arange(0,1280)  ,   imarray[com[0]], animated = True)  #plt.plot(Xcutx, Xcuty,  linewidth =2, color = 'black')
#    Xcut.set_xlim(0,1280)
#               
#    #cut along Y axis 
#    Ycut = plt.subplot(gs[1,1])
#    an2, =Ycut.plot( imarray[:,com[1]],  np.arange(0,960), animated = True) #plt.plot(Ycuty, Ycutx, linewidth =2, color = 'black')
#    Ycut.set_ylim(0,960)
#   
#  
#
##    #colorbar :
#    cbaxes = fig1.add_axes([0.915, 0.125, 0.03, 0.605]) #left, bottom, width, height
#    cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
#    labels = np.arange(0, 1.2, 0.2)
#    loc = np.arange(0, imarray.max() + imarray.max()/5, imarray.max()/5)
#    cbaxes.set_ticks(loc)
#    cbaxes.set_ticklabels(labels)
 
      #focal spot image
    pic = plt.subplot(gs[1,0])
    an3 = pic.imshow( imarray, animated = True, cmap=custom_cmap) # plt.pcolor(imarray, cmap=custom_cmap)
    an4 = pic.scatter(com[1], 960-com[0], color='black')
    
      
    def plot_func(filename):
        im = Image.open(str(filepath) +'//' +  str(filename))
        imarray = np.array(im)
          ###Find the center of mass :                        
        for k in np.arange(0,np.size(imarray[0])):
            for j in np.arange(0,np.size(imarray[:,0])):
                if imarray [j,k] < 15:
                    imarray [j,k] = 0
                else :
                    imarray[j,k] = imarray[j,k] - 15
        
        com = ndimage.measurements.center_of_mass(imarray)
        #an1.set_data(np.arange(0,1280)  ,   imarray[com[0]])
        #an2.set_data( imarray[:,com[1]],  np.arange(0,960)  )
        an3.set_array(imarray[::-1])
        an4.set_offsets([com[1], 960-com[0]])
        #return an1,
        #return an2,
        return an3,
        return an4,

    ani = FuncAnimation(fig1, plot_func, frames = list_im, interval=200)
    ani.save('test_DM50_beam_pointing\pointing_video_2.mp4', dpi = 100)
    
    
    
    