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
filepath = r'Z:\Laser\Beamprofiles\DM50\20181211'
filename = 'DM50-25A-1st-focus-bis.tiff'
outdir = r'Z:\Laser\Beamprofiles\DM50\20181211'

scaling_x = 6 #µm per pixel
scaling_y = 6
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
##plt.scatter(com[1], com[0], color = 'black')



#Big figure parameters
fig = plt.figure(figsize=(9,7))
gs = fig.add_subplot(111)
gs = GridSpec(2,2, width_ratios=[5,1], height_ratios=[1,4], wspace =0, hspace=0)

#
#        
#cut along X axis
Xcut = plt.subplot(gs[0,0])
Xcutx = np.arange(0,np.size(imarray[0])) * scaling_x
Xcuty = imarray[com[0]]
Xcut = plt.plot(Xcutx, Xcuty,  linewidth =2, color = 'black')
Xcut = plt.tick_params(axis='both', left='off', bottom='off', labelleft='off', labelbottom='off', labelsize=12)
fx = interpolate.interp1d(Xcutx, Xcuty)    #fonction "lissée"
 #   Xcut = plt.plot( Xcutx, fx (Xcutx), color = 'blue', alpha = 0.5, linestyle = '--')  #make sure the interpolation is ok
Xcut = plt.xlim(0,Xcutx[np.size(Xcutx)-1])
       

#cut along Y axis with Gaussian fit to estimate the dimension
Ycut = plt.subplot(gs[1,1])
Ycutx = np.arange(0,np.size(imarray[:,0])) * scaling_y
Ycuty = imarray[:,com[1]]
Ycut = plt.plot(Ycuty, Ycutx, linewidth =2, color = 'black')
Ycut = plt.tick_params(axis='both', left='off', bottom='off', labelleft='off', labelbottom='off', labelsize=12)
fy = interpolate.interp1d(Ycutx, Ycuty)
  #  Ycut = plt.plot( fy (Ycutx), Ycutx, color = 'blue', alpha = 0.5, linestyle = '--')  #make sure the interpolation is ok
Ycut = plt.ylim(0,Ycutx[np.size(Ycutx)-1])
                       

#focal spot image
ax1 = plt.subplot(gs[1,0])
ax1 = plt.pcolor(imarray, cmap=custom_cmap)
ax1 = plt.axvline(np.mean(com[1]),color='black')               #to visualize the axis for the Y cut
ax1 = plt.axhline(np.mean(com[0]),color='black')               #to visualize the axis for the X cut    

                   
#estimation of the FWHM from interpolation :
ax1 = plt.text( 50, 50, str(np.round(FWHM (Xcutx, fx (Xcutx)),1))+ ' µm x ' + str(np.round(FWHM (Ycutx, fy (Ycutx)),1))+ ' µm' , color='black', fontsize = s+2) 


#colorbar :
cbaxes = fig.add_axes([0.915, 0.125, 0.03, 0.605]) #left, bottom, width, height
cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
labels = np.arange(0, 1.2, 0.2)
loc = np.arange(0, imarray.max() + imarray.max()/5, imarray.max()/5)
cbaxes.set_ticks(loc)
cbaxes.set_ticklabels(labels)



fig.text(0.45, 0.94, str(filepath) + str(filename), ha="center", va="bottom", size=s-5,color="black")
#save the figure
plt.savefig(str(outdir) + '\_' + str(filename[0:np.size(filename)-6]),dpi=100,  bbox_inches='tight')
fig = plt.close()


