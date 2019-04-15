# -*- coding: utf-8 -*-
"""
Created on Mon May 14 10:40:15 2018

@author: ouille
"""

from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap



## FILEPATH for the focal spot image
filepath = r'Z:\Laser\Beamprofiles\DM50\20181012\Chiller 19°C\seq6_black_and_white'
filename = 'Image0911 12-40-45.tiff'

outdir = r'C:\Users\ouille\Desktop\Article_LASER\Figure7'
s = 20 #text and labels size
im = Image.open(str(filepath) +'//' +  str(filename))
imarray = np.array(im)
#imarray = imarray[:,:,0]  #pas besoin si N&B !!!! mieux !!!



#Point the center :
plt.figure()
plt.imshow(imarray, cmap='jet')
plt.title('Choose a point')
point = plt.ginput(1,10)
plt.close()
xmax = point[0][0]
ymax=point[0][1]



#Big figure parameters
fig = plt.figure()
gs = fig.add_subplot(111)
gs = GridSpec(2,2, width_ratios=[4,1], height_ratios=[1,4], wspace =0, hspace=0)

#cut along X axis
Xcut = plt.subplot(gs[0,0])
Xcutx = np.arange(0,np.size(imarray[0])) 
Xcuty = imarray[xmax]
Xcut = plt.plot(Xcutx, Xcuty,  linewidth =1, color = 'black')
Xcut = plt.xlim(0,np.size(Xcutx))
Xcut = plt.axis('off')

#cut along Y axis with Gaussian fit to estimate the dimension
Ycut = plt.subplot(gs[1,1])
Ycutx = np.arange(0,np.size(imarray[:,0])) 
Ycuty = imarray[:,ymax]
Ycut = plt.plot(Ycuty, Ycutx, linewidth =1, color = 'black')
Ycut = plt.ylim(0,np.size(Ycutx))
Ycut = plt.axis('off')

#focal spot image
ax1 = plt.subplot(gs[1,0])
ax1 = plt.pcolor(imarray, cmap='jet')
ax1 = plt.axis('off')

#show the chosen point
plt.scatter(xmax,ymax,s=1, color='black')   #to see where the hottest pixel found is
ax1 = plt.axvline(xmax,linewidth=1,color='black')               #to visualize the axis for the Y cut
ax1 = plt.axhline(ymax,linewidth=1, color='black')               #to visualize the axis for the X cut    
           
                 
                 
#save the figure
#plt.savefig(str(outdir) + '\_' + str(filename[0:np.size(filename)-6]),dpi=200,  bbox_inches='tight')
