# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 13:56:52 2018

@author: ouille
"""


import os
from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap



## FILEPATH for the focal spot image
filepath = r'Z:\Laser\Beamprofiles\DM50\20181012\Chiller 19°C\seq6_black_and_white'
filename = 'Image0908 12-40-44.tiff'


scaling_x = 3.75 #um per pixel
scaling_y = 3.75
s = 16 #text and labels size

#Create your own color map :
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'blue', 'green', 'yellow', 'orange', 'red', 'darkred'])


im = Image.open(str(filepath) +'//' +  str(filename))
imarray = np.array(im)
        
#sort of delete the background (not in a proper way though)        
for i in np.arange(0,np.size(imarray[0])):
    for j in np.arange(0,np.size(imarray[:,0])):
        if imarray [j,i] < 15:
            imarray [j,i] = 0
        else :
            imarray[j,i] = imarray[j,i] - 15

plt.figure()
plt.pcolor(imarray, cmap = custom_cmap)


# point the center of mass
from scipy import ndimage
coo = ndimage.measurements.center_of_mass(imarray)
plt.scatter(coo[1], coo[0], color = 'black')