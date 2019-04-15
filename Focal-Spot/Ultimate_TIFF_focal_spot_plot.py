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
filepath = r'Z:\SHHG\Stefan_and_Marie\ExperimentsWithDan_February2019\20190213'
filename = 'FocalSpot0.tiff'

outdir = r'C:\Users\ouille\Desktop\codes_python_divers'
E = 2.5
tau = 4.
scaling_x =  0.028 #43./768. #um per pixel
scaling_y =   0.028 #59./1024.
h = 7.  #half the image size in µm
s = 20 #text and labels size
im = Image.open(str(filepath) +'//' +  str(filename))
imarray = np.array(im)

#Create your own color map :
custom_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'blue', 'green', 'yellow', 'orange', 'red', 'darkred'])

##check what your color map looks like :
#fig, ax = plt.subplots()
#im = ax.imshow(np.random.random((10, 10)), cmap=custom_cmap, interpolation='nearest')
#fig.colorbar(im)
#plt.show()


##Function to obtain the FWHM :
def FWHM (x, y):
    ymax = np.max(y)
    print('\n', ymax*0.5 , 'max*0.5')
    for i, value in enumerate (y):
        if value == ymax:
            index_max = i
   #         print(index_max)
    for i, value in enumerate ( y[0:index_max] ):
        if value >= 0.5*y[index_max]:
            low = i
            break
    print(y[low] , 'y[low]')
    for i, value in enumerate (y):
        if i > index_max:
            if value <= 0.5*y[index_max]:
                high = i
                break    
    print(y[high] , 'y[high]', '\n')
    return (x[high] - x[low])


##if the image is in RGB
#def rgb2gray(rgb):
#    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
#    gray = 0.3* r + 0.59* g + 0.11 * b
#    return gray




#Big figure parameters
fig = plt.figure(figsize=(7,7))
gs = fig.add_subplot(111)
gs = GridSpec(2,2, width_ratios=[4,1], height_ratios=[1,4], wspace =0, hspace=0)



#Find the hottest pixel :
coordinates = np.where(imarray == imarray.max())
x1 = coordinates[1][0] - h/scaling_x ;x2= coordinates[1][0] + h/scaling_x   ;y1=  coordinates[0][0] -  h/scaling_y   ;y2= coordinates[0][0] + h/scaling_y

                

#cut along X axis
Xcut = plt.subplot(gs[0,0])
Xcutx = np.arange(0,np.size(imarray[0])) * scaling_x
Xcuty = imarray[coordinates[0][0]]
Xcut = plt.plot(Xcutx, Xcuty,  linewidth =2, color = 'black')
Xcut = plt.xlim(x1*scaling_x,x2*scaling_x)
## with Gaussian fit to estimate the dimension :
#mean = sum(Xcutx * Xcuty) / sum(Xcuty)
#sigma = np.sqrt(sum(Xcuty * (Xcutx - mean)**2) / sum(Xcuty))
#def Gauss(x, a, x0, sigma):
#    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
#popt1,pcov = curve_fit(Gauss, Xcutx, Xcuty, p0=[max(Xcuty), mean, sigma])
X_newX = np.arange (x1*scaling_x,x2*scaling_x, 0.0001)
#Xdim = abs(round(2.3548*popt1[2], 1))
#Xcut = plt.plot(X_newX, Gauss(X_newX, *popt1), linestyle = '--', color = 'black')
Xcut = plt.tick_params(axis='both', left='off', bottom='off', labelleft='off', labelbottom='off', labelsize=12)
fx = interpolate.interp1d(Xcutx, Xcuty)    #fonction "lissée"
Xcut = plt.plot( X_newX, fx (X_newX), color = 'blue', alpha = 0.5, linestyle = '--')  #make sure the interpolation is ok
#print (FWHM (X_newX, fx (X_newX)) , 'X FWHM')




#cut along Y axis with Gaussian fit to estimate the dimension
Ycut = plt.subplot(gs[1,1])
Ycutx = np.arange(0,np.size(imarray[:,0])) * scaling_y
Ycuty = imarray[:,coordinates[1][0]]
Ycut = plt.plot(Ycuty, Ycutx, linewidth =2, color = 'black')
Ycut = plt.ylim(y1*scaling_y,y2*scaling_y)
#mean = sum(Ycutx * Ycuty) / sum(Ycuty)
#sigma = np.sqrt(sum(Ycuty * (Ycutx - mean)**2) / sum(Ycuty))
#popt,pcov = curve_fit(Gauss, Ycutx, Ycuty, p0=[max(Ycuty), mean, sigma])
X_newY = np.arange (y1*scaling_y,y2*scaling_y, 0.0001)
#Ydim = abs(round(2.3548*popt[2], 1))
#Ycut = plt.plot(Gauss(X_newY, *popt), X_newY ,linestyle = '--', color = 'black')
Ycut = plt.tick_params(axis='both', left='off', bottom='off', labelleft='off', labelbottom='off', labelsize=12)
fy = interpolate.interp1d(Ycutx, Ycuty)
#print (FWHM (X_newY, fy (X_newY)) , 'Y FWHM')
Ycut = plt.plot( fy (X_newY), X_newY, color = 'blue', alpha = 0.5, linestyle = '--')  #make sure the interpolation is ok



#focal spot image
ax1 = plt.subplot(gs[1,0])
ax1 = plt.pcolor(imarray, cmap=custom_cmap)
ax1 = plt.ylabel('y (µm)', size=s)
ax1 = plt.xlabel('x (µm)', size=s)
ax1 = plt.xlim(x1,x2)
ax1 = plt.ylim(y1,y2)
xmax=coordinates[1][0]
ymax = coordinates[0][0]
listx = [xmax - h/scaling_x ,xmax - h/2/scaling_x, xmax , xmax + h/2/scaling_x  , xmax + h/scaling_x] ; listy = [ymax - h/scaling_y ,ymax - h/2/scaling_y, ymax , ymax + h/2/scaling_y  , ymax + h/scaling_y]
ax1 = plt.xticks(listx, [-h,-h/2,0,h/2,h], size=s) ; ax1 = plt.yticks(listy, [-h,-h/2,0,h/2,h], size=s)
#plt.scatter(coordinates[1][0],coordinates[0][0], color='r')   #to see where the hottest pixel found is
#ax1 = plt.axvline(coordinates[1][0],color='r')               #to visualize the axis for the Y cut
#ax1 = plt.axhline(coordinates[0][0],color='r')               #to visualize the axis for the X cut    

                   
                   
##estimation of the size from gaussian fit :
#ax1 = plt.text(x1 + 7,y1 + 7, str(Ydim)+ ' µm x ' + str(Xdim)+ ' µm' , color='black', fontsize = 22) 

#estimation of the FWHM from interpolation :
ax1 = plt.text(x1 + 7,y1 + 7, str(np.round(FWHM (X_newX, fx (X_newX)),1))+ ' µm x ' + str(np.round(FWHM (X_newY, fy (X_newY)),1))+ ' µm' , color='black', fontsize = s+2) 



#colorbar :
cbaxes = fig.add_axes([0.915, 0.125, 0.03, 0.605]) #left, bottom, width, height
cbaxes = plt.colorbar(cax=cbaxes,orientation='vertical')
labels = np.arange(0, 1.2, 0.2)
loc = np.arange(0, imarray.max() + imarray.max()/5, imarray.max()/5)
cbaxes.set_ticks(loc)
cbaxes.set_ticklabels(labels)



## What is the peak intensity ? (Stefan's method)
#on se limite au carré affiché sur l'image pour répartir l'énergie et calculer l'intensité pic :
Joules_per_sq_cm = imarray / np.sum (imarray[y1:y2,  x1:x2])  * E * 10**(-3) / (scaling_x*10**(-4)) / (scaling_y*10**(-4)) 
intensity_methode1 = Joules_per_sq_cm / (tau * 10**(-15))
Ipeak = intensity_methode1.max()
#print('The peak intensity is = ' + str(Ipeak))
#print('with E =', str(E), 'mJ and T = ', str(tau), 'fs')
#print('using this focal spot picture from  ', str(filepath))
fig.text( 0.45, 0.90, 'I = ' + str("%.2e"%Ipeak) + ' W/cm² with E = ' + str(E) + ' mJ and T = ' + str(tau) + ' fs  ',ha="center", va="bottom", size=s-7 ,color="black"  )
fig.text(0.45, 0.96, str(filepath) + str(filename), ha="center", va="bottom", size=s-12,color="black")

#save the figure
plt.savefig(str(outdir) + '\_' + str(filename[0:np.size(filename)-6]),dpi=200,  bbox_inches='tight')
























#
## Peak intensity taking the temporal profile from the dscan into account :
#filepath = r'Z:\Laser\D-Scan\HCF\20171110'
#filename= '1200mbar_3.8mJ'
#f_pulse     =  filepath +  '\\' + filename + '_retrieved_pulse.csv'
#raw_data = np.genfromtxt(f_pulse, delimiter=',', dtype=float)  
#pulse_T = raw_data[:,0] #list of times
#pulse_Env = raw_data[:,2]   #Electric field enveloppe
#temp_I = (pulse_Env[235:278])**2
#
#    
#facteur = temp_I / np.sum(temp_I)
#time_interval = pulse_T[100] - pulse_T[99]
#intensity = Joules_per_sq_cm.max() * facteur.max() / (time_interval*10**(-15))
#
#print ( 'Peak intensity taking the temporal intensity profile into account : ' , str(intensity))
#
#
#
#plt.figure()
#plt.plot(pulse_T[235:278], facteur)