# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 11:57:40 2020

@author: Stefan's Matlab script for SDI analysis

TRANSLATION TO PYTHON NOT FINISHED

"""

#libraries
import numpy as np
import os
import time
from PIL import Image
from matplotlib import pyplot as plt
import sys
from operator import itemgetter


lamb = 800;  # laser wavelength (nm) for SDI (possibly selected by filter)
thetai = np.deg2rad(55); #incidence angle in target
n0overnc = 50;  # max density of target in units of nc
                # Ludo's PRX recommends to use 10 instead of 400, but that gives way too long gradients...
                # Using barrier suppression intensities and a 1e16 W/cm^2 prepulse intensity, we can estimate
                # that we should create a plasma of 102 n_cr (see evernote).
                # On the other hand, if the probe-pulse intensity is higher than that, I guess n0 would also
                # increase during the reflection. So, this values is very uncertain.
                # Neil recommends using nc=50 when comparing to his hydro-dynamics sims
                
intensityratio = 0.13; # ratio of prepulse intensity at 1st order diffraction spots and at 0th order of main pulse,
                       # measured on average on focal spots from 14/06/2019

#%%
savefig=1;
remoteDrive = 'Z:\\'



#%%
#figSavePath = r'Z:\SHHG\Analyzed Data\20200612\\'

dataPath = 'SHHG\Jaismeen_Marie\\'
#dateCode = '20200612'
#
#dataPath = r'Z:\SHHG\Stefan_and_Marie\2nd-campaign-May2019'
#dateCode='20190521'

#dataPath = r'Z:\SHHG\Jaismeen_Marie'




#
##28-fs scan:
#fileNumbers = np.arange(371,427+1) 
#brightLowerLeftPx = [290,650]  #theses pixels are to be measured 
#brightUpperRightPx = [330,690]
#darkLowerLeftPx = [340, 650]
#darkUpperRightPx = [380, 690]
##
##4-fs scan:
#fileNumbers = np.arange(172,279+1) 
#brightLowerLeftPx = [290,650]  #theses pixels are to be measured 
#brightUpperRightPx = [330,690]
#darkLowerLeftPx = [340, 650]
#darkUpperRightPx = [380, 690]

##24fs scan 
#fileNumbers = np.arange( 317,395 )


#4fs scan  20200626
dateCode='20200626'
fileNumbers = np.arange( 1001,1098 )
brightLowerLeftPx = [330,360]  #theses pixels are to be measured 
brightUpperRightPx = [370,390]
darkLowerLeftPx = [290,360]
darkUpperRightPx = [320,400]



##27fs scan   20200612
#dateCode='20200612'  #in my thesis
#fileNumbers=np.arange(371,427)
#brightLowerLeftPx = [290,660]  #theses pixels are to be measured 
#brightUpperRightPx = [330,690]
#darkLowerLeftPx = [250,615]
#darkUpperRightPx = [280,645]
#

#
### 2019.06.05 SDI
#dateCode='20190605'
#fileNumbers=np.arange(110,142)
#brightLowerLeftPx = [450,630] 
#brightUpperRightPx = [485,660]
#darkLowerLeftPx = [400, 630]
#darkUpperRightPx = [430, 660]
#
#



#brightLowerLeftPx = [345,675]  #theses pixels are to be measured 
#brightUpperRightPx = [395,710]
#
#darkLowerLeftPx = [320, 705]
#darkUpperRightPx = [360, 755]



##20200702_993_1102     3.6fs from 20200702
#dateCode='20200702'
#fileNumbers=np.arange(995,1102)
#brightLowerLeftPx = [420,660]  #theses pixels are to be measured 
#brightUpperRightPx = [455,690]
#darkLowerLeftPx = [370,660]
#darkUpperRightPx = [400,690]  #####  V, H


#%%
dataDirectory = os.path.join(remoteDrive,dataPath)
slitDirectory = os.path.join(dataDirectory, dateCode + '\slit')


DataLog = os.path.join(dataDirectory, dateCode + '\DataLog.txt')
delays =  np.loadtxt(DataLog, skiprows = (2), usecols=(6,) ) 

                       
                       
                       
     
#%% to measure the pixels for the initially bright and dark regions
Figure1 = plt.figure(figsize=(10,10))

SDIim0 = Image.open(slitDirectory + '\\tir' + str(fileNumbers[0]) + '.tiff')
SDIimarr0 = np.array(SDIim0)

plt.imshow(SDIimarr0,vmin=0, vmax=2500)
#plt.xlim(650,850)

#plt.xlim(400,900)
#plt.ylim(100,600)
##
#plt.scatter( 630,450 , linewidth=1, color='red')  ## h, v
#plt.scatter( 660,485 , linewidth=1, color='red')  ## h, v
#
#plt.scatter( 630,400 , linewidth=1, color='green')
#plt.scatter( 660,430 , linewidth=1, color='green')

plt.scatter(darkLowerLeftPx[1], darkLowerLeftPx[0]  , linewidth=1, color='red')
plt.scatter(darkUpperRightPx[1], darkUpperRightPx[0]  , linewidth=1, color='red')

plt.scatter(brightUpperRightPx[1], brightUpperRightPx[0]  , linewidth=1, color='green')
plt.scatter(brightLowerLeftPx[1], brightLowerLeftPx[0]  , linewidth=1, color='green')


 
#%% Treat data
delay = np.zeros(np.size(fileNumbers)) ; SDIContrast = np.zeros(np.size(fileNumbers))
for k, i in enumerate (fileNumbers) :
    delay[k] = -delays[i-1]
    
    temp = np.array(Image.open(slitDirectory + '\\tir' + str(i) + '.tiff'))
    
    bright = sum(temp[brightLowerLeftPx[0],brightLowerLeftPx[1]:brightUpperRightPx[1]])
    for m, n in enumerate (np.arange(brightLowerLeftPx[0]+1,brightUpperRightPx[0]+1,1)):
        bright = bright + sum (temp[n,brightLowerLeftPx[1]:brightUpperRightPx[1]])

    dark = sum(temp[darkLowerLeftPx[0],darkLowerLeftPx[1]:darkUpperRightPx[1]])
    for m, n in enumerate (np.arange(darkLowerLeftPx[0]+1,darkUpperRightPx[0]+1,1)):
        dark= dark + sum (temp[n,darkLowerLeftPx[1]:darkUpperRightPx[1]])
        
    SDIContrast[k] = (bright - dark)/(bright + dark) 

       
    print("{} shots to go ".format((fileNumbers[-1] - i)) )
    sys.stdout.flush()

    
print('\n')
    


    

#%% sort acording to delay
if dateCode == '20200612':
    delay = delay + 1.77 #added line to account for the offset on the 12th of June, 2020

#only keep positive delays
pos_delay = []; posindcs = [] ; pos_SDIContrast = []
for i, j in enumerate(delay) :
    if j >= 0:
       pos_delay = np.append( pos_delay , j)
       pos_indcs = np.append( posindcs , i)
       pos_SDIContrast = np.append( pos_SDIContrast , SDIContrast[i]) 
       
#sort :
sortedindcs, delay = zip(*sorted(enumerate(pos_delay), key=itemgetter(1))) 
SDIContrast = [pos_SDIContrast[k] for k in sortedindcs]




#%% find zero crossings and smooth

#smoothing method
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


#smoothby = 3
#[fooXindcs,fooTi]=crossings(smooth(diff(SDIContrast),smoothby),delay(1:end-1),0,'int') #inversion delays
#if fooTi(1)> 1 
#    fooXindcs = [1; fooXindcs]
#    fooTi = [0; fooTi]
#
##for noisy data this search for extrema sometmes finds several "little extrema" nearby
##whose positions should be averaged to a single one. That's what we do here:
#checkdiff = diff(sign(SDIContrast(fooXindcs))); #when this is zero, we're still sitting on the same extremum,
#                                                #when it's -2 we've gone from a max to a min, 
#                                                #when it's +2, we've gone from a min to a max
#foo = -2; #the first sign switch we're looking for is from max to min                  
#
#lastgood =1; #when averaging the positions of "little extrema", start from this indx and go until where you found the sign switch
#             #, will be reset in the loop each time we find a sign switch
#k=1; #we start writing the proper Ti vector from indx 1 on, obviously, will be incremented in the loop
#for i=1:length(fooTi)-1
#    if checkdiff (i) == foo #if you've arrived at an indx of fooTi where SDIContrast has changed its sign (first max to min, then min to max, back and forth)
#        Ti(k) = mean(fooTi(lastgood:i)); #delay for inversion is averaged over all the last "little extrema" 
#        Xindcs(k) = round(mean(fooXindcs(lastgood:i))); #same for the indices
#        k=k+1;          #increment indx for proper Ti
#        lastgood = i+1; #for next sign switch start averageing at this next 
#        foo = -foo;     #change back and forth between "max to min" and "min to max"
#
#    #if at the next indx of fooTi, SDIContrast has not changed sign, just go to the next indx of fooTi
#


#%%


s=18
smoothed=smooth(np.append(SDIContrast[0:6],SDIContrast), 7)[6::]
#smoothed=np.append(smoothed[0], smoothed)
plt.figure(figsize=(5,5))
plt.scatter(delay,SDIContrast, color='navy', label = 'data')
#plt.xlim(np.min(delay), np.max(delay))
plt.xlim(0, np.max(delay))
plt.xlabel(r'Prepulse lead $\tau$ (ps)', fontsize=s)
plt.ylabel(r'$SDI_{contrast}$', fontsize=s)
plt.tick_params (labelsize=s-3)
plt.plot(delay, smoothed, color='0.4')

sign_change=[]
signs = np.sign(smoothed)
extrema_value = []
extrema_index=[]
for j, k in enumerate (smoothed) :
    if j > 0 and j<np.size(signs)-1:
        if signs [j+1] != signs[j]:
            sign_change.append(j)
for n, p in enumerate (sign_change):
    if n==0:
        extrema_value.append(np.max(np.abs(smoothed[0:p]))   *   signs[sign_change[0]]  )       
    else :
        extrema_value.append(np.max(np.abs(smoothed[sign_change[n-1]:p]))  *   signs[sign_change[n]]   )
extrema_idx = [np.where(smoothed == i) for i in extrema_value]
extrema_delay = [  delay[ extrema_idx[i][0][0] ] for i, j in enumerate (extrema_idx) ]
plt.scatter(extrema_delay, extrema_value, color='red',linewidth=s-10, label='extrema')
plt.legend(fontsize=s-4)
plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\SDI_contrast_'+str(dateCode)+'_'+str(fileNumbers[0]) +'-'+str(fileNumbers[-1]) + '.png', bbox_inches='tight')


##smoothed curve
#smoothed=smooth(SDIContrast, 8)
#plt.plot(delay, smoothed  , color='light blue')
##plt.legend()
##plt.plot(delay(Xindcs),SDIContrast(Xindcs),'o')
##




#phase vs delay
plt.figure(figsize=(5,5))
plt.scatter(extrema_delay, np.arange(0, np.pi*np.size(extrema_delay),np.pi), color='red', label='extrema')
plt.xlim(0, np.max(delay))
plt.ylim(ymin=0)
plt.xlabel(r'Prepulse lead $\tau$ (ps)', fontsize=s)
plt.ylabel(r'$\Delta \phi (rad)$', fontsize=s)
plt.tick_params (labelsize=s-3)


#fit =  np.polyfit(extrema_delay, np.arange(0, np.pi*np.size(extrema_delay),np.pi), 1)


from scipy.optimize import curve_fit

def func(x,a):
    return a*x

fit = curve_fit(func, extrema_delay, np.arange(0, np.pi*np.size(extrema_delay),np.pi) )
fit_ex = curve_fit(func, extrema_delay[0:4], np.arange(0, np.pi*np.size(extrema_delay[0:4]),np.pi) )


# let's find cs now!!

############### NEIL PHD :
intensityratio = 0.13 # ratio of prepulse intensity at 1st order diffraction spots and at 0th order of main pulse,
                       # measured on average on focal spots from 14/06/2019
lambd0 = 800 #central wavelength (nm)
thetai = np.deg2rad(55) #deg
ne0overnc = 50   # ne0 / nc ratio
#cs =fit[0][0]* lambd0 * np.cos(thetai) / (4*np.pi) *  (  (1-np.sqrt(intensityratio)) * np.log(50 / np.exp(1)* (np.cos(thetai))**2    )**(-1))
#cs_ex =fit_ex[0][0]* lambd0 * np.cos(thetai) / (4*np.pi) *  (  (1-np.sqrt(intensityratio)) * np.log(50 / np.exp(1)* (np.cos(thetai))**2    )**(-1))

cs =fit[0][0]* lambd0 / np.cos(thetai) / (4*np.pi) *  (  (1-np.sqrt(intensityratio)) * np.log(50 / np.exp(1)* (np.cos(thetai))**2    )**(-1))
cs_ex =fit_ex[0][0]* lambd0 / np.cos(thetai) / (4*np.pi) *  (  (1-np.sqrt(intensityratio)) * np.log(50 / np.exp(1)* (np.cos(thetai))**2    )**(-1))



print('cs = ' + str(np.round(cs,2)) + 'nm/ps')



#plt.plot(  np.arange(np.min(delay), np.max(delay), 1)     ,     fit[1] + fit[0]* np.arange(np.min(delay), np.max(delay), 1)    , color='navy' , label='fit : coeff = ' + str(np.round(fit[0],3)))

plt.plot(  np.arange(0, np.max(delay), 1)     ,      fit[0][0]* np.arange(0, np.max(delay), 1)    , color='navy' , label='fit : coeff = ' + str(np.round(fit[0][0],3)))
plt.plot(  np.arange(0, np.max(delay), 1)     ,      fit_ex[0][0]* np.arange(0, np.max(delay), 1)    , color='green' , label='fit : coeff = ' + str(np.round(fit_ex[0][0],3)))
plt.legend(fontsize=s-4)
plt.text(10,np.pi*np.size(extrema_delay),'$c_s = $' + str(np.round(cs,1))+ ' nm/ps', color='navy', fontsize=s)
plt.text(10,np.pi*np.size(extrema_delay)+1.5,'$c_s = $' + str(np.round(cs_ex,1))+ ' nm/ps', color='green', fontsize=s)

plt.savefig(r'C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\SHHG_python\generated_files\SDI_linear_fit_'+str(dateCode)+'_'+str(fileNumbers[0]) +'-'+str(fileNumbers[-1]) + '.png', bbox_inches='tight')








#%% Improved/Corrected analysis according to Neïl (thesis p.116)
#excludelast =0 #how many of the last inversions to exclude from fit (often the osciallations get slower for long delays)

#dPhi = (0:1:length(Ti)-1)*np.pi
#linfit=fit(Ti(1:end-excludelast)',dPhi(1:end-excludelast)','poly1','Lower',[0,0],'Upper',[Inf,0])
#dPhidT = linfit.p1 #in rad/ps
#
#vPNeil = dPhidT *np.cos(thetai)*lam/4/np.pi * ((1-np.sqrt(intensityratio)) * np.log(n0overnc *np.exp(-1) / np.cos(thetai)^2))^(-1)
#
#fprintf ('*********************************************\n');
#fprintf ('Neïl`s plasma expansion speed is %4.2f nm/ps.\n',vPNeil);
#fprintf ('*********************************************\n')




#%%
#plotName = ['SDIScanFigure_',num2str(fileNumbers(1)),'-',num2str(fileNumbers(end))];
#f = myFigure(plotName,20,8,[]);
#
#subplot(1,2,1)
#plot(delay, SDIContrast,'o')
##plot(delay, SDIContrast,'o-')
#hold on
#plot(delay, smooth(SDIContrast,smoothby) ,'-')
#xlim([0,max(delay)])
#
#xlabel('Prepulse lead (ps)')
#ylabel('SDI Contrast')
#
#subplot(1,2,2)
#plot(Ti, dPhi,'o')
#hold on
#plot(linfit ,'-')
#
#ylim([0,0.1*ceil((10*max(dPhi)))])
#legend('experimental',['lin. fit v_p =',num2str(round(vPNeil*10)/10),' nm/ps'])
#
#xlabel('Prepulse lead (ps)')
#ylabel('Dephasing (rad)')




#%%
#if savefig==1
#    if ~exist(fullfile(figSavePath,dateCode),'dir')
#        mkdir(fullfile(figSavePath,dateCode));
#    saveFigure(f, fullfile(figSavePath,dateCode,['SDIScanFigure_',dateCode,'_shots_',num2str(fileNumbers(1)),'-',num2str(fileNumbers(end))]), 'pdf');