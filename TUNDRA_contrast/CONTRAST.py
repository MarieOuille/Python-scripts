# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:36:56 2018

@author: ouille
"""

#%% librairies
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import os
import scipy



#%% indicate the file and file path (to be changed)
filepath = r'Z:\Laser\Contrast\TUNDRA\2020\20200115'
filename0='AC20200115-1841-'
save = 1  ## save = 1 if you want to save the figure as .png , 0 if not.

#%% load data 
file = str(filepath) + '/' + str(filename0) + '.dat'

#filepath2 = r'Z:\Laser\Contrast\TUNDRA\20191030'
#filename2='AC20191030-1409-'
#file2 = str(filepath2) + '/' + str(filename2) + '.dat'
          

##Read infos contained in the first lines
data1 = open(file, 'r')
line = data1.readlines()
k=[]
j=[]
for i in (np.arange(0,np.size(line))): 
    if 'END_HEADER' in line[i]:
        k=i-1
    if 'measurements/point' in line[i]:
        j.append(i)
print(line[max(j)+1 : k+1])
            
            
            
            
data0 = np.loadtxt(str(file), skiprows=k+2)
t0 = data0[:,0]
#data2 = np.loadtxt(str(file2), skiprows=20)
#t2 = data2[:,0]




#%% Figure
fig = plt.figure(figsize=[10,6])
plt.semilogy(t0, (data0[:,1])/max(data0[:,1]), label= 'original', color='navy')
#plt.semilogy(t2, (data2[:,1])/max(data2[:,1]), label= 'GRISMs', color='red')

##to plot the mirrored trace :
plt.semilogy(-t0[::-1], (data0[:,1][::-1])/max(data0[:,1]), label= 'mirrored trace', color='navy', alpha = 0.2)



#plt.xlim(-20,0.0)          ## change the x-axis limits
#plt.ylim(1e-13, 1e+0)      ## change the y-axis limits  
plt.xlabel('time (ps)')
plt.ylabel('Normalized intensity (log scale, a.u.)')
plt.legend()
fig.suptitle(str(line[max(j)+1 : k+1]) + '\n' + 'file = ' + str(filename0)  )       ## title with comments and filename
if save == 1 :
    plt.savefig(str(filepath)+ '\\' + str(filename0) + '.png'  , dpi=500)            
































#%%



#cwd = os.getcwd()









#plt.figure()
#data1 = np.loadtxt(str(filepath) + '/' + str(filename1), skiprows=21)
#t1 = data1[:,0]
#plt.semilogy(t1, (data1[:,1])/max(data1[:,1]), label ='Femtopower 20180604-#1645')
##plt.xlim(-100,100)
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend()
#plt.savefig(str(filepath)+ 'Femtopower 20180604-#1645.png'  , dpi=500)
#
#plt.figure()
#data2 = np.loadtxt(str(filepath) + '/' + str(filename2), skiprows=8)
#t2 = data2[:,0]
#plt.semilogy(t2, (data2[:,1])/max(data2[:,1]), label ='Femtopower 20180604-#1752')
##plt.xlim(-100,100)
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend()
#plt.xlim([-17,10])
##plt.text(-51, 0.000000000007, 'detection threshold', bbox=dict(boxstyle='square', fc='white', ec='white'), color = 'darkorange')
#plt.savefig(str(filepath)+ 'Femtopower 20180604-#1752.png'  , dpi=500)
#
#
#
#
#
#
##XPW
#
#plt.figure()
#filepath = r'C:\Users\ouille\Desktop\CONTRASTE\last-tundra-meas\xpw'
#filename10='AC20180606-1659-.dat'
#filename1='AC20180606-1803-.dat'
#filename12='AC20180606-1841-.dat'
#filename3='AC20180606-1902-.dat'
#
#data10 = np.loadtxt(str(filepath) + '/' + str(filename10), skiprows=14)
#t10 = data10[:,0]
#plt.figure()
#plt.semilogy(t10, (data10[:,1])/max(data10[:,1]), label= 'XPW 20180606-#1659')
##plt.xlim(-100,100)
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend()
#plt.xlim([-40,20])
#plt.savefig(str(filepath)+ 'XPW 20180606-#1659.png'  , dpi=500)
#
#
#plt.figure()
#data1 = np.loadtxt(str(filepath) + '/' + str(filename1), skiprows=14)
#t1 = data1[:,0]
#plt.semilogy(t1, (data1[:,1])/max(data1[:,1]), linewidth = 0.6, label ='XPW 20180606-#1803')
#data12 = np.loadtxt(str(filepath) + '/' + str(filename12), skiprows=21)
#t12 = data12[:,0]
#plt.semilogy(t12[250:np.size(t12)], (data12[:,1][250:np.size(t12)])/max(data12[:,1]),  linewidth = 0.6, label ='XPW 20180606-#1841')
#plt.xlim(50,400)
#plt.ylim(10**-(12), 10**(-5))
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend()
#plt.savefig(str(filepath)+ 'XPW 20180606-#1841.png'  , dpi=500)
#
#
#plt.figure()
#data12 = np.loadtxt(str(filepath) + '/' + str(filename12), skiprows=21)
#t12 = data12[:,0]
#plt.semilogy(t12, (data12[:,1])/max(data12[:,1]), label ='XPW 20180606-#1841')
#
#data3 = np.loadtxt(str(filepath) + '/' + str(filename3), skiprows=14)
#t3 = data3[:,0]
#plt.semilogy(t3, (data3[:,1])/max(data3[:,1]), label ='XPW 20180606-#1902')
##plt.xlim(-100,100)
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend()
#plt.savefig(str(filepath)+ 'XPW 20180606-#1841_XPW 20180606-#1902.png'  , dpi=500)
#
#
#
#
#
#
#
#
##CPA2
#
#
#filepath = r'C:\Users\ouille\Desktop\CONTRASTE\last-tundra-meas\cpa2'
#filename0='AC20180607-1240-.dat'
#filename1='AC20180607-1351-.dat'
#filename2='AC20180607-1459-.dat'
#
##data0 = np.loadtxt(str(filepath) + '/' + str(filename0), skiprows=14)
##t0 = data0[:,0]
##plt.figure()
##plt.semilogy(t0, (data0[:,1])/max(data0[:,1]), label= 'CPA2 20180607-1240')
##plt.xlim(-40,40)
##plt.xlabel('time (ps)')
##plt.ylabel('Normalized intensity (log scale, a.u.)')
##plt.legend()
#
##plt.figure()
##data1 = np.loadtxt(str(filepath) + '/' + str(filename1), skiprows=21)
##t1 = data1[:,0]
##plt.semilogy(t1, (data1[:,1])/max(data1[:,1]), label ='CPA2 20180607-1351')
##plt.xlim(-40,40)
##plt.xlabel('time (ps)')
##plt.ylabel('Normalized intensity (log scale, a.u.)')
##plt.legend()
##plt.ylim(10**-(12), 5*10**(0))
##plt.savefig(str(filepath)+ 'CPA2 20180607-1351_pic.png'  , dpi=500)
##
#
#plt.figure()
#data1 = np.loadtxt(str(filepath) + '/' + str(filename1), skiprows=21)
#t1 = data1[:,0]
#plt.semilogy(t1, (data1[:,1])/max(data1[:,1]), label ='CPA2 20180607-#1351')
##plt.xlim(-40,40)
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend()
#plt.savefig(str(filepath)+ 'CPA2 20180607-1351_4ns.png'  , dpi=500)
#
#plt.figure()
#data1 = np.loadtxt(str(filepath) + '/' + str(filename1), skiprows=21)
#t1 = data1[:,0]
#plt.semilogy(t1, (data1[:,1])/max(data1[:,1]), label ='CPA2 20180607-#1351')
#plt.xlim(-900,-820)
#plt.ylim(10**-(12), 10**(-8))
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend()
#plt.savefig(str(filepath)+ 'CPA2 20180607-1351_prepulse800ps.png'  , dpi=500)
#
#
#plt.figure()
#data2 = np.loadtxt(str(filepath) + '/' + str(filename2), skiprows=14)
#t2 = data2[:,0]
#plt.semilogy(t2, (data2[:,1])/max(data2[:,1]), color='navy', linewidth =0.7, label ='CPA2 20180607-#1351')
#plt.semilogy(-t2, (data2[:,1])/max(data2[:,1]), color='green', linewidth =0.5, alpha=0.33, label ='CPA2 20180607-#1351 mirrored')
#plt.xlim(-50,50)
##plt.ylim(10**-(12), 10**(-8))
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend(fontsize=7)
#plt.savefig(str(filepath)+ 'CPA2 20180607-#1351_.png'  , dpi=500)
#
#
#plt.figure()
#data2 = np.loadtxt(str(filepath) + '/' + str(filename2), skiprows=14)
#t2 = data2[:,0]
#plt.semilogy(t2, (data2[:,1])/max(data2[:,1]), color='navy', linewidth =0.7, label ='CPA2 20180607-#1351')
#plt.semilogy(-t2, (data2[:,1])/max(data2[:,1]), color='green', linewidth =0.5, alpha=0.4, label ='CPA2 20180607-#1351 mirrored')
#plt.xlim(-30,0.1)
##plt.ylim(10**-(12), 10**(-8))
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend(fontsize=7)
#plt.savefig(str(filepath)+ 'CPA2 20180607-#1351_zoom.png'  , dpi=500)
#
#
#plt.figure()
#data2 = np.loadtxt(str(filepath) + '/' + str(filename2), skiprows=14)
#t2 = data2[:,0]
#plt.semilogy(t2, (data2[:,1])/max(data2[:,1]), color='navy', linewidth =0.7, label ='CPA2 20180607-#1351')
#plt.semilogy(-t2, (data2[:,1])/max(data2[:,1]), color='green', linewidth =0.5, alpha=0.4, label ='CPA2 20180607-#1351 mirrored')
#plt.xlim(-60,-55)
#plt.ylim(10**-(12), 10**(-8))
#plt.xlabel('time (ps)')
#plt.ylabel('Normalized intensity (log scale, a.u.)')
#plt.legend(fontsize=7)
#plt.savefig(str(filepath)+ 'CPA2 20180607-1351_zoom60ps.png'  , dpi=500)
