# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:36:56 2018

@author: ouille
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import os
import scipy

cwd = os.getcwd()



#FEMTOPOWER

filepath = r'Z:\Laser\TUNDRA\20190612-SN2-CC'
filename0='AC20190612-1538-'
file = str(filepath) + '/' + str(filename0) + '.dat'

data0 = np.loadtxt(str(file), skiprows=20)
t0 = data0[:,0]
fig = plt.figure(figsize=[10,6])
plt.semilogy(t0, (data0[:,1])/max(data0[:,1]), label= 'original trace', color='navy')
#plt.semilogy(t0, (data0[:,1][::-1])/max(data0[:,1]), label= 'mirrored trace', color='navy', alpha = 0.2)
#plt.xlim(-900,-800)
#plt.ylim(1e-13, 1e-6)
plt.xlabel('time (ps)')
plt.ylabel('Normalized intensity (log scale, a.u.)')
plt.legend()
fig.suptitle('Compression Chamber, 0 mbar' + '\n' + 'file = ' + str(filename0))
plt.savefig(str(filepath)+ '\\' + str(filename0) + '_plot.png'  , dpi=500)




#Read infos contained in the first lines
data1 = open(file, 'r')
line = data1.readlines()


















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
