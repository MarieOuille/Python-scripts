# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 17:38:52 2017

@author: ouille
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import os

## FIRST : EXTRACT DATA 
#print ('Where is the file located? e.g : Z:\Laser\CEP\20190311')
#filepath = input('The path of the file is :')
#print ('What is the name of the .dat file? e.g : 20190301_CEPdata ')
#filename = input('The name of the file is :')

filepath = r'Z:\Laser\CEP\20180615'
filename ='CEP stab_30ms'


#paramaters to be chosen :
Pi_window = 1 # 1 to plot the CEP vs Time between -pi and +pi ; 0 if not
lines2pi = 1  #1 to plot horizontal grey lines every 2Pi, 0 if not




file= str(filename) + '.dat'
path = os.getcwd()
if not os.path.exists('analyzed_data'):
    os.mkdir('analyzed_data')
output_dir=path+'\\analyzed_data'    #output directory = new directory created where this python file is located
files = str(filepath) + "\\" + str(file)

#information contained in the file header :
with open(files) as fp:
    infos = fp.readlines (600)
    print ('\n','\n', '\n', infos[0][0:32], '\n')
    
    param3 = infos[3].split(',')
    print (param3[2])
    print (param3[0])
    print (param3[1], '\n')
    
    
    param4 = infos[4].split(',')
    print(param4[0])
    print(param4[1])
    print(param4[2])    
    print(param4[3], '\n')
    
    param = infos[5].split(',')
    print(param[0])
    print(param[1])
    print(param[2], '\n')
    
    print(param[3])
    print(param[4])
    print(param[5])
    print(param[6], '\n')
    
    print(param[7])
    print(param[8])
    print(param[9], '\n', '\n')
    
    
#load data :    
times = np.loadtxt(files, dtype='str', usecols = (0,),skiprows=8) # probe timing (ns)
phases = np.loadtxt(files, dtype='str', usecols = (1,),skiprows=8) # Phase (radian)
amps = np.loadtxt(files, dtype='str', usecols = (3,),skiprows=8)
voltages =  np.genfromtxt(files, dtype='float', skip_header=8, skip_footer=1, usecols = (6,))


#loop in case there are several measurements :
nb = np.size(times) 
indices = [0]
for i in np.arange(0,nb):
    if 'Meas' in times[i]:
        print ('line', i, ' = end of a measurement')  
        indices.append(i)

for k in np.arange(0,np.size(indices)-1) :
    if k==0:
        b = int(indices[k])
    elif k>0:
        b = int(indices[k])+8
    e = int(indices[k+1]-1)   
    voltage =[] ; time = [] ; phase = [] ; squared_phase = [] ; amp = [] #creation of the lists 
    
    for j in np.arange(b,e):
        time.append(float(times[j][2:7])/1.)     #time is in sec (/1.) or min (/60.)
        phase.append(float(phases[j][2:7]))
        squared_phase.append((float(phases[j][2:7]))**2)
        amp.append(float(amps[j][2:7]))
        voltage.append(float(voltages[j]))
        
        
        
    #Phase deviation rms in mrad
    mean_phase = sum(phase)/np.size(phase)
    mean_squared_phase = sum(squared_phase)/np.size(squared_phase)
    ecart_type = int(round(np.sqrt(mean_squared_phase - mean_phase*mean_phase),3)*1000)
    

    #center values on 0 
    for i in np.arange(0,np.size(phase)):
        phase[i]=phase[i]-mean_phase

       
    #Figure parameters
    fig=plt.figure(figsize=(12,7))
    gs = gridspec.GridSpec(3, 2, wspace = 0,  width_ratios=[2, 1])
    s=12


    
    #CEP as a function of time
    CEPvsT = plt.subplot(gs[1,0])
    plt.plot(time,phase, '.', markersize=1.2,  color='navy')
    plt.ylabel('relative CEP (rad)', size=s, color='navy')
    plt.tick_params(axis='both',labelsize=s)
    plt.tick_params(axis='y',labelcolor ='navy', color = 'navy')
    if lines2pi == 1: # Puts some horizontal lines every 2Pi starting from target value
        target = 0   #CEP target value
        plt.axhline(target, color = 'grey', alpha = 0.8, linewidth = 1)  #horizontal line = target value
        upper = max(round(max(phase)/(2*np.pi)),  abs(round(min(phase)/(2*np.pi))) )
        for n in np.arange(1,upper):  
            plt.axhline(target + n*2*np.pi, color = 'grey', alpha = 0.2, linewidth = 1)
            plt.axhline(target - n*2*np.pi, color = 'grey', alpha = 0.2, linewidth = 1)
    if Pi_window == 1:
        plt.ylim([-np.pi,np.pi])
        plt.yticks([-np.pi, -np.pi/2,0,np.pi/2, np.pi], ['-$\pi$','-$\pi$/2', '0', '$\pi$/2','$\pi$' ], size=s)
    elif Pi_window == 0:
        plt.ylim(min(phase)-0.1*abs(min(phase)), max(phase)+0.1*max(phase) )
    
    
    #Histogram
    histo = plt.subplot(gs[1,1])
    plt.hist(phase, bins=round(ecart_type/20), orientation='horizontal', color='navy')  #bins =.... number of columns
    plt.yticks([])
    plt.axis('off')
    plt.xlabel('occurrence', size=s)
    plt.tick_params(axis='both',labelsize=s)
    if Pi_window == 1:
        plt.ylim([-np.pi,np.pi])
    elif Pi_window == 0:
        plt.ylim(min(phase)-0.1*abs(min(phase)), max(phase)+0.1*max(phase) )
    d,f=np.histogram(phase,bins=round(ecart_type/20))
    plt.text(0.1*max(d), max(phase)-0.1*max(phase) ,'$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)
      
    
    #FFT amplitude = f(time)
    ampli= plt.subplot(gs[0,0])
    plt.plot(time, amp, '.',markersize=2.5, color = 'red', alpha = 0.2)
    plt.ylabel('FFT amplitude', color='red', size=s)
    plt.tick_params(axis='y', labelcolor='red',color='red',labelsize=s)


    # Output voltage = f(time)
    volt = plt.subplot(gs[2,0])
    plt.plot(time, voltage, '.',markersize=2.5, color = 'green', alpha = 0.2)
    plt.xlabel('time(sec)',size=s)
    plt.ylabel('Output voltage (V)', color='green', size=s)
    plt.tick_params(axis='y', labelcolor='green',color='green',labelsize=s)



    # title and save
    fig.suptitle('file : ' + str(files), size = s -4)
    plt.savefig( str(output_dir) + '\\' + filepath[-8::] + '_' + str(filename)+ '_measnum' +str(k) + '_FFTamplitude.png', dpi=300,  bbox_inches='tight')