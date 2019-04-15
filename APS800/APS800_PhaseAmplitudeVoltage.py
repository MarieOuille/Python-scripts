# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 17:38:52 2017

@author: ouille
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import os
from itertools import takewhile

## FIRST : EXTRACT DATA 
#print ('Where is the file located? e.g : Z:\Laser\CEP\20190311')
#filepath = input('The path of the file is :')
#print ('What is the name of the .dat file? e.g : 20190301_CEPdata ')
#filename = input('The name of the file is :')

filepath = r'Z:\Laser\CEP\20180619'
filename ='HCF-4fs-1500mbar-ter-single'

file= str(filename) + '.dat'

#output directory = new directory created where this python file is located
path = os.getcwd()
if not os.path.exists('analyzed_data'):
    os.mkdir('analyzed_data')
output_dir=path+'\\analyzed_data'

files = str(filepath) + "\\" + str(file)
times = np.loadtxt(files, dtype='str', usecols = (0,),skiprows=8) # probe timing (ns)
phases = np.loadtxt(files, dtype='str', usecols = (1,),skiprows=8) # Phase (radian)
amps = np.loadtxt(files, dtype='str', usecols = (3,),skiprows=8)
voltages =  np.genfromtxt(files, dtype='str', skip_header=8, skip_footer=1, usecols = (6,))

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
    print(b,e)    
    voltage =[] ; time = [] ; phase = [] ; squared_phase = [] ; amp = [] #creation of the lists 
    
    for j in np.arange(b,e):
        time.append(float(times[j][2:7])/1.)     #time is in sec (/1.) or min (/60.)
        phase.append(float(phases[j][2:7]))
        squared_phase.append((float(phases[j][2:7]))**2)
        amp.append(float(amps[j][2:7]))
        voltage.append(float(voltages[j][2:7]))
        
        

        
    #Phase deviation rms in mrad
    mean_phase = sum(phase)/np.size(phase)
    mean_squared_phase = sum(squared_phase)/np.size(squared_phase)
    ecart_type = int(round(np.sqrt(mean_squared_phase - mean_phase*mean_phase),3)*1000)


    #center values on 0 (optional)
    for i in np.arange(0,np.size(phase)):
        phase[i]=phase[i]-mean_phase

       
    #Figure parameters
    fig=plt.figure(figsize=(10,6))
    gs = gridspec.GridSpec(3, 2, wspace = 0,  width_ratios=[2, 1])
    s=12

    
    #CEP as a function of time
    CEPvsT = plt.subplot(gs[0,0])
    plt.plot(time,phase, '.', markersize=1.2,  color='navy')
    plt.xlabel('time (sec)', size=s)
    plt.ylabel('relative CEP (rad)', size=s)
    plt.tick_params(axis='both',labelsize=s)
    
    ## Put some horizontal lines every 2Pi starting from target value
    target = 0
    plt.axhline(target, color = 'grey', alpha = 0.8, linewidth = 1)
    for n in np.arange(1,11):
        plt.axhline(target + n*2*np.pi, color = 'grey', alpha = 0.2, linewidth = 1)
        plt.axhline(target - n*2*np.pi, color = 'grey', alpha = 0.2, linewidth = 1)
    ##remove the "##" of the 2 following line to center around 0 and display from -Pi to +Pi :
    #plt.ylim([-np.pi,np.pi])
    #plt.yticks([-np.pi, -np.pi/2,0,np.pi/2, np.pi], ['-$\pi$','-$\pi$/2', '0', '$\pi$/2','$\pi$' ], size=s)
        
    
    #Histogram
    histo = plt.subplot(gs[0,1])
    plt.hist(phase, bins=round(ecart_type/20), orientation='horizontal', color='navy')  #bins =.... number of columns
    plt.yticks([])
    plt.axis('off')
    plt.xlabel('occurrence', size=s)
    plt.tick_params(axis='both',labelsize=s)
    #plt.ylim([-np.pi,np.pi])
    d,f=np.histogram(phase,bins=round(ecart_type/20))
    plt.text(0.1*max(d), max(phase)-0.1*max(phase) ,'$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)
    
    









#    ##optional bis : plot just one part of the graph 
#    #interval selection :
#    start = 500 #s
#    end = 600 #s
#    def st(l, b): return len([x for x in takewhile(lambda x: x[1] <= b, enumerate(l))])
#    a = st(time,start)
#    b = st(time, end)
#    
#    
#    #Phase deviation rms in mrad
#    mean_phase2 = sum(phase[a:b])/np.size(phase[a:b])
#    
#    for i in np.arange(0,np.size(phase)):
#        squared_phase[i]=phase[i]*phase[i]
#             
#    mean_squared_phase2 = sum(squared_phase[a:b])/np.size(squared_phase[a:b])
#    ecart_type = int(round(np.sqrt(mean_squared_phase2 - mean_phase2*mean_phase2),3)*1000)
#    
#    
#           
#    #Figure parameters
#    fig=plt.figure(figsize=(10,3))
#    gs = gridspec.GridSpec(1, 2, wspace = 0,  width_ratios=[2, 1])
#    s=12
#    
#    #Histogram
#    histo = plt.subplot(gs[0,1])
#    plt.hist(phase[a:b], bins=round(ecart_type/20), orientation='horizontal', color='navy')  #bins =.... number of columns
#    plt.yticks([])
#    plt.xlabel('occurrence', size=s)
#    #plt.ylim([-np.pi,np.pi])
#    plt.tick_params(axis='both',labelsize=s)
#    d,f=np.histogram(phase[a:b], bins=round(ecart_type/20))
#    plt.axis('off')
#    plt.text(0.1*max(d), max(phase[a:b])-0.1*abs(max(phase[a:b])-np.mean(phase[a:b])), '$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)
#    
#    #CEP as a function of time
#    CEPvsT = plt.subplot(gs[0,0])
#    plt.plot(time[a:b],phase[a:b], '.', markersize=1.2,  color='navy')
#    plt.xlabel('time (sec)', size=s)
#    plt.ylabel('relative CEP (rad)', size=s)
#    plt.tick_params(axis='both',labelsize=s)
#    #plt.ylim([-np.pi,np.pi])
#    #plt.yticks([-np.pi, -np.pi/2,0,np.pi/2, np.pi], ['-$\pi$','-$\pi$/2', '0', '$\pi$/2','$\pi$' ], size=s)
#    
#    #Save figure
#    fig.text(0.1,1.05,'file : ' + str(files), size = s -4)
#    plt.savefig( str(output_dir) + '\\' + str(filename) + '_phase_histo_' + str(start) + '_' + str(end) +  '.png', dpi=300,  bbox_inches='tight')
#




    
    
    
    #FFT amplitude = f(time)
    ampli= plt.subplot(gs[1,0])
    plt.plot(time, amp, '.',markersize=2.5, color = 'red', alpha = 0.2)
    plt.xlabel('time(sec)',size=s)
    plt.ylabel('FFT amplitude', color='red', size=s)
    plt.tick_params(axis='y', labelcolor='red',color='red',labelsize=s)


    # Output voltage = f(time)
    volt = plt.subplot(gs[2,0])
    plt.plot(time, voltage, '.',markersize=2.5, color = 'green', alpha = 0.2)
    plt.xlabel('time(sec)',size=s)
    plt.ylabel('Output voltage', color='green', size=s)
    plt.tick_params(axis='y', labelcolor='green',color='green',labelsize=s)
#    plt.ylim(-10,10)


    # title and save
    fig.suptitle('file : ' + str(files), size = s -4)
    plt.savefig( str(output_dir) + '\\' + str(filename)+ '_measnum' +str(k) + '_FFTamplitude.png', dpi=300,  bbox_inches='tight')