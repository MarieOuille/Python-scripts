# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 17:38:52 2017

@author: ouille
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import os

#Script pour traiter les données de l'APS800

# FIRST : EXTRACT DATA 
filepath = r'Z:\Laser\CEP\20181010'
filename ='power_amp_no_pump_PC2at7-4kV'
file= str(filename) + '.dat'

#output directory = new directory created where this python file is located
path = os.getcwd()
if not os.path.exists('analyzed_data'):
    os.mkdir('analyzed_data')
output_dir=path+'\\analyzed_data'

files = str(filepath) + "\\" + str(file)
times = np.loadtxt(files, dtype='str', usecols = (0,),skiprows=9) # probe timing (ns)
phases = np.loadtxt(files, dtype='str', usecols = (1,),skiprows=9) # Phase (radian)
amps = np.loadtxt(files, dtype='str', usecols = (3,),skiprows=9)


nb = np.size(times) 
indices = [-8]
for i in np.arange(0,nb):
    if 'Meas' in times[i]:
        print ('line', i, ' = end of a measurement')  
        indices.append(i)

for k in np.arange(0,np.size(indices)-1) :
    b = int(indices[k]+8)
    e = int(indices[k+1]-1) 
    print(b,e)    
    time = [] ; phase = [] ; squared_phase = [] ; amp = [] #creation of the lists 
    
    for j in np.arange(b,e):
        time.append(float(times[j][2:7])/1.)     #time is in sec (/1.) or min (/60.)
        phase.append(float(phases[j][2:7]))
        squared_phase.append((float(phases[j][2:7]))**2)
        amp.append(float(amps[j][2:7]))

        
        

        
    #Phase deviation rms in mrad
    mean_phase = sum(phase)/np.size(phase)
    mean_squared_phase = sum(squared_phase)/np.size(squared_phase)
    ecart_type = int(round(np.sqrt(mean_squared_phase - mean_phase*mean_phase),3)*1000)


    #center values on 0 (optional)
    for i in np.arange(0,np.size(phase)):
        phase[i]=phase[i]-mean_phase

       
    #Figure parameters
    fig=plt.figure(figsize=(10,3))
    gs = gridspec.GridSpec(1, 2, wspace = 0,  width_ratios=[2, 1])
    s=12

    
    #CEP as a function of time
    CEPvsT = plt.subplot(gs[0,0])
    plt.plot(time,phase, '.', markersize=1.2,  color='navy')
    plt.xlabel('time (sec)', size=s)
    plt.ylabel('relative CEP (rad)', size=s)
    plt.tick_params(axis='both',labelsize=s)
    target = 0
    #plt.axhline(target, color = 'grey', alpha = 0.8, linewidth = 1)
    #plt.axhline(1.6, color = 'grey', alpha = 0.8, linewidth = 1)
    #for n in np.arange(2,11):
    #    plt.axhline(target + n*2*np.pi, color = 'grey', alpha = 0.8, linewidth = 1)
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
    plt.text(2,np.mean(phase)+1,'$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)
    
    
    #Save figure
    fig.text(0.1,1.05,'file : ' + str(files), size = s -4)
    plt.savefig( str(output_dir) + '\\' + str(filename)+ '_measnum' +str(k) + '_phase_histo.png', dpi=300,  bbox_inches='tight')
   
    



##optional : amplitude and CEP vs time
#plt.figure(figsize=(7,3))
#plt.plot(time, amp, '.',markersize=2.5, color = 'red', alpha = 0.2)
#plt.ylim(np.min(amp)-1, np.max(amp)+0.1)
#plt.xlabel('time(sec)',size=s)
#plt.ylabel('FFT amplitude', color='red', size=s)
#plt.tick_params(axis='y', labelcolor='red',color='red',labelsize=s)
#B = plt.twinx()
#B.plot(time, phase, '.',markersize=2.5, color = 'navy')
#B.set_ylim(-np.pi/2, np.pi)
#B.set_ylabel('CEP (rad)', color = 'navy', size=s)
#B.tick_params(axis='y', labelcolor='navy',color='navy',labelsize=s)






##optional bis : plot just one part of the graph 
#interval selection :
start=10
end = 30
from itertools import takewhile
def st(l, b): return len([x for x in takewhile(lambda x: x[1] <= b, enumerate(l))])
a = st(time,start)
b = st(time, end)


#Phase deviation rms in mrad
mean_phase2 = sum(phase[a:b])/np.size(phase[a:b])

for i in np.arange(0,np.size(phase)-1):
    squared_phase[i]=phase[i]*phase[i]
         
mean_squared_phase2 = sum(squared_phase[a:b])/np.size(squared_phase[a:b])
ecart_type = int(round(np.sqrt(mean_squared_phase2 - mean_phase2*mean_phase2),3)*1000)

#center values on 0
for i in np.arange(0,np.size(phase)-1):
    phase[i]=phase[i]-mean_phase

       
#Figure parameters
fig=plt.figure(figsize=(10,3))
gs = gridspec.GridSpec(1, 2, wspace = 0,  width_ratios=[2, 1])
s=12

#Histogram
histo = plt.subplot(gs[0,1])
plt.hist(phase[a:b], bins=round(ecart_type/20), orientation='horizontal', color='navy')  #bins =.... number of columns
plt.yticks([])
plt.xlabel('occurrence', size=s)
plt.tick_params(axis='both',labelsize=s)
#plt.ylim([-np.pi,np.pi])
plt.title('$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)

#CEP as a function of time
CEPvsT = plt.subplot(gs[0,0])
plt.plot(time[a:b],phase[a:b], '.', markersize=1.2,  color='navy')
plt.xlabel('time (sec)', size=s)
plt.ylabel('relative CEP (rad)', size=s)
plt.tick_params(axis='both',labelsize=s)
#plt.ylim([-np.pi,np.pi])
#plt.yticks([-np.pi, -np.pi/2,0,np.pi/2, np.pi], ['-$\pi$','-$\pi$/2', '0', '$\pi$/2','$\pi$' ], size=s)

#Save figure
fig.text(0.1,1.05,'file : ' + str(files), size = s -4)
plt.savefig( str(output_dir) + '\\' + str(filename) + '_phase_histo_' + str(start) + '_' + str(end) +  '.png', dpi=300,  bbox_inches='tight')