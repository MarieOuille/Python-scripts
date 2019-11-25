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
#print ('\n','Where is the file located?  (e.g : Z:\Laser\CEP\20190311)')
#filepath = input('The path of the file is :  ')
#print ('\n', 'What is the name of the .dat file?  (e.g : 20190301_CEPdata) ')
#filename = input('The name of the file is :  ')

filepath = r'Z:\Laser\CEP\20190430'
filename ='30ms'




        
def APS (filepath, filename, n, save=1, picenter=0, lines2pi=0, partial=0, start=0, end=10):
           
    path = os.getcwd()
    if not os.path.exists('analyzed_data'):
        os.mkdir('analyzed_data')
    output_dir=path+'\\analyzed_data'    #output directory = new directory created where this python file is located
    files = str(filepath) + "\\" + str(filename) + '.dat'                  
               
    #split the .dat file into different measurements
    fins = [] ; debuts = []
    data = open(files, 'r')
    line = data.readlines()
            

    for i, value in enumerate(line):
        if 'Measurement finished' in line[i]:
            fins.append(i)
        if 'APS800 log file' in line[i]:
            debuts.append(i)
    print ('\n', 'This file contains ', str(np.size(debuts)), ' measurement(s)', '\n', '\n')


    

    #subdata = data of measurement number n only
    subdata = line[debuts[n]+9 :fins[n]]
    time = [];phase=[];amp=[];voltage=[]; squared_phase=[] ; target = []; squared_amp=[]
    for k, p in enumerate(subdata):
        time.append(float(p.split('\t')[0]))
        phase.append(float(p.split('\t')[1]))
        amp.append(float(p.split('\t')[3]))
        voltage.append(float(p.split('\t')[6]))
        target.append(float(p.split('\t')[2]))
        squared_phase.append(phase[k]*phase[k])
        squared_amp.append(amp[k]*amp[k])

        
        
        
        
        
        
        
        
        
        
        
        

    def index(l, p): return len([x for x in takewhile(lambda x: x[1] <= p, enumerate(l))])
    
    if partial == 0:
        a = index(time, time[0])-1
        b = index(time, time[-1])

    elif partial == 1:  
        a = index(time, start)
        b = index(time, end)

    
    #Phase deviation rms in mrad
    mean_phase = sum(phase[a:b])/np.size(phase[a:b])
    mean_squared_phase = sum(squared_phase[a:b])/np.size(squared_phase[a:b])
    ecart_type = int(round(np.sqrt(mean_squared_phase - mean_phase*mean_phase),3)*1000)
    #FFT amplitude deviation
    mean_amp = sum(amp[a:b])/np.size(amp[a:b])          #amplitude deviation rms in mrad
    mean_squared_amp = sum(squared_amp[a:b])/np.size(squared_amp[a:b])           #amplitude deviation rms in mrad
    ecart_type_amp = (round(np.sqrt(mean_squared_amp - mean_amp*mean_amp),3))    
    
    #center values on 0 
    for i in np.arange(0,np.size(phase)):
        phase[i]=phase[i]-mean_phase
        target[i]=target[i] - mean_phase
       
        
              
              
              
              
              
              
              
              
              
              
              
              
              
              
    #Figure parameters
    fig=plt.figure(figsize=(12,9))
    gs = gridspec.GridSpec(3, 2, wspace = 0,  width_ratios=[3, 1], top = 0.75)
    s=12
    
         
    #CEP as a function of time
    plt.subplot(gs[1,0])
    plt.plot(time[a:b],phase[a:b], '.', markersize=1,  color='navy', label='measured phase')
    plt.plot(time[a:b],target[a:b], '--', markersize=0.5,  color='grey', alpha =0.8, label = 'target phase')
    plt.ylabel('relative CEP (rad)', size=s, color='navy')
    plt.tick_params(axis='both',labelsize=s)
    plt.tick_params(axis='y',labelcolor ='navy', color = 'navy')
    if lines2pi == 1: # Puts some horizontal lines every 2Pi starting from target value
        targetp = target[a]   #CEP target value
        plt.axhline(targetp, color = 'grey', alpha = 0.8, linewidth = 1)  #horizontal line = target value
        upper = max(round(max(phase[a:b])/(2*np.pi)),  abs(round(min(phase[a:b])/(2*np.pi))) )
        for k in np.arange(1,upper+1):  
            plt.axhline(targetp + k*2*np.pi, color = 'grey', alpha = 0.3, linewidth = 1, label='$\pm 2n\pi$, n integer')
            plt.axhline(targetp - k*2*np.pi, color = 'grey', alpha = 0.3, linewidth = 1)
    if picenter == 1:
        plt.ylim([-np.pi,np.pi])
        plt.yticks([-np.pi, -np.pi/2,0,np.pi/2, np.pi], ['-$\pi$','-$\pi$/2', '0', '$\pi$/2','$\pi$' ], size=s)
    elif picenter == 0:
        plt.ylim(min(phase[a:b])-0.1*abs(min(phase[a:b])), max(phase[a:b])+0.1*max(phase[a:b]) )
    plt.legend(loc= 'best', fontsize=s-4)
      
    #Histogram CEP
    plt.subplot(gs[1,1])
    plt.hist(phase[a:b], bins=15, orientation='horizontal', color='navy')  #bins =.... number of columns
    plt.yticks([])
    plt.axis('off')
    plt.xlabel('occurrence', size=s)
    plt.tick_params(axis='both',labelsize=s)
    if picenter== 1:
        plt.ylim([-np.pi,np.pi])
    elif picenter == 0:
        plt.ylim(min(phase[a:b])-0.1*abs(min(phase[a:b])), max(phase[a:b])+0.1*max(phase[a:b]) )
    plt.title('$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)
         
    #FFT amplitude = f(time)
    plt.subplot(gs[0,0])
    plt.plot(time[a:b], amp[a:b], '.',markersize=1, color = 'red')
    plt.ylabel('FFT amplitude', color='red', size=s)
    plt.tick_params(axis='y', labelcolor='red',color='red')
    plt.tick_params(axis='both',labelsize=s)         
    
    #Historgram Amplitude
    plt.subplot(gs[0,1])
    plt.hist(amp, bins=12, orientation='horizontal', color='red')  #bins =.... number of columns
    plt.axis('off')
    plt.title('$\sigma_{rms}$ = ' + str(round(ecart_type_amp/mean_amp*100,1)) + '% = ' + str(round(ecart_type_amp,2)), size=s+1)
                        
    # Output voltage = f(time)
    plt.subplot(gs[2,0])
    plt.plot(time[a:b], voltage[a:b], '.',markersize=1, color = 'green')
    plt.xlabel('time(sec)',size=s)
    plt.ylabel('Output voltage (V)', color='green', size=s)
    plt.tick_params(axis='y', labelcolor='green',color='green')
    plt.tick_params(axis='both',labelsize=s)
    
    # title and save
    subheader = line[debuts[n]:debuts[n]+9]     #information contained in the file header :
    param3 = subheader[3].split(',')
    param4 =subheader[4].split(',')
    param = subheader[5].split(',')
    fig.suptitle('file : ' + str(files) + '\n' + 'Measurement n°' + str(n+1) + '/' +str(np.size(debuts)) + '\n' + ' \n' + str(param3[2]) + ', ' + str(param3[0]) + ', ' + str(param3[1]) + '\n' + str(param4[0]) + ', ' + str(param4[1]) + ', ' + str(param4[2])+ ', ' + str(param4[3]) + '\n' + str(param[0]) + ', ' + str(param[1]) + ', ' + str(param[2]) + '\n' + str(param[3]) + ', ' +str(param[4]) + ', ' +str(param[5]) + ', ' + str(param[6]) + '\n' + str(param[7]) +', ' + str(param[8]) +' ,' + str(param[9]), size = s-1)
    
    
    if save == 1:
        if partial == 1:
            plt.savefig( str(output_dir) + '\\' + str(filename)+ '_measnum' +str(n) + '_interval_' + str(start) +'_'+str(end) + '.png', dpi=900,  bbox_inches='tight')
        elif partial == 0 :
            plt.savefig( str(output_dir) + '\\' + str(filename)+ '_measnum' +str(n) + '.png', dpi=900,  bbox_inches='tight')
        print ('A .png file has been saved in the folder "analyzed_data" with 3 plots : FFT amplitude, CEP and Output voltage versus time','\n','\n','You will find the "analyzed_data" folder in your computer at the same lacoation as this .py file you are executing right now!' )
        


    

    
    



APS (filepath, filename, 0, save=1, picenter=0, lines2pi=0, partial = 0, start=0 , end=10)   
#the 3rd argument is the measurement number you want to plot ; Be careful, it starts with 0 !!!
#save = 1 to save the .png file
#picenter = 1 : YES, plot between -pi and +pi
# lines2pi = 1 : ADDS some horizontal lines every 2pi
# partial = 1 : YES, you want to make a partial plot   
# start and end : beginning and end of the partial plot (seconds)