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




        
def APS (filepath, filename, n, save=1, picenter=0, lines2pi=0):
       
    
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
    print ('\n', 'This file contains ', str(np.size(debuts)), ' measurements', '\n', '\n')

#    info(n)
#    APSplot(n)





#def info (n): #information contained in the file header :
    subheader = line[debuts[n]:debuts[n]+9]
    print ('\n','\n', '\n','Measurement n°',str(n+1), '/', str(np.size(debuts)),'\n','Here are some basic information on your ', subheader[0][0:32], ': \n')
    param3 = subheader[3].split(',')
    print (param3[2])
    print (param3[0])
    print (param3[1], '\n')
    param4 =subheader[4].split(',')
    print(param4[0])
    print(param4[1])
    print(param4[2])    
    print(param4[3], '\n')
    param = subheader[5].split(',')
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
    







#def APSplot (n):

    subdata = line[debuts[n]+9 :fins[n]]
    time = [];phase=[];amp=[];voltage=[]; squared_phase=[]
    for k, p in enumerate(subdata):
        time.append(float(p.split('\t')[0]))
        phase.append(float(p.split('\t')[1]))
        amp.append(float(p.split('\t')[3]))
        voltage.append(float(p.split('\t')[6]))
        squared_phase.append(phase[k]*phase[k])

        
        
             
        
        
        
    if partial == 0 :    
        #Phase deviation rms in mrad
        mean_phase = sum(phase)/np.size(phase)
        mean_squared_phase = sum(squared_phase)/np.size(squared_phase)
        ecart_type = int(round(np.sqrt(mean_squared_phase - mean_phase*mean_phase),3)*1000)
        
        #center values on 0 
        for i in np.arange(0,np.size(phase)):
            phase[i]=phase[i]-mean_phase
           
        #Figure parameters
        fig=plt.figure(figsize=(12,7))
        gs = gridspec.GridSpec(3, 2, wspace = 0,  width_ratios=[3, 1])
        s=12
            
        #CEP as a function of time
        plt.subplot(gs[1,0])
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
        if picenter == 1:
            plt.ylim([-np.pi,np.pi])
            plt.yticks([-np.pi, -np.pi/2,0,np.pi/2, np.pi], ['-$\pi$','-$\pi$/2', '0', '$\pi$/2','$\pi$' ], size=s)
        elif picenter == 0:
            plt.ylim(min(phase)-0.1*abs(min(phase)), max(phase)+0.1*max(phase) )
        
        #Histogram
        plt.subplot(gs[1,1])
        plt.hist(phase, bins=round(ecart_type/20), orientation='horizontal', color='navy')  #bins =.... number of columns
        plt.yticks([])
        plt.axis('off')
        plt.xlabel('occurrence', size=s)
        plt.tick_params(axis='both',labelsize=s)
        if picenter== 1:
            plt.ylim([-np.pi,np.pi])
        elif picenter == 0:
            plt.ylim(min(phase)-0.1*abs(min(phase)), max(phase)+0.1*max(phase) )
        d,f=np.histogram(phase,bins=round(ecart_type/20))
        plt.title('$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)
        
        #FFT amplitude = f(time)
        plt.subplot(gs[0,0])
        plt.plot(time, amp, '.',markersize=1.2, color = 'red')
        plt.ylabel('FFT amplitude', color='red', size=s)
        plt.tick_params(axis='y', labelcolor='red',color='red')
        plt.tick_params(axis='both',labelsize=s)    
        
        # Output voltage = f(time)
        plt.subplot(gs[2,0])
        plt.plot(time, voltage, '.',markersize=1.2, color = 'green')
        plt.xlabel('time(sec)',size=s)
        plt.ylabel('Output voltage (V)', color='green', size=s)
        plt.tick_params(axis='y', labelcolor='green',color='green')
        plt.tick_params(axis='both',labelsize=s)
        
        # title and save
        fig.suptitle('file : ' + str(files) + '\n' + 'Measurement n°' + str(n+1) + '/' +str(np.size(debuts)), size = s -4)
        if save == 1:
            plt.savefig( str(output_dir) + '\\' + str(filename)+ '_measnum' +str(n) + '.png', dpi=900,  bbox_inches='tight')
            print ('A .png file has been saved in the folder "analyzed_data" with 3 plots : FFT amplitude, CEP and Output voltage versus time','\n','\n','You will find the "analyzed_data" folder in your computer at the same lacoation as this .py file you are executing right now!' )
            
    
  
            
            
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
        

    #FOR PARTIAL PLOTS
    elif partial == 1:
        
        def index(l, p): return len([x for x in takewhile(lambda x: x[1] <= p, enumerate(l))])
        a = index(time, start)
        b = index(time, end)
    
        
            
        #Phase deviation rms in mrad
        mean_phase = sum(phase[a:b])/np.size(phase[a:b])
        mean_squared_phase = sum(squared_phase[a:b])/np.size(squared_phase[a:b])
        ecart_type = int(round(np.sqrt(mean_squared_phase - mean_phase*mean_phase),3)*1000)
        
        #center values on 0 
        for i in np.arange(0,np.size(phase)):
            phase[i]=phase[i]-mean_phase
           
        #Figure parameters
        fig=plt.figure(figsize=(12,7))
        gs = gridspec.GridSpec(3, 2, wspace = 0,  width_ratios=[3, 1])
        s=12
            
        #CEP as a function of time
        plt.subplot(gs[1,0])
        plt.plot(time[a:b],phase[a:b], '.', markersize=1.2,  color='navy')
        plt.ylabel('relative CEP (rad)', size=s, color='navy')
        plt.tick_params(axis='both',labelsize=s)
        plt.tick_params(axis='y',labelcolor ='navy', color = 'navy')
        if lines2pi == 1: # Puts some horizontal lines every 2Pi starting from target value
            target = 0   #CEP target value
            plt.axhline(target, color = 'grey', alpha = 0.8, linewidth = 1)  #horizontal line = target value
            upper = max(round(max(phase[a:b])/(2*np.pi)),  abs(round(min(phase[a:b])/(2*np.pi))) )
            for n in np.arange(1,upper):  
                plt.axhline(target + n*2*np.pi, color = 'grey', alpha = 0.2, linewidth = 1)
                plt.axhline(target - n*2*np.pi, color = 'grey', alpha = 0.2, linewidth = 1)
        if picenter == 1:
            plt.ylim([-np.pi,np.pi])
            plt.yticks([-np.pi, -np.pi/2,0,np.pi/2, np.pi], ['-$\pi$','-$\pi$/2', '0', '$\pi$/2','$\pi$' ], size=s)
        elif picenter == 0:
            plt.ylim(min(phase[a:b])-0.1*abs(min(phase[a:b])), max(phase[a:b])+0.1*max(phase[a:b]) )
        
        #Histogram
        plt.subplot(gs[1,1])
        plt.hist(phase[a:b], bins=round(ecart_type/20), orientation='horizontal', color='navy')  #bins =.... number of columns
        plt.yticks([])
        plt.axis('off')
        plt.xlabel('occurrence', size=s)
        plt.tick_params(axis='both',labelsize=s)
        if picenter== 1:
            plt.ylim([-np.pi,np.pi])
        elif picenter == 0:
            plt.ylim(min(phase[a:b])-0.1*abs(min(phase[a:b])), max(phase[a:b])+0.1*max(phase[a:b]) )
        d,f=np.histogram(phase,bins=round(ecart_type/20))
        plt.title('$\sigma_{rms}$ = ' + str(ecart_type) + ' mrad', size=s+1)
        
        #FFT amplitude = f(time)
        plt.subplot(gs[0,0])
        plt.plot(time[a:b], amp[a:b], '.',markersize=1.2, color = 'red')
        plt.ylabel('FFT amplitude', color='red', size=s)
        plt.tick_params(axis='y', labelcolor='red',color='red')
        plt.tick_params(axis='both',labelsize=s)    
        
        # Output voltage = f(time)
        plt.subplot(gs[2,0])
        plt.plot(time[a:b], voltage[a:b], '.',markersize=1.2, color = 'green')
        plt.xlabel('time(sec)',size=s)
        plt.ylabel('Output voltage (V)', color='green', size=s)
        plt.tick_params(axis='y', labelcolor='green',color='green')
        plt.tick_params(axis='both',labelsize=s)
        
        # title and save
        fig.suptitle('file : ' + str(files) + '\n' + 'Measurement n°' + str(n+1) + '/' +str(np.size(debuts)), size = s -4)
        if save == 1:
            plt.savefig( str(output_dir) + '\\' + str(filename)+ '_measnum' +str(n) + '_interval_' + str(start) +'_'+str(end) + '.png', dpi=900,  bbox_inches='tight')
            print ('A .png file has been saved in the folder "analyzed_data" with 3 plots : FFT amplitude, CEP and Output voltage versus time','\n','\n','You will find the "analyzed_data" folder in your computer at the same lacoation as this .py file you are executing right now!' )
            
    
    
    
    
    
    
    
  
    
    
    
    
    
    

#PARTIAL PLOT PARAMATERS 
partial = 0 #1 if YES, you want to make a partial plot   
#0 if NO you don't want to make a partial plot
start = 0  #begining of the plot in seconds
end = 10 #end of the plot in seconds 



APS (filepath, filename, 0)   #the last argument is the measurement number ; Be careful, it starts with 0 !!!
