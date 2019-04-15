# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 12:30:50 2018

@author: ouille
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import os

cwd = os.getcwd()
filepath = r'Z:\Laser\D-Scan\HCF\20180305'

filename0='0mbar-5.00mJ-33880fs2'
filename1='300mbar-4.84mJ-33880fs2'
filename2='600mbar-4.86mJ-33880fs2'
filename3='900mbar-4.96mJ-33880fs2'
filename4='1200mbar-4.87mJ-33880fs2'
filename5='1400mbar-4.65mJ-33880fs2'
filenames = [filename0, filename1, filename2, filename3, filename4, filename5]

E0=5.00 ; E1=4.84 ;E2=4.86 ;E3=4.96 ;E4=4.87 ;E5=4.65
E = [E0,E1,E2,E3,E4,E5]

E0c=4.97 ;E1c=4.92 ;E2c=4.94 ;E3c=4.98 ;E4c=5.00 ;E5c=5.00
Ec=[E0c,E1c,E2c,E3c,E4c,E5c]

P0 = 0 ; P1=300 ;P2 = 600 ; P3=900 ;P4=1200 ;P5=1400
P =[P0,P1,P2,P3,P4,P5]



##############################################

output_dir=filepath
f_spectrum=[] ; f_pulse=[] ; f_retrieved=[] ; f_measured=[] ; spectrum_wl=[]
ret_SHG_wl = [] ; ret_insertions = [] ; ret_dscan = [] ; spectrum_wl=[]; spectrum_I=[] ;spectrum_phase=[]
spectrum_freq=[] ; freq_width=[] ; freq_width2=[] 

for i, value in enumerate(filenames):
    f_spectrum.append(filepath + '\\' + filenames[i] + '_spectrum.csv')
    f_pulse.append(filepath +  '\\' + filenames[i] + '_retrieved_pulse.csv')
    f_retrieved.append( filepath +  '\\' + filenames[i] + '_retrieved_dscan.csv')
    f_measured.append(filepath +  '\\' + filenames[i] + '_measured_dscan.csv')
    
    #Retrieved trace
    data = np.genfromtxt(f_retrieved[i], dtype=float)
    data=data[2:np.size(data)]
    fheader = np.genfromtxt(f_retrieved[i], dtype=None)
    nWL = float(fheader[0][4:7]) #number of wavelengths
    nI= float(fheader[1][3:6])  #number of insertion values
    ret_SHG_wl.append(data[0:nWL]) #list of SHG wavelengths(x-axis)
    ret_insertions.append( data[nWL:(nWL+nI)] ) #list of insertion values(y-axis)
    ret_dscan .append( np.reshape(data[nWL+nI:np.size(data)],(nI, nWL)))  #matrix where each line corresponds to the spectrum for a given insertion

    #Spectrum
    raw_data = np.genfromtxt(f_spectrum[i], delimiter=',', dtype=float)
    spectrum_wl.append(raw_data[:,0])   #list of wavelengths
    spectrum_I.append(raw_data[:,1])  #spectral intensity
    spectrum_phase.append(raw_data[:,2])  #spectral phase
    spectrum_freq.append(2*3.14159265359*299.792458/spectrum_wl[i])

    #boundaries for the spectrum :
    wrmin=[]
    wrmax=[]
    for p, value in enumerate (spectrum_wl[i][0:1400]):
        if spectrum_I[i][p]<=0.0001 :
            wrmin.append(p)
    for p, value in enumerate (np.arange(1400,np.size(spectrum_wl[i])) ): 
        if spectrum_I[i][value]<=0.0001 :
            wrmax.append(value)
    wmin = max(wrmin)
    wmax=min(wrmax)
    
#    #to check the remaining spectrum if needed :
#    plt.figure()
#    plt.plot(spectrum_wl[i][wmin:wmax], spectrum_I[i][wmin:wmax])
    
    #compute the spectral width according to the definition in Pinault's article :
    d_omega = -1 * np.diff(spectrum_freq[i][wmin:wmax])
    spec =abs( spectrum_I[i][wmin+1:wmax])
    specint = np.sum(spec*d_omega)
    omega = spectrum_freq[i][wmin+1:wmax]
    freq_width2.append( np.sum(omega**2 * spec * d_omega)/specint  - (  np.sum(omega * spec  * d_omega)/specint  )**2 )
    freq_width.append(np.sqrt( freq_width2[i]))




########################################

fig = plt.figure(figsize=(11,8))
gs0 = GridSpec(2, 1, height_ratios=[3,2], hspace=0.4, left=0.07, top = 0.98, bottom = 0.08, right = 0.98)

gs00 = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0], wspace = 0.4)
gs01 = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs0[1], wspace= 0.1)


Escan = plt.subplot(gs00[0,0])
Escan.tick_params(axis='both', labelsize=16)
Escan = plt.ylabel('Energy (mJ)', size='17')
Escan = plt.xlabel('Pressure (mbar)',size='17')
Escan = plt.plot(P,E, '--o--', color='navy', label='compressed (0 fs²)')
Escan = plt.plot(P,Ec, '--o--', color= 'r', label='chirped (275 fs²)')
Escan = plt.ylim([0,5.5])
Escan = plt.legend(loc='lower right', fontsize='xx-large')


sigmascan = plt.subplot(gs00[0, 1])
sigmascan.tick_params(axis='both', labelsize=16)
sigmascan=plt.ylabel('Spectral width ($fs^{-1}$)',size='17')
sigmascan=plt.xlabel('Pressure (mbar)',size='17')
for i, value in enumerate(freq_width):
    freq_width [i] = 2*freq_width[i]
sigmascan = plt.ylim([0,1.1])
sigmascan=plt.plot(P, freq_width, '--o--', color='navy', label='experiments')
#add data points from numerical simulations :
simulation = np.loadtxt(r'C:\Users\ouille\Desktop\Article_LASER\Figure4\\pressuredependence.dat', skiprows=(1))
X = simulation[:,0]
Y= simulation [:,2]
sigmascan=plt.plot(X*1000,Y, '--o--', color='r', label='simulation')
sigmascan=plt.legend(loc='upper left', fontsize='xx-large')



###################################################################################

## if an extra trace is wanted : 
#trace0 = plt.subplot(gs01[0,0])
#trace0.set_ylabel('Insertion (mm)',size='15')
#trace0.set_xlabel('Lambda (nm)',size='15')
##trace0.yaxis.set_visible(False)
#trace0.text(280,3.4,'vacuum', color='white', size='12')
#trace0 = plt.pcolor(ret_SHG_wl[0], ret_insertions[0], ret_dscan[0], cmap='jet')

trace1 = plt.subplot(gs01[0,0])
trace1.set_ylabel('Insertion (mm)',size='17')
trace1.set_xlabel('$\lambda$ (nm)',size='17')
trace1.text(280,3.4,'0.3 bar', color='white', size='15', weight='bold')
trace1.tick_params(axis='both', labelsize=16)
trace1 = plt.pcolor(ret_SHG_wl[1], ret_insertions[1], ret_dscan[1], cmap='jet')


trace2 = plt.subplot(gs01[0,1])
trace2.yaxis.set_visible(False)
trace2.set_xlabel('$\lambda$ (nm)',size='17')
trace2.tick_params(axis='both', labelsize=16)
trace2.text(280,3.4,'0.6 bar', color='white', size='15', weight='bold')
trace2 = plt.pcolor(ret_SHG_wl[2], ret_insertions[2], ret_dscan[2], cmap='jet')

trace3 = plt.subplot(gs01[0,2])
trace3.set_xlabel('$\lambda$ (nm)',size='17')
trace3.yaxis.set_visible(False)
trace3.tick_params(axis='both', labelsize=16)
trace3.text(280,3.4,'0.9 bar', color='white', size='15', weight='bold')
trace3 = plt.pcolor(ret_SHG_wl[3], ret_insertions[3], ret_dscan[3], cmap='jet')

trace4 = plt.subplot(gs01[0,3])
trace4.set_xlabel('$\lambda$ (nm)',size='17')
trace4.yaxis.set_visible(False)
trace4.tick_params(axis='both', labelsize=16)
trace4.text(280,3.4,'1.2 bar', color='white', size='15', weight='bold')
trace4 = plt.pcolor(ret_SHG_wl[4], ret_insertions[4], ret_dscan[4], cmap='jet')

## if an extra trace is wanted : 
#trace5 = plt.subplot(gs01[0,4])
#trace5.set_xlabel('Lambda (nm)',size='15')
#trace5.yaxis.set_visible(False)
#trace5.text(280,3.4,'P=1400 mbar', color='white', size='12')
#trace5 = plt.pcolor(ret_SHG_wl[5], ret_insertions[5], ret_dscan[5], cmap='jet')

plt.savefig(str(cwd) + '\\' + 'Figure4'+ '.png')
