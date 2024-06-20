# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 09:47:32 2022

@author: marie
"""

import numpy as np
import scipy
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import interpolate
from scipy import optimize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
s=12


 
def lam2freq(x):   #lambda in µm and Freq in THz
    return c/x*1000


def freq2lam(x):
    return c/x*1000



#%% bbo shg  1030 (o) + 1030 (o) = 515 (e)
## WORK IN PROGRESS, TO BE CONFIRMED

# Ref: Kato 1986
def noBBO(l):
    return np.sqrt(2.7359 + 0.01878/(l**2 - 0.01822) - 0.01354*l**2  )

def neBBO(l):
    return np.sqrt(2.3753 + 0.01224/(l**2 - 0.01667) - 0.01516*l**2  )

def netBBO(l,tet): 
    return np.sqrt(1/(np.cos(tet)**2/noBBO(l)**2 + np.sin(tet)**2/neBBO(l)**2 ))


#%% 
## WORK IN PROGRESS, TO BE CONFIRMED

Angle_FWHM=[]
Duration_FWHM=[]

C_thicks = [500, 1000 , 2000]
for it, C_thick in enumerate (C_thicks):
    
    Angle = np.arange(15,30,0.01)
    c=0.299792458   #µm/fs
    freq_min=50#freq2lam(3.4+1)
    freq_max=400#freq2lam(3.4-1)
    THz_freq=np.arange(freq_min,freq_max,0.04)      #
    Lpump=1.03 #µm   
    LTHz = 1/(THz_freq/c)*1000
    
    ## frequency (THz) range for the FWHM calculation:
    r1 =250
    r2= 350
    
    # angle (deg) range for the FWHM calculation:
    ra1= 20
    ra2= 28    
    
    def FWHM(x,y,a=freq_min, b=freq_max):
        ia = np.where(np.isclose(x,a))[0][0] 
        ib = np.where(np.isclose(x,b))[0][0] 
        xnew= x [ia : ib]
        ynew= y [ia:ib]
        ymax = np.max(ynew)
    #print('\n', ymax*0.5 , 'max*0.5')
        for i, value in enumerate (ynew):
            if value==ymax:
                index_max=i
                print(index_max)
                
        for i, value in enumerate ( ynew[0:index_max] ):
            if value >= 0.5*ynew[index_max]:
                low = i
                break
        print(low)
        print(ynew[low] , 'y[low]')
        
        for i, value in enumerate (ynew):
            if i > index_max+1:
                if value <= 0.5*ynew[index_max]:
                    high = i
                    break
        print(high)
        print(ynew[high] , 'y[high]', '\n')
        
        return (xnew[high] - xnew[low])
    
    
    
    
    #%Type 1 phase matching
    print('o + o --> e')
    
    Pump_freq=c/Lpump*1000
    idler_freq=Pump_freq+THz_freq
    Lidler=1/(idler_freq/c)*1000
    PM_BBO = np.zeros((np.size(Angle), np.size(LTHz)))
    
    for i, b in enumerate (Angle):
          CA = b*np.pi/180
          Dk = 2*np.pi*noBBO(Lpump)/Lpump - 2*np.pi*netBBO(Lidler,CA)/Lidler + 2*np.pi*noBBO(LTHz)/LTHz
          #Dk = 2*np.pi*net_GS(Lpump,CA)/Lpump - 2*np.pi*ne_GS(Lidler)/Lidler - 2*np.pi*net_GS(LTHz,CA)/LTHz
          PM_BBO[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
          for j, f in enumerate (THz_freq):
              if math.isnan(PM_BBO[i,j]):
                  PM_BBO[i,j]=0
         
    #### COLOR PLOT (optional)
    if it==0:
        fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
        im = ax.pcolor(THz_freq, Angle, PM_BBO, cmap='jet')
        ax.set_ylabel('angle (deg)', fontsize=s) 
        ax.set_xlabel('signal frequency (THz)', fontsize=s)   
        ax.tick_params(labelsize=s)
        secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
        secax.set_xlabel('signal wavelength (µm)', fontsize=s)
        secax.set_xticks((0.5, 1 ,1.5, 2, 3, 4))
        secax.tick_params(labelsize=s)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        fig.suptitle('BBO, SHG, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
        fig.savefig(r'D:\phase matching OPA\python\results\BBO_SHG_1030nm_'+ str(C_thick)+'_plot1.png', dpi=300, bbox_inches='tight')
        


        fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
        im = ax.pcolor(idler_freq, Angle, PM_BBO, cmap='jet')
        ax.set_ylabel('angle (deg)', fontsize=s) 
        ax.set_xlabel('idler frequency (THz)', fontsize=s)   
        ax.tick_params(labelsize=s)
        secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
        secax.set_xlabel('idler wvl (µm)', fontsize=s)
        secax.set_xticks((0.2, 0.5, 1 ,1.5, 2, 3, 4))
        secax.tick_params(labelsize=s)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        fig.suptitle('BBO, SHG, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
        fig.savefig(r'D:\phase matching OPA\python\results\BBO_SHG_1030nm_'+ str(C_thick)+'_plot1_bis.png', dpi=300, bbox_inches='tight')
        
    

    
    # spectrum at the angle which maximizes the total signal (over the full frequency range)
    Profile_THz=np.zeros(np.size(Angle))
    for i, a in enumerate (Angle):
        Profile_THz[i] = np.sum(PM_BBO[i,:])
    Max_angle=np.where(Profile_THz==np.max(Profile_THz))
    print('index ' + str(Max_angle[0][0]))
    print('angle ' + str(Angle[Max_angle[0][0]]) + '°')
    
    
    fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
    ax.plot(THz_freq,PM_BBO[Max_angle[0][0],:], color='navy', label= str(round(Angle[Max_angle[0][0]],2)) +'°')
    ax.legend()
    ax.set_ylabel('Relative intensity', fontsize=s) 
    ax.set_xlabel('signal frequency (THz)', fontsize=s) 
    fig.text(0,1.1,'BBO SHG: pump 1030nm (o) + signal (o) --> idler (e) \n Spectrum of the amplified signal for the optimum angle \n which gives the highest integrated signal')
    secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((0.5, 1, 1.5,  2))
    fig.savefig(r'D:\phase matching OPA\python\results\BBO_SHG_1030nm_'+ str(C_thick)+'_plot2.png', dpi=300, bbox_inches='tight')
    
    
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
    #ax.set_xticks((1,2,3,4,5,6))
    #ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
    fig.text(0,1.1,'BBO SHG: pump 3.4µm (o) + signal (o) --> idler (e) ')
    fig.savefig(r'D:\phase matching OPA\python\results\BBO_SHG_1030nm_'+ str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    ## Ici, on s'intéresse à 1.03µm
    ## Quel est l'angle optimum pour 1.03µm ?
    
    Target = round(lam2freq(1.03),1)  #target frequency in THz
    idx=[]
    prec=1e-4
    idx = np.where(np.isclose(THz_freq, Target,prec))
    if np.size(idx) == 0:
        while np.size (idx) ==0 :
            prec = prec*10
            idx = np.where(np.isclose(THz_freq, Target,prec))
            print(idx)
            if prec >= 1:
                print('error, cannot find the index')
                break
    
    idx=int(np.mean(idx))
    # print(THz_freq[idx])
    # print(Target)
    idxmax = np.where(PM_BBO[:, idx] == np.max(PM_BBO[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM_BBO[:, idx], color='red')
    #ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
    fig.text(0,1.2,'BBO SHG: pump 1.03µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 1.03µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM_BBO[:, idx], ra1,ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Output signal')
    fig.savefig(r'D:\phase matching OPA\python\results\AGS_SHG_3400nm_'+ str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
    
    
    
    
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM_BBO[idxmax, :], color='red')
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    # ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
    fig.text(0,1.2,'BBO SHG: pump 1.03µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
    ax.set_xlabel('Signal frequency (THz)')
    ax.set_ylabel('Output signal')
    secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((0.5, 1 ,1.5, 2, 3, 4 , 5))
    
    
    
    FWHM_freq = round(FWHM(THz_freq, PM_BBO[idxmax, :],r1,r2),4)
    print(str(FWHM_freq) + ' THz FWHM' )
    fig.text(0,1.4,str(FWHM_freq) + ' THz FWHM', fontsize=s)
    
    #assuming gaussian shape
    print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')
    fig.text(0,1.3,'Assuming Gaussian shape: ' + str("{:10.4}".format(2*np.log(2)/(np.pi*FWHM_freq*1e12))) + ' s FWHM', fontsize=s)
    
    fig.text(0,1.5,'Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    
    fig.savefig(r'D:\phase matching OPA\python\results\BBO_SHG_1030nm_'+ str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
    
    
    
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM(Angle, PM_BBO[:, idx], ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    
        
        
    
    
    
fig_F=plt.figure()
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
fig_F.savefig(r'D:\phase matching OPA\python\results\BBO_SHG_1030nm_plotFINAL.png', dpi=300, bbox_inches='tight')
  





