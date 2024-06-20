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



#FWHM function
def FWHM (x, y):
    ymax = np.max(y)
    print('\n', ymax*0.5 , 'max*0.5')
    for i, value in enumerate (y):
        if value==ymax:
            index_max=i
            print(index_max)
    for i, value in enumerate ( y[0:index_max] ):
        if value >= 0.5*y[index_max]:
            low = i
            break
    print(low)
    print(y[low] , 'y[low]')
    for i, value in enumerate (y):
        if i > index_max+1:
            if value <= 0.5*y[index_max]:
                high = i
                break
    print(high)
    print(y[high] , 'y[high]', '\n')
    return (x[high] - x[low])



freq_min=20
freq_max=70





def FWHM_res(x,y,a=freq_min, b=freq_max):
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




#%% CASE 1 : LGS crystal LiGaS2, 10µm generation using 1.03µm (and 1.1483µm)
#biaxial crystal
#collinear geometry
folder=r'D:\phase matching OPA\python\results\LGS_10um_option1'

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print(' LGS crystal 10µm (o) + 1.1483µm (o) = 1.03µm (e), collinear geometry, Plane XZ')

c=0.299792458   #µm/fs
C_thick=500 #µm
THz_freq=np.arange(freq_min,freq_max,0.02)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


# Optics Letters 42, 21 (2017), Kato et al:
#nx
def nxLGS(l):
    return np.sqrt(5.40285 + 0.11050/(l**2 - 0.03238) + 460.11/(l**2 - 422.74))
#ny = no
def noLGS(l): 
    return np.sqrt(5.62849 + 0.12052/(l**2-0.03429) + 485.31/(l**2 - 424.93) )
#nz
def nzLGS(l):
    return np.sqrt(5.77851 + 0.12325/(l**2 - 0.03465) + 580.13/(l**2 - 451.62)  )

def neLGS(l,a):
    return nxLGS(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (nxLGS(l)*np.tan(a)/nzLGS(l))**2 )  )


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
Angle = np.arange(0,90,0.02)

PM1 = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
#    CA=b*np.pi/180
  #   print(a)
     CA = b*np.pi/180
     Dk = 2*np.pi*neLGS(Lpump,CA)/Lpump - 2*np.pi*noLGS(Lidler)/Lidler - 2*np.pi*noLGS(LTHz)/LTHz
     PM1[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
    
    
  

fig, ax = plt.subplots(figsize=(4,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM1, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(10,65)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.tick_params(labelsize=s)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.suptitle('LGS plane XZ, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm \n signal 10µm (o) + idler 1.1483µm (o) = pump 1.03µm (e)', fontsize=s)
fig=plt.savefig( folder + '\\LGS_DFG_10000nm_1148nm_PMcolor.png', dpi=300, bbox_inches='tight')
  


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM1[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')

plt.figure()
plt.plot(THz_freq,PM1[Max_angle[0][0],:], label= str(Angle[Max_angle[0][0]]) +'°')
plt.legend()
plt.title('LGS plane XZ, Spectrum of the amplified signal for a given angle \n (angle for which the max integrated signal is obtained) \n Effective thickness =' + str(C_thick) + 'µm')
plt.savefig( folder + '\\LGS_DFG_10000nm_1148nm_broadest_spec.png', dpi=300, bbox_inches='tight')
  

################################################################################################################
#%##optimum angle for 10µm?

# # # angle (deg) range for the FWHM calculation:
# ra1= 53.5
# ra2= 55.5

Target = round(lam2freq(10),1)  #target frequency in THz
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
idxmax = np.where(PM1[:, idx] == np.max(PM1[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM1[:, idx], color='red')
#fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
fig.text(0,1.1,'Angular profile at 10µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°) \n Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM1[:, idx]))*1000,3)) + 'mrad FWHM', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
plt.savefig(folder + '\\LGS_DFG_10000nm_1148nm_angular_profile_for_10um.png', dpi=300, bbox_inches='tight')
  



#%
fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
ax.plot(THz_freq, PM1[idxmax, :], color='red')
#ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
fig.text(0,1.1,'LGS 10µm signal with 1.03µm pump \n Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'° \n Effective thickness =' + str(C_thick) + 'µm ', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((2,4,6,8,10,12,14))
plt.savefig( folder + '\\LGS_DFG_10000nm_1148nm_spectrum_for_10um.png', dpi=300, bbox_inches='tight')
  


###############################################################################% thickness scan
## frequency (THz) range for the FWHM calculation:
r1 =22
r2=65

# angle (deg) range for the FWHM calculation:
ra1= 40
ra2= 60



Angle_FWHM=[]
Duration_FWHM=[]

C_thicks = [500, 1000,2000]
for it, C_thick in enumerate (C_thicks):

    if C_thick == 500:
        r2=65
    else:
        r2=40
             
    for i, b in enumerate (Angle):
         CA = b*np.pi/180
         Dk = 2*np.pi*neLGS(Lpump,CA)/Lpump - 2*np.pi*noLGS(Lidler)/Lidler - 2*np.pi*noLGS(LTHz)/LTHz
         PM1[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
        
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
#    fig.text(0,1.1,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e) ')
    fig.savefig(folder + '\\LGS_10mic_'+str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    Target = round(lam2freq(10),1)  #target frequency in THz
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
    idxmax = np.where(PM1[:, idx] == np.max(PM1[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM1[:, idx], color='red')
 #   fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 10µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM_res(Angle, PM1[:, idx], ra1, ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_xlim(ra1,ra2)
    ax.set_ylabel('Output signal')
    fig.savefig(folder +'\\LGS_10mic_'+str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
  
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM1[idxmax, :], color='red')
 #   if C_thick != 500:
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    # ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
 #   fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
    ax.set_xlabel('Signal frequency (THz)')
    ax.set_ylabel('Output signal')
    secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    
    
    # if C_thick == 500:
    #     FWHM_freq = round(FWHM(THz_freq, PM1[idxmax, :]),4)
    # else:
    FWHM_freq = round(FWHM_res(THz_freq, PM1[idxmax, :],r1,r2),4)
    print(str(FWHM_freq) + ' THz FWHM' )
    fig.text(0,1.4,str(FWHM_freq) + ' THz FWHM', fontsize=s)
    
    #assuming gaussian shape
    print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')
    fig.text(0,1.3,'Assuming Gaussian shape: ' + str("{:10.4}".format(2*np.log(2)/(np.pi*FWHM_freq*1e12))) + ' s FWHM', fontsize=s)
    
    fig.text(0,1.5,'Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    
    fig.savefig(folder + '\\LGS_10mic_'+str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
    
    
    
    
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM_res(Angle, PM1[:, idx],ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    

    
fig_F=plt.figure(figsize=(4,3))
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
fig_F.savefig(folder + '\\LGS_10mic_plotFINAL.png', dpi=300, bbox_inches='tight')
  

















































#%% CASE 2 : LGS crystal LiGaS2, 10µm generation using 1.03µm (and 1.1483µm) NOT FINISHED!!!!!!!!!!!!!!
### theta = 90° (meaning the propagation direction is in the plane x,y)
### phi = 37.3° 
### the idler (1.1483µm) is polarized along z --> sees the index nz only  --> nz=no here
### the pump (1.03µm) and the signal (10µm) are polarized in the (x,y) plane --> see a mix of nx and ny with angle phi  (here CA)  --> nx,ny,phi = ne(phi)  here


folder=r'D:\phase matching OPA\python\results\LGS_10um_thick500um_option2'

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print(' LGS crystal 10µm (e) + 1.1483µm (o) = 1.03µm (e), collinear geometry, Plane XY')

c=0.299792458   #µm/fs
C_thick=500 #µm
THz_freq=np.arange(freq_min,freq_max,0.02)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


# Optics Letters 42, 21 (2017), Kato et al:
#nx
def nxLGS(l):
    return np.sqrt(5.40285 + 0.11050/(l**2 - 0.03238) + 460.11/(l**2 - 422.74))
#ny = no
def nyLGS(l): 
    return np.sqrt(5.62849 + 0.12052/(l**2-0.03429) + 485.31/(l**2 - 424.93) )
#nz
def noLGS(l):
    return np.sqrt(5.77851 + 0.12325/(l**2 - 0.03465) + 580.13/(l**2 - 451.62)  )

def neLGS(l,a):
    return nyLGS(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (nyLGS(l)*np.tan(a)/nxLGS(l))**2 )  )


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
Angle = np.arange(0,90,0.02)

PM1 = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
#    CA=b*np.pi/180
  #   print(a)
     CA = b*np.pi/180
     Dk = 2*np.pi*neLGS(Lpump,CA)/Lpump - 2*np.pi*noLGS(Lidler)/Lidler - 2*np.pi*neLGS(LTHz,CA)/LTHz
     PM1[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
    
    
    


fig, ax = plt.subplots(figsize=(4,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM1, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(10,65)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.tick_params(labelsize=s)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.suptitle('LGS plane XY, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm \n signal 10µm (e) + idler 1.1483µm (o) = pump 1.03µm (e)', fontsize=s)
fig=plt.savefig( folder + '\\LGS_DFG_10000nm_1148nm_PMcolor.png', dpi=300, bbox_inches='tight')
  





Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM1[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')

plt.figure()
plt.plot(THz_freq,PM1[Max_angle[0][0],:], label= str(Angle[Max_angle[0][0]]) +'°')
plt.legend()
plt.title('LGS plane XY, Spectrum of the amplified signal for a given angle \n (angle for which the max integrated signal is obtained) \n Effective thickness =' + str(C_thick) + 'µm')
plt.savefig( folder + '\\LGS_DFG_10000nm_1148nm_broadest_spec.png', dpi=300, bbox_inches='tight')
  
####################################################% optimum angle for 10µm?

# # # angle (deg) range for the FWHM calculation:
# ra1= 53.5
# ra2= 55.5

Target = round(lam2freq(10),1)  #target frequency in THz
idx=[]; idxmax=[];
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
idxmax = np.where(PM1[:, idx] == np.max(PM1[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM1[:, idx], color='red')
#fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
fig.text(0,1.1,'Angular profile at 10µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°) \n Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM1[:, idx]))*1000,3)) + 'mrad FWHM', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
plt.savefig(folder + '\\LGS_DFG_10000nm_1148nm_angular_profile_for_10um.png', dpi=300, bbox_inches='tight')
  


#############################################################################################

fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
ax.plot(THz_freq, PM1[idxmax, :], color='red')
#ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
fig.text(0,1.1,'LGS 10µm signal with 1.03µm pump \n Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'° \n Effective thickness =' + str(C_thick) + 'µm ', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((2,4,6,8,10,12,14))
plt.savefig( folder + '\\LGS_DFG_10000nm_1148nm_spectrum_for_10um.png', dpi=300, bbox_inches='tight')


###########################################################################################%% thickness scan
## frequency (THz) range for the FWHM calculation:
r1 =22
r2=60

# angle (deg) range for the FWHM calculation:
ra1= 25
ra2= 45



Angle_FWHM=[]
Duration_FWHM=[]

C_thicks = [500,1000,2000]
for it, C_thick in enumerate (C_thicks):
    if C_thick != 2000:
        r2=60
    else:
        r2=35

             
    for i, b in enumerate (Angle):
         CA = b*np.pi/180
         Dk = 2*np.pi*neLGS(Lpump,CA)/Lpump - 2*np.pi*noLGS(Lidler)/Lidler - 2*np.pi*neLGS(LTHz,CA)/LTHz
         PM1[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
        
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
#    fig.text(0,1.1,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e) ')
    fig.savefig(folder + '\\LGS_10mic_'+str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    Target = round(lam2freq(10),1)  #target frequency in THz
    idx=[]; idxmax=[];
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
    idxmax = np.where(PM1[:, idx] == np.max(PM1[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM1[:, idx], color='red')
 #   fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 10µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM_res(Angle, PM1[:, idx], ra1, ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_xlim(ra1,ra2)
    ax.set_ylabel('Output signal')
    fig.savefig(folder +'\\LGS_10mic_'+str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
  
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM1[idxmax, :], color='red')
   # if C_thick != 500:
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    # ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
 #   fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
    ax.set_xlabel('Signal frequency (THz)')
    ax.set_ylabel('Output signal')
    secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    
    
    # if C_thick == 500:
    #     FWHM_freq = round(FWHM(THz_freq, PM1[idxmax, :]),4)
    # else:
    #     
    FWHM_freq = round(FWHM_res(THz_freq, PM1[idxmax, :],r1,r2),4)
    print(str(FWHM_freq) + ' THz FWHM' )
    fig.text(0,1.4,str(FWHM_freq) + ' THz FWHM', fontsize=s)
    
    #assuming gaussian shape
    print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')
    fig.text(0,1.3,'Assuming Gaussian shape: ' + str("{:10.4}".format(2*np.log(2)/(np.pi*FWHM_freq*1e12))) + ' s FWHM', fontsize=s)
    
    fig.text(0,1.5,'Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    
    fig.savefig(folder + '\\LGS_10mic_'+str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
    
    
    
    
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM_res(Angle, PM1[:, idx],ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    
    
      
    
fig_F=plt.figure(figsize=(4,3))
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
fig_F.savefig(folder + '\\LGS_10mic_plotFINAL.png', dpi=300, bbox_inches='tight')
  
































#%% SHG 10µm + 10µm = 5µm using GaSe (GS)
# SHG for 10000 nm --> 5000nm
## GaSe (called GS in SNLO) just for ref, bc absorption is too high <-- YUYA
# refractiveindex.info:
print ('GaSe/GS, SHG @10000nm')
    

folder = r'D:\phase matching OPA\python\results\GaSe_SHG_10000nm\\'


Angle = np.arange(1,25,0.01)
c=0.299792458   #µm/fs
C_thick=4000 #µm
freq_min=10#freq2lam(3.4+1)
freq_max=180#freq2lam(3.4-1)
THz_freq=np.arange(freq_min,freq_max,0.01)      #
Lpump=10 #µm   
LTHz = 1/(THz_freq/c)*1000

## frequency (THz) range for the FWHM calculation:
r1 =20
r2=60

# angle (deg) range for the FWHM calculation:
ra1= 8
ra2= 20


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




def no_GS(l):
    return np.sqrt(10.6409 + 0.3788/((l**2) - 0.1232) + 7090.7/(l**2 -2216.3))
def ne_GS(l):
    return np.sqrt(8.2477 + 0.2881/((l**2) - 0.1669) +4927.5 /((l)**2-1990.1))
def net_GS(l,tet): 
    return np.sqrt(1/( np.cos(tet)**2/no_GS(l)**2 + np.sin(tet)**2/ne_GS(l)**2 ))



#%Type 1 phase matching for LN
print('o + o --> e')

Pump_freq=c/Lpump*1000
idler_freq=Pump_freq+THz_freq
Lidler=1/(idler_freq/c)*1000
PM_GS = np.zeros((np.size(Angle), np.size(LTHz)))




Angle_FWHM=[]
Duration_FWHM=[]

C_thicks = [500,1000,2000]
for it, C_thick in enumerate (C_thicks):

    
    for i, b in enumerate (Angle):
          CA = b*np.pi/180
          Dk = 2*np.pi*no_GS(Lpump)/Lpump - 2*np.pi*net_GS(Lidler,CA)/Lidler + 2*np.pi*no_GS(LTHz)/LTHz
          #Dk = 2*np.pi*net_GS(Lpump,CA)/Lpump - 2*np.pi*ne_GS(Lidler)/Lidler - 2*np.pi*net_GS(LTHz,CA)/LTHz
          PM_GS[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
          for j, f in enumerate (THz_freq):
              if math.isnan(PM_GS[i,j]):
                  PM_GS[i,j]=0
      
    if it==0:
        
        ##" OPTIONAL COLOR PLOT
        fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
        im = ax.pcolor(THz_freq, Angle, PM_GS, cmap='jet')
        ax.set_ylabel('angle (deg)', fontsize=s) 
        ax.set_xlabel('signal frequency (THz)', fontsize=s)   
        ax.tick_params(labelsize=s)
        secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
        secax.set_xlabel('signal wavelength (µm)', fontsize=s)
        secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
        secax.tick_params(labelsize=s)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        fig.suptitle('GaSe (GS), SHG, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
        fig.savefig(folder +str(C_thick)+'_plot1.png', dpi=300, bbox_inches='tight')
        
        
    # spectrum at the angle which maximizes the total signal (over the full frequency range)
    Profile_THz=np.zeros(np.size(Angle))
    for i, a in enumerate (Angle):
        Profile_THz[i] = np.sum(PM_GS[i,:])
    Max_angle=np.where(Profile_THz==np.max(Profile_THz))
    print('index ' + str(Max_angle[0][0]))
    print('angle ' + str(Angle[Max_angle[0][0]]) + '°')
    
    
    fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
    ax.plot(THz_freq,PM_GS[Max_angle[0][0],:], color='navy', label= str(round(Angle[Max_angle[0][0]],2)) +'°')
    ax.legend()
    ax.set_ylabel('Relative intensity', fontsize=s) 
    ax.set_xlabel('signal frequency (THz)', fontsize=s) 
    fig.text(0,1.1,'GaSe (GS) SHG: pump' + str(Lpump) + 'µm (o) + signal (o) --> idler (e) \n Spectrum of the amplified signal for the optimum angle \n which gives the highest integrated signal')
    secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    fig.savefig(folder +str(C_thick)+'_plot2.png', dpi=300, bbox_inches='tight')
    
    
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
    #ax.set_xticks((1,2,3,4,5,6))
    #ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
    fig.text(0,1.1,'GaSe (GS) SHG: pump 10µm (o) + 10µm signal (o) --> 5µm idler (e) ')
    fig.savefig(folder +str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    ## Ici, on s'intéresse à 10µm
    ## Quel est l'angle optimum pour 10µm ?
    
    Target = round(lam2freq(10),1)  #target frequency in THz
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
    idxmax = np.where(PM_GS[:, idx] == np.max(PM_GS[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM_GS[:, idx], color='red')
    fig.text(0,1.2,'GaSe (GS) SHG: pump '+ str(Lpump) +'µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 3.4µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM_GS[:, idx], ra1,ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Output signal')
    fig.savefig(folder +str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
    
    
    
    
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM_GS[idxmax, :], color='red')
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    fig.text(0,1.2,'GaSe (GS) SHG: pump '+ str(Lpump) + 'µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
    ax.set_xlabel('Signal frequency (THz)')
    ax.set_ylabel('Output signal')
    secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    
    FWHM_freq = round(FWHM(THz_freq, PM_GS[idxmax, :],r1,r2),4)
    print(str(FWHM_freq) + ' THz FWHM' )
    fig.text(0,1.4,str(FWHM_freq) + ' THz FWHM', fontsize=s)
    
    #assuming gaussian shape
    print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')
    fig.text(0,1.3,'Assuming Gaussian shape: ' + str("{:10.4}".format(2*np.log(2)/(np.pi*FWHM_freq*1e12))) + ' s FWHM', fontsize=s)
    fig.text(0,1.5,'Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    fig.savefig(folder +str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM(Angle, PM_GS[:, idx], ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    
        
        
    
    
    
fig_F=plt.figure()
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
fig_F.savefig(folder + 'plotFINAL.png', dpi=300, bbox_inches='tight')
  





#%% SHG 10µm + 10µm = 5µm using GaSe (GS) TYPE 2
# SHG for 10000 nm --> 5000nm
## GaSe (called GS in SNLO) just for ref, bc absorption is too high <-- YUYA
# refractiveindex.info:
print ('GaSe/GS, SHG @10000nm')
    

folder = r'D:\phase matching OPA\python\results\GaSe_SHG_10000nm_type2\\'


Angle = np.arange(1,40,0.04)
c=0.299792458   #µm/fs
C_thick=4000 #µm
freq_min=10#freq2lam(3.4+1)
freq_max=180#freq2lam(3.4-1)
THz_freq=np.arange(freq_min,freq_max,0.04)      #
Lpump=10 #µm   
LTHz = 1/(THz_freq/c)*1000

## frequency (THz) range for the FWHM calculation:
r1 =11
r2=120

# angle (deg) range for the FWHM calculation:
ra1= 10
ra2= 30


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




def no_GS(l):
    return np.sqrt(10.6409 + 0.3788/((l**2) - 0.1232) + 7090.7/(l**2 -2216.3))
def ne_GS(l):
    return np.sqrt(8.2477 + 0.2881/((l**2) - 0.1669) +4927.5 /((l)**2-1990.1))
def net_GS(l,tet): 
    return np.sqrt(1/( np.cos(tet)**2/no_GS(l)**2 + np.sin(tet)**2/ne_GS(l)**2 ))



#%Type 1 phase matching for LN
print('o + e --> e')

Pump_freq=c/Lpump*1000
idler_freq=Pump_freq+THz_freq
Lidler=1/(idler_freq/c)*1000
PM_GS = np.zeros((np.size(Angle), np.size(LTHz)))




Angle_FWHM=[]
Duration_FWHM=[]

C_thicks = [500,1000,2000]
for it, C_thick in enumerate (C_thicks):

    if C_thick == 2000:
        r2=46
    else:
        r2=120

    for i, b in enumerate (Angle):
          CA = b*np.pi/180
          Dk = 2*np.pi*no_GS(Lpump)/Lpump - 2*np.pi*net_GS(Lidler,CA)/Lidler + 2*np.pi*net_GS(LTHz,CA)/LTHz
          #Dk = 2*np.pi*net_GS(Lpump,CA)/Lpump - 2*np.pi*ne_GS(Lidler)/Lidler - 2*np.pi*net_GS(LTHz,CA)/LTHz
          PM_GS[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
          for j, f in enumerate (THz_freq):
              if math.isnan(PM_GS[i,j]):
                  PM_GS[i,j]=0
      
    if it==0:
        
        ##" OPTIONAL COLOR PLOT
        fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
        im = ax.pcolor(THz_freq, Angle, PM_GS, cmap='jet')
        ax.set_ylabel('angle (deg)', fontsize=s) 
        ax.set_xlabel('signal frequency (THz)', fontsize=s)   
        ax.tick_params(labelsize=s)
        secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
        secax.set_xlabel('signal wavelength (µm)', fontsize=s)
        secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
        secax.tick_params(labelsize=s)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        fig.suptitle('GaSe (GS), SHG, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
        fig.savefig(folder +str(C_thick)+'_plot1.png', dpi=300, bbox_inches='tight')
        
        
    # spectrum at the angle which maximizes the total signal (over the full frequency range)
    Profile_THz=np.zeros(np.size(Angle))
    for i, a in enumerate (Angle):
        Profile_THz[i] = np.sum(PM_GS[i,:])
    Max_angle=np.where(Profile_THz==np.max(Profile_THz))
    print('index ' + str(Max_angle[0][0]))
    print('angle ' + str(Angle[Max_angle[0][0]]) + '°')
    
    
    fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
    ax.plot(THz_freq,PM_GS[Max_angle[0][0],:], color='navy', label= str(round(Angle[Max_angle[0][0]],2)) +'°')
    ax.legend()
    ax.set_ylabel('Relative intensity', fontsize=s) 
    ax.set_xlabel('signal frequency (THz)', fontsize=s) 
    fig.text(0,1.1,'GaSe (GS) SHG: pump' + str(Lpump) + 'µm (o) + signal (e) --> idler (e) \n Spectrum of the amplified signal for the optimum angle \n which gives the highest integrated signal')
    secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    fig.savefig(folder +str(C_thick)+'_plot2.png', dpi=300, bbox_inches='tight')
    
    
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
    #ax.set_xticks((1,2,3,4,5,6))
    #ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
    fig.text(0,1.1,'GaSe (GS) SHG: pump 10µm (o) + 10µm signal (e) --> 5µm idler (e) ')
    fig.savefig(folder +str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    ## Ici, on s'intéresse à 10µm
    ## Quel est l'angle optimum pour 10µm ?
    
    Target = round(lam2freq(10),1)  #target frequency in THz
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
    idxmax = np.where(PM_GS[:, idx] == np.max(PM_GS[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM_GS[:, idx], color='red')
    fig.text(0,1.2,'GaSe (GS) SHG: pump '+ str(Lpump) +'µm (o) + signal (e) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 3.4µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM_GS[:, idx], ra1,ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Output signal')
    fig.savefig(folder +str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
    
    
    
    
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM_GS[idxmax, :], color='red')
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    fig.text(0,1.2,'GaSe (GS) SHG: pump '+ str(Lpump) + 'µm (o) + signal (e) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
    ax.set_xlabel('Signal frequency (THz)')
    ax.set_ylabel('Output signal')
    secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    
    FWHM_freq = round(FWHM(THz_freq, PM_GS[idxmax, :],r1,r2),4)
    print(str(FWHM_freq) + ' THz FWHM' )
    fig.text(0,1.4,str(FWHM_freq) + ' THz FWHM', fontsize=s)
    
    #assuming gaussian shape
    print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')
    fig.text(0,1.3,'Assuming Gaussian shape: ' + str("{:10.4}".format(2*np.log(2)/(np.pi*FWHM_freq*1e12))) + ' s FWHM', fontsize=s)
    fig.text(0,1.5,'Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    fig.savefig(folder +str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM(Angle, PM_GS[:, idx], ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    
        
        
    
    
    
fig_F=plt.figure()
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
fig_F.savefig(folder + 'plotFINAL.png', dpi=300, bbox_inches='tight')



#%%  spatial walk off
###########    tan rho = - 1 /ne dne/dtheta



  