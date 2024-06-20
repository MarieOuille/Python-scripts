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


 
def lam2freq(x):
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



#%% KTA literature overview
# KTA_thick=(2,2,10,1.4,8,5,2,2) #mm
# pump_duration=(70,53,85,6000,1000,1400,1100,100)  #fs
# tum_duration=(115,180,350,52,125,79,73,70) #fs

# plt.figure(figsize=(4,4))
# plt.scatter(KTA_thick, tum_duration)
# plt.xlabel('KTA thickness (mm)')
# plt.ylabel('3µm pulse duration (fs)')
# plt.savefig(r'D:\phase matching OPA\python\results\KTA_literature\\graph1_bis.png', bbox_inches='tight')



# plt.figure(figsize=(4,4))
# plt.scatter(KTA_thick, pump_duration, color='red')
# plt.xlabel('KTA thickness (mm)')
# plt.ylabel('pump pulse duration (fs)')
# plt.savefig(r'D:\phase matching OPA\python\results\KTA_literature\\graph2_bis.png', bbox_inches='tight')



# # #%% case 1 : AGS crystal, 10µm generation using a 0.8µm laser
# # collinear geometry
# # uniaxial crystal

# print('\n AGS crystal (collinear geometry, uniaxial crystal)')
# Lp = 0.8 #pump wavelength (µm)        extraordinary
# Li = 10 #idler wavelength (µm)        ordinary
# Ls = (1/Lp - 1/Li)**(-1) #signal wavelength (µm)          ordinary
# print('Lp=' + str(Lp) +'µm, Ls=' + str(round(Ls,4)) + 'µm, Li=' + str(Li) + 'µm')
# #print(r'1/Ls + 1/Li - 1/Lp =' + str(1/Ls + 1/Li - 1/Lp)) #check

# #refractive index of AGS depending on the angle (refractiveindex.info, AgGaS2, 293.15 Kelvin)
# def no(l):
#     return np.sqrt(5.7975 + 0.2311/((l**2) - 0.0688) - 0.00257*(l)**2)
# def ne(l):
#     return np.sqrt(5.5436 + 0.2230/((l**2) - 0.0946) - 0.00261*(l)**2)
# lams = np.arange(0.6,12,0.1)
# plt.figure()
# plt.title(r'AGS crystal, $\theta$ = 0°')
# plt.plot(lams, no(lams), label='no', color='black',linestyle='--')
# plt.plot(lams, ne(lams), label='ne', color='black')
# plt.legend()
# plt.xlabel(r'$\lambda$ (µm)')
# plt.ylabel('refractive index')
# plt.axvline(Lp, color='green')
# plt.axvline(Ls, color='red')
# plt.axvline(Li, color='purple')

# def net(l,tet):
#     return np.sqrt(1/( np.cos(tet)**2/no(l)**2 + np.sin(tet)**2/ne(l)**2 ))



# # teta=np.deg2rad(42.6)   #angle

# # plt.figure()
# # plt.plot(lams, no(lams), label='o')
# # plt.plot(lams, net(lams,teta), label='e')
# # plt.title(r'$\theta = $' + str(np.rad2deg(teta)) +'°')

# #Type 1 phase matching : e (pump, 0.8µm) --> o + o 
# def y1(x) :
#     return net(Lp,x)/Lp - no(Ls)/Ls - no(Li)/Li 
# # plt.figure()
# # plt.plot(np.rad2deg(np.arange(-3.14/2,3.14/2,0.001)), y1(np.arange(-3.14/2,3.14/2,0.001)))
# result1 = scipy.optimize.fsolve(y1,0.6)
# #print(result)
# print('Type 1 phase matching angle: ' + str(np.round(np.rad2deg(result1[0]),1)) + '°')



# print('net(Lp,theta)=' + str(np.round(net(Lp,result1[0]),3)))
# print('no(Li)=' + str(np.round(no(Li),3)))
# print('no(Ls)=' + str(np.round(no(Ls),3)))


# #Type 2 phase matching : e --> o + e (10µm)
# def y2(x) :
#     return net(Lp,x)/Lp - no(Ls)/Ls - net(Li,x)/Li 
# # plt.figure()
# # plt.plot(np.rad2deg(np.arange(-3.14/2,3.14/2,0.001)), y2(np.arange(-3.14/2,3.14/2,0.001)))
# result2 = scipy.optimize.fsolve(y2,0.69)
# #print(result)
# print('Type 2 phase matching angle: ' + str(round(np.rad2deg(result2[0]),1)) +'°' )

# print('net(Lp,theta)=' + str(np.round(net(Lp,result2[0]),3)))
# print('net(Li, theta)=' + str(np.round(net(Li, result2[0]),3)))
# print('no(Ls)=' + str(np.round(no(Ls),3))+'\n')


# # plt.figure()
# # plt.title(r'AGS crystal, $\theta$ =' + str(np.rad2deg(result[0]))+'°')
# # plt.plot(lams, no(lams), label='no', color='black',linestyle='--')
# # plt.plot(lams, net(lams,np.rad2deg(result[0])), label='ne', color='black')
# # plt.legend()
# # plt.xlabel(r'$\lambda$ (µm)')
# # plt.ylabel('refractive index')
# # plt.axvline(Lp, color='green')
# # plt.axvline(Ls, color='red')
# # plt.axvline(Li, color='purple')




#%% CASE 2 : LGS crystal LiGaS2
#biaxial crystal
#collinear geometry

print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print(' LGS crystal type 1 (collinear geometry)')

c=0.299792458   #µm/fs
C_thick=500 #µm
freq_min=15
freq_max=50
THz_freq=np.arange(freq_min,freq_max,0.02)
Lpump=1.55 #µm
LTHz = 1/(THz_freq/c)*1000

LGS_o_A1 = 4.478907;
LGS_o_B1 = 0.120426;
LGS_o_C1 = 0.034616;
LGS_o_D1 = 0.0035119;

LGS_X_A1 = 4.326834;
LGS_X_B1 = 0.1030907;
LGS_X_C1 = 0.0309876;
LGS_X_D1 = 0.0037015;

LGS_Z_A1 = 4.493881;
LGS_Z_B1 = 0.1177452;
LGS_Z_C1 = 0.0337004;
LGS_Z_D1 = 0.0037767;


def noLGS(l):
    return np.sqrt(LGS_o_A1 + LGS_o_B1/(l**2 -LGS_o_C1) - LGS_o_D1*l**2)
# np.sqrt(5.40285 + 0.11050/(l**2 - 0.03238) + 460.11/(l**2 - 422.74))
def nzLGS(l): 
    return np.sqrt(LGS_Z_A1 + LGS_Z_B1/(l**2 -LGS_Z_C1) - LGS_Z_D1*l**2)
#np.sqrt(5.62849 + 0.12052/(l**2-0.03429) + 485.31/(l**2 - 424.93) )
def nxLGS(l):
    return np.sqrt(LGS_X_A1 + LGS_X_B1/(l**2 -LGS_X_C1) - LGS_X_D1*l**2)
#np.sqrt(5.77851 + 0.12325/(l**2 - 0.03465) + 580.13/(l**2 - 451.62)  )

def neLGS(l,a):
    return nxLGS(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (nxLGS(l)*np.tan(a)/nzLGS(l))**2 )  )

Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
Angle = np.arange(0,90,0.05)

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
fig.suptitle('LGS type 1, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM1[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')

plt.figure()
plt.plot(THz_freq,PM1[Max_angle[0][0],:], label= str(Angle[Max_angle[0][0]]) +'°')
plt.legend()
plt.title('LGS type 1, Spectrum of the amplified signal for a given angle')





















#%%Lithium Niobate LN  undoped   negative unixial
print ('\n %%%%%%%%%%%%%%%%%%%%%%%%% \n Lithium Niobate ')

def ne_LN(l):
    return np.sqrt(1+ 2.9804*l**2/(l**2 -0.02047) + 0.5981*l**2/(l**2 -0.0666)+8.9543*l**2/(l**2 -416.08))

def no_LN(l):
    return np.sqrt(1 + 2.6734*l**2/(l**2 - 0.01764) + 1.2290*l**2/(l**2 - 0.05914) + 12.614*l**2/(l**2 - 474.60))

def ne_LN_tet (a,l):
    return np.sqrt( 1/(np.cos(a)**2/no_LN(l)**2  + np.sin(a)**2/ne_LN(l)**2))


c=0.299792458   #µm/fs
C_thick=500 #µm
freq_min=50
freq_max=400
THz_freq=np.arange(freq_min,freq_max,0.04)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


#%Type 1 phase matching for LN
print('Type 1 phase matching')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM1_LN = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
#    CA=b*np.pi/180
  #   print(a)
     CA = b*np.pi/180
     Dk = 2*np.pi*ne_LN_tet(Lpump,CA)/Lpump - 2*np.pi*no_LN(Lidler)/Lidler - 2*np.pi*no_LN(LTHz)/LTHz
     PM1_LN[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
     for j, f in enumerate (THz_freq):
         if math.isnan(PM1_LN[i,j]):
             PM1_LN[i,j]=0
    

fig, ax = plt.subplots(figsize=(8,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM1_LN, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(0,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.suptitle('Undoped LN type 1, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM1_LN[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')

plt.figure()
plt.plot(THz_freq,PM1_LN[Max_angle[0][0],:], label= str(Angle[Max_angle[0][0]]) +'°')
plt.legend()
plt.title('Undoped LN type 1, Spectrum of the amplified signal for a given angle')



print('\n Type 2 phase matching')

PM2_LN = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
     CA = b*np.pi/180
     Dk = 2*np.pi*ne_LN_tet(Lpump,CA)/Lpump - 2*np.pi*ne_LN_tet(LTHz,CA)/LTHz - 2*np.pi*no_LN(Lidler)/Lidler
     PM2_LN[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
     for j, f in enumerate (THz_freq):
         if math.isnan(PM2_LN[i,j]):
             PM2_LN[i,j]=0
    

fig, ax = plt.subplots(figsize=(8,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM2_LN, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(0,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.suptitle('Undoped LN type 2, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM2_LN[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')

plt.figure()
plt.plot(THz_freq,PM2_LN[Max_angle[0][0],:], label= str(Angle[Max_angle[0][0]]) +'°')
plt.legend()
plt.title('Undoped LN type 2, Spectrum of the amplified signal for a given angle')




















#%%%%Lithium Niobate LN  5% MgO doped   negative unixial
## PROBLEM : only up to 5µm
print ('\n %%%%%%%%%%%%%%%%%%%%%%%%% \n 5% MgO doped Lithium Niobate ')

Ae =2.4272
Be =0.01478
Ce = 1.4617
De = 0.05612
Ee = 9.6536
Fe = 371.216

Ao =  2.2454
Bo =  0.01242
Co = 1.3005
Do = 0.05313
Eo = 6.8972
Fo = 331.33
 

def no_LNd(l):
    return np.sqrt(Ae*l**2/(l**2 -Be)+Ce*l**2/(l**2 -De) + Ee*l**2/(l**2 -Fe) + 1)

def ne_LNd(l):
    return np.sqrt(Ao*l**2/(l**2 -Bo)+Co*l**2/(l**2 -Do) + Eo*l**2/(l**2 -Fo) + 1)

def ne_LNd_tet (a,l):
    return np.sqrt( 1/(np.cos(a)**2/no_LNd(l)**2  + np.sin(a)**2/ne_LNd(l)**2))


c=0.299792458   #µm/fs
C_thick=500 #µm
freq_min=50
freq_max=260
THz_freq=np.arange(freq_min,freq_max,0.04)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


#%Type 1 phase matching for LN
print('Type 1 phase matching')


# plt.figure()
# plt.plot(LTHz, ne_LNd(LTHz), color='red', label='ne')
# plt.plot(LTHz, no_LNd(LTHz), color='blue', label='no')
# plt.legend()
# plt.title('5% Doped LN crystal refractive indices')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM1_LNd = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
#    CA=b*np.pi/180
  #   print(a)
      CA = b*np.pi/180
      Dk = 2*np.pi*ne_LNd_tet(Lpump,CA)/Lpump - 2*np.pi*no_LNd(Lidler)/Lidler - 2*np.pi*no_LNd(LTHz)/LTHz
     # Dk = 2*np.pi*no_LNd(Lpump)/Lpump - 2*np.pi*ne_LNd_tet(Lidler,CA)/Lidler - 2*np.pi*ne_LNd_tet(LTHz,CA)/LTHz
      PM1_LNd[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM1_LNd[i,j]):
              PM1_LNd[i,j]=0
     

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM1_LNd, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(41,48)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.suptitle('5% MgO doped LN type 1, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)


Profile_THzd=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THzd[i] = np.sum(PM1_LNd[i,:])
Max_angle=np.where(Profile_THzd==np.max(Profile_THzd))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(4,4),constrained_layout=True)
ax.plot(THz_freq,PM1_LNd[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'5% MgO doped LN type 1, \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))


fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))









#%Type 2 phase matching for doped LN
print('\n Type 2 phase matching')

PM2_LNd = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
      CA = b*np.pi/180
      Dk = 2*np.pi*ne_LNd_tet(Lpump,CA)/Lpump - 2*np.pi*ne_LNd_tet(Lidler, CA)/Lidler - 2*np.pi*no_LNd(LTHz)/LTHz
    #  Dk = 2*np.pi*ne_LNd_tet(Lpump,CA)/Lpump - 2*np.pi*ne_LNd_tet(LTHz, CA)/LTHz - 2*np.pi*no_LNd(Lidler)/Lidler
      PM2_LNd[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM2_LNd[i,j]):
              PM2_LNd[i,j]=0
     

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM2_LNd, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(10,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.suptitle('5% MgO doped LN type 2, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)


Profile_THzd=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THzd[i] = np.sum(PM2_LNd[i,:])
Max_angle=np.where(Profile_THzd==np.max(Profile_THzd))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(4,4),constrained_layout=True)
ax.plot(THz_freq,PM2_LNd[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'5% MgO doped LN type 2, \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))


fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))



























#%% KTA Type 1, PLAN XZ
### refractive indices: Journal of the Optical Society of America B Vol. 12, Issue 5, pp. 794-796 (1995) •https://doi.org/10.1364/JOSAB.12.000794
print('\n KTA crystal')

def nx_KTA(l):
    return np.sqrt( 1.90713 + 1.23522*l**2 /(l**2 - 0.19692**2)-0.01025*l**2)
def no_KTA(l):
    return np.sqrt( 2.15912 + 1.00099*l**2 /(l**2 - 0.21844**2)-0.01096*l**2)
def nz_KTA(l):
    return np.sqrt(2.14768 + 1.29559*l**2/(l**2 - 0.22719**2)-0.01436*l**2)


def net_KTA(l,a):
    return nx_KTA(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (nx_KTA(l)*np.tan(a)/nz_KTA(l))**2 )
                             )
Angle = np.arange(0,90,0.1)
c=0.299792458   #µm/fs
C_thick=20000 #µm
freq_min=50
freq_max=260
THz_freq=np.arange(freq_min,freq_max,0.02)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


#%Type 1 phase matching for KTA
print('Type 1 phase matching: pump (o) --> idler (e) + signal (e)')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM1_KTA = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
      CA = b*np.pi/180
      Dk = 2*np.pi*no_KTA(Lpump)/Lpump - 2*np.pi*net_KTA(Lidler, CA)/Lidler - 2*np.pi*net_KTA(LTHz, CA)/LTHz
      PM1_KTA[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM1_KTA[i,j]):
              PM1_KTA[i,j]=0

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM1_KTA, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(32,38)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.text(0,1.05,'KTA type 1: pump (o) --> idler (e) + signal (e) \n pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching.png', dpi=300, bbox_inches='tight')



Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM1_KTA[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(4,4),constrained_layout=True)
ax.plot(THz_freq,PM1_KTA[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'KTA type 1: pump (o) --> idler (e) + signal (e) \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_2.png', dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_3.png', dpi=300, bbox_inches='tight')



## On cherche maintenant pour une longueur d'onde de signal précise
## Ici, on s'intéresse à 3.4µm
## Quel est l'angle optimum pour 3.4µm ?

Target = round(lam2freq(3.4),1)  #target frequency in THz
idx = np.where(np.isclose(THz_freq, Target))[0][0]
# print(idx)
# print(THz_freq[idx])
# print(Target)
idxmax = np.where(PM1_KTA[:, idx] == np.max(PM1_KTA[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM1_KTA[:, idx], color='red')
ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
fig.text(0,1.1,'Angular profile at 3.4 µm  ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_4.png', dpi=300, bbox_inches='tight')







fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(THz_freq, PM1_KTA[idxmax, :], color='red')
# ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,2,3,4,5,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_5.png', dpi=300, bbox_inches='tight')



FWHM_freq = round(FWHM(THz_freq, PM1_KTA[idxmax, :]),0)
print(str(FWHM_freq) + ' THz FWHM' )
#assuming gaussian shape
print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')



#%% KTA Type 2, option 2  PLAN XZ
print('\n KTA crystal')

def nx_KTA(l):
    return np.sqrt( 1.90713 + 1.23522*l**2 /(l**2 - 0.19692**2)-0.01025*l**2)
def no_KTA(l):
    return np.sqrt( 2.15912 + 1.00099*l**2 /(l**2 - 0.21844**2)-0.01096*l**2)
def nz_KTA(l):
    return np.sqrt(2.14768 + 1.29559*l**2/(l**2 - 0.22719**2)-0.01436*l**2)


def net_KTA(l,a):
    return nx_KTA(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (nx_KTA(l)*np.tan(a)/nz_KTA(l))**2 )
                             )
Angle = np.arange(0,90,0.1)
c=0.299792458   #µm/fs
C_thick=500 #µm
freq_min=50
freq_max=260
THz_freq=np.arange(freq_min,freq_max,0.2)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


#%Type 2 phase matching for KTA
print('Type 2 phase matching: pump (o) --> idler (o) + signal (e)')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM2_KTA = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
      CA = b*np.pi/180
      Dk = 2*np.pi*no_KTA(Lpump)/Lpump - 2*np.pi*no_KTA(Lidler)/Lidler - 2*np.pi*net_KTA(LTHz, CA)/LTHz
      PM2_KTA[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM2_KTA[i,j]):
              PM2_KTA[i,j]=0

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM2_KTA, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(0,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.text(0,1.05,'KTA type 2: pump (o) --> idler (o) + signal (e) \n pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching.png', dpi=300, bbox_inches='tight')


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM2_KTA[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(4,4),constrained_layout=True)
ax.plot(THz_freq,PM2_KTA[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'KTA type 2: pump (o) --> idler (o) + signal (e) \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_2.png', dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_3.png', dpi=300, bbox_inches='tight')



## On cherche maintenant pour une longueur d'onde de signal précise
## Ici, on s'intéresse à 3.4µm
## Quel est l'angle optimum pour 3.4µm ?

Target = round(lam2freq(3.4),1)  #target frequency in THz
idx = np.where(np.isclose(THz_freq, Target))[0][0]
# print(idx)
# print(THz_freq[idx])
# print(Target)
idxmax = np.where(PM2_KTA[:, idx] == np.max(PM2_KTA[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM2_KTA[:, idx], color='red')
#ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
fig.text(0,1.1,'Angular profile at 3.4 µm  ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_4.png', dpi=300, bbox_inches='tight')







fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(THz_freq, PM2_KTA[idxmax, :], color='red')
# ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,2,3,4,5,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_5.png', dpi=300, bbox_inches='tight')



FWHM_freq = round(FWHM(THz_freq, PM2_KTA[idxmax, :]),0)
print(str(FWHM_freq) + ' THz FWHM' )
#assuming gaussian shape
print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')








#%% KTA Type 2, option 1  PLAN XZ
print('\n KTA crystal')

def nx_KTA(l):
    return np.sqrt( 1.90713 + 1.23522*l**2 /(l**2 - 0.19692**2)-0.01025*l**2)
def no_KTA(l):
    return np.sqrt( 2.15912 + 1.00099*l**2 /(l**2 - 0.21844**2)-0.01096*l**2)
def nz_KTA(l):
    return np.sqrt(2.14768 + 1.29559*l**2/(l**2 - 0.22719**2)-0.01436*l**2)


def net_KTA(l,a):
    return nx_KTA(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (nx_KTA(l)*np.tan(a)/nz_KTA(l))**2 )
                             )
Angle = np.arange(37,45,0.01)
c=0.299792458   #µm/fs
C_thick=2000 #µm
freq_min=55
freq_max=125
THz_freq=np.arange(freq_min,freq_max,0.01)
#THz_freq=np.arange(50,150,0.01)
Lpump=0.8 #µm 1.03
LTHz = 1/(THz_freq/c)*1000


#%Type 2 phase matching for KTA
print('Type 2 phase matching: 1.03µm pump (o) --> idler (e) + 3µm signal (o)')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM2_KTA = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
      CA = b*np.pi/180
      Dk = 2*np.pi*no_KTA(Lpump)/Lpump - 2*np.pi*net_KTA(Lidler, CA)/Lidler - 2*np.pi*no_KTA(LTHz)/LTHz
      PM2_KTA[i,:]=(np.sinc(Dk*C_thick/2/np.pi))**2  #(np.sinc(Dk*C_thick/2/np.pi))**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM2_KTA[i,j]):
              PM2_KTA[i,j]=0

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM2_KTA, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
#ax.set_ylim(0,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.text(0,1.05,'KTA type 2: pump (o) --> idler (o) + signal (e) \n pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching.png', dpi=300, bbox_inches='tight')


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM2_KTA[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(4,4),constrained_layout=True)
ax.plot(THz_freq,PM2_KTA[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'KTA type 2: pump (o) --> idler (o) + signal (e) \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_2.png', dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_3.png', dpi=300, bbox_inches='tight')



## On cherche maintenant pour une longueur d'onde de signal précise
## Ici, on s'intéresse à 3.4µm
## Quel est l'angle optimum pour 3.4µm ?
Target = round(lam2freq(3.86),1)  #target frequency in THz
idx = np.where(np.isclose(THz_freq, Target))[0][0]
# print(idx)
# print(THz_freq[idx])
# print(Target)
idxmax = np.where(PM2_KTA[:, idx] == np.max(PM2_KTA[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM2_KTA[:, idx], color='red')
ax.set_xlim(Angle[idxmax]-0.5, Angle[idxmax]+0.5)
fig.text(0,1.1,'Angular profile at 3.4 µm  ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_4.png', dpi=300, bbox_inches='tight')




# #FWHM function
# def FWHM (x, y):
#     limit = 80
#     ymax = np.max(y[limit:-1])
#     print('\n', ymax*0.5 , 'max*0.5')
#     for i, value in enumerate (y):
#         if value==ymax:
#             index_max=i
#             print(index_max)
#     for i, value in enumerate ( y[0:index_max] ):
#         if i>limit:
#             if value >= 0.5*y[index_max]:
#                 low = i
#                 break
#     print(low)
#     print(y[low] , 'y[low]')
#     for i, value in enumerate (y):
#         if i > index_max+1:
#             if value <= 0.5*y[index_max]:
#                 high = i
#                 break
#     print(high)
#     print(y[high] , 'y[high]', '\n')
#     return (x[high] - x[low])





fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(THz_freq, PM2_KTA[idxmax, :], color='red')
#ax.set_xlim(THz_freq[idx-500], THz_freq[idx+500])
fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((2.9,2.95,3,3.05,3.1))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_5.png', dpi=300, bbox_inches='tight')



FWHM_freq = round(FWHM(THz_freq, PM2_KTA[idxmax, :]),1)
print(str(FWHM_freq) + ' THz FWHM' )
#assuming gaussian shape
print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')





#%% KTA Type 2, option 2  PLAN YZ
print('\n KTA crystal')

def no_KTA(l):
    return np.sqrt( 1.90713 + 1.23522*l**2 /(l**2 - 0.19692**2)-0.01025*l**2)
def ny_KTA(l):
    return np.sqrt( 2.15912 + 1.00099*l**2 /(l**2 - 0.21844**2)-0.01096*l**2)
def nz_KTA(l):
    return np.sqrt(2.14768 + 1.29559*l**2/(l**2 - 0.22719**2)-0.01436*l**2)


def net_KTA(l,a):
    return ny_KTA(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (ny_KTA(l)*np.tan(a)/nz_KTA(l))**2 )
                             )
Angle = np.arange(0,90,0.1)
c=0.299792458   #µm/fs
C_thick=2000 #µm
freq_min=50
freq_max=260
THz_freq=np.arange(freq_min,freq_max,0.2)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


#%Type 2 phase matching for KTA
print('Type 2 phase matching: pump (o) --> idler (o) + signal (e)')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM2_KTA = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
      CA = b*np.pi/180
      Dk = 2*np.pi*no_KTA(Lpump)/Lpump - 2*np.pi*no_KTA(Lidler)/Lidler - 2*np.pi*net_KTA(LTHz, CA)/LTHz
      PM2_KTA[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM2_KTA[i,j]):
              PM2_KTA[i,j]=0

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM2_KTA, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(0,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.text(0,1.05,'KTA type 2: pump (o) --> idler (o) + signal (e) \n pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching.png', dpi=300, bbox_inches='tight')


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM2_KTA[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(4,4),constrained_layout=True)
ax.plot(THz_freq,PM2_KTA[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'KTA type 2: pump (o) --> idler (o) + signal (e) \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_2.png', dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_3.png', dpi=300, bbox_inches='tight')



## On cherche maintenant pour une longueur d'onde de signal précise
## Ici, on s'intéresse à 3.4µm
## Quel est l'angle optimum pour 3.4µm ?

Target = round(lam2freq(3.4),1)  #target frequency in THz
idx = np.where(np.isclose(THz_freq, Target))[0][0]
# print(idx)
# print(THz_freq[idx])
# print(Target)
idxmax = np.where(PM2_KTA[:, idx] == np.max(PM2_KTA[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM2_KTA[:, idx], color='red')
#ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
fig.text(0,1.1,'Angular profile at 3.4 µm  ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_4.png', dpi=300, bbox_inches='tight')







fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(THz_freq, PM2_KTA[idxmax, :], color='red')
# ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,2,3,4,5,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_5.png', dpi=300, bbox_inches='tight')



FWHM_freq = round(FWHM(THz_freq, PM2_KTA[idxmax, :]),1)
print(str(FWHM_freq) + ' THz FWHM' )
#assuming gaussian shape
print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')


#%% KTA Type 2, option 1 Plane YZ
print('\n KTA crystal')

def no_KTA(l):
    return np.sqrt( 1.90713 + 1.23522*l**2 /(l**2 - 0.19692**2)-0.01025*l**2)
def ny_KTA(l):
    return np.sqrt( 2.15912 + 1.00099*l**2 /(l**2 - 0.21844**2)-0.01096*l**2)
def nz_KTA(l):
    return np.sqrt(2.14768 + 1.29559*l**2/(l**2 - 0.22719**2)-0.01436*l**2)


def net_KTA(l,a):
    return ny_KTA(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (ny_KTA(l)*np.tan(a)/nz_KTA(l))**2 )
                             )
Angle = np.arange(0,90,0.1)
c=0.299792458   #µm/fs
C_thick=7000 #µm
freq_min=10
freq_max=150
THz_freq=np.arange(freq_min,freq_max,0.01)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


#%Type 2 phase matching for KTA
print('Type 2 phase matching: pump (o) --> idler (e) + signal (o)')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM2_KTA = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
      CA = b*np.pi/180
      Dk = 2*np.pi*no_KTA(Lpump)/Lpump - 2*np.pi*net_KTA(Lidler,CA)/Lidler - 2*np.pi*no_KTA(LTHz)/LTHz
      PM2_KTA[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM2_KTA[i,j]):
              PM2_KTA[i,j]=0

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM2_KTA, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.set_ylim(0,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.text(0,1.05,'KTA type 2: pump (o) --> idler (o) + signal (e) \n pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching.png', dpi=300, bbox_inches='tight')


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM2_KTA[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
ax.plot(THz_freq,PM2_KTA[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'KTA type 2: pump (o) --> idler (o) + signal (e) \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_2.png', dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_3.png', dpi=300, bbox_inches='tight')



## On cherche maintenant pour une longueur d'onde de signal précise
## Ici, on s'intéresse à 3.4µm
## Quel est l'angle optimum pour 3.4µm ?


Target = round(lam2freq(3.4),1)  #target frequency in THz
idx = np.where(np.isclose(THz_freq, Target))[0][0]
# print(idx)
# print(THz_freq[idx])
# print(Target)
idxmax = np.where(PM2_KTA[:, idx] == np.max(PM2_KTA[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM2_KTA[:, idx], color='red')
#ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
fig.text(0,1.1,'Angular profile at 3.4 µm  ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_4.png', dpi=300, bbox_inches='tight')







fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(THz_freq, PM2_KTA[idxmax, :], color='red')
# ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,2,3,4,5,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_5.png', dpi=300, bbox_inches='tight')


#FWHM function
def FWHM (x, y):
    limit = 6000
    ymax = np.max(y[limit:-1])
    print('\n', ymax*0.5 , 'max*0.5')
    for i, value in enumerate (y):
        if value==ymax:
            index_max=i
            print(index_max)
    for i, value in enumerate ( y[0:index_max] ):
        if i>limit:
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


# #FWHM function
# def FWHM (x, y):
#     ymax = np.max(y)
#     print('\n', ymax*0.5 , 'max*0.5')
#     for i, value in enumerate (y):
#         if value==ymax:
#             index_max=i
#             print(index_max)
#     for i, value in enumerate ( y[0:index_max] ):
#         if value >= 0.5*y[index_max]:
#             low = i
#             break
#     print(low)
#     print(y[low] , 'y[low]')
#     for i, value in enumerate (y):
#         if i > index_max+1:
#             if value <= 0.5*y[index_max]:
#                 high = i
#                 break
#     print(high)
#     print(y[high] , 'y[high]', '\n')
#     return (x[high] - x[low])


FWHM_freq = round(FWHM(THz_freq, PM2_KTA[idxmax, :]),1)
print(str(FWHM_freq) + ' THz FWHM' )
#assuming gaussian shape
print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')



#%% KTA Type 1   Plane YZ
print('\n KTA crystal')

def no_KTA(l):
    return np.sqrt( 1.90713 + 1.23522*l**2 /(l**2 - 0.19692**2)-0.01025*l**2)
def ny_KTA(l):
    return np.sqrt( 2.15912 + 1.00099*l**2 /(l**2 - 0.21844**2)-0.01096*l**2)
def nz_KTA(l):
    return np.sqrt(2.14768 + 1.29559*l**2/(l**2 - 0.22719**2)-0.01436*l**2)


def net_KTA(l,a):
    return ny_KTA(l)*np.sqrt(  (1+(np.tan(a))**2)/(1+ (ny_KTA(l)*np.tan(a)/nz_KTA(l))**2 )
                             )
Angle = np.arange(26,30,0.005)
c=0.299792458   #µm/fs
C_thick=10000 #µm
freq_min=50
freq_max=150
THz_freq=np.arange(freq_min,freq_max,0.01)
Lpump=1.03 #µm
LTHz = 1/(THz_freq/c)*1000


#%Type 2 phase matching for KTA
print('Type 1 phase matching: pump (o) --> idler (e) + signal (e)')


Pump_freq=c/Lpump*1000
idler_freq=Pump_freq-THz_freq
Lidler=1/(idler_freq/c)*1000
PM2_KTA = np.zeros((np.size(Angle), np.size(LTHz)))

for i, b in enumerate (Angle):
      CA = b*np.pi/180
      Dk = 2*np.pi*no_KTA(Lpump)/Lpump - 2*np.pi*net_KTA(Lidler,CA)/Lidler - 2*np.pi*net_KTA(LTHz,CA)/LTHz
      PM2_KTA[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
      for j, f in enumerate (THz_freq):
          if math.isnan(PM2_KTA[i,j]):
              PM2_KTA[i,j]=0

fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
im = ax.pcolor(THz_freq, Angle, PM2_KTA, cmap='jet')
ax.set_ylabel('angle (deg)', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
#ax.set_ylim(0,90)  
ax.tick_params(labelsize=s)
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,2.5,3,4,5,6))
secax.tick_params(labelsize=s)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
fig.text(0,1.05,'KTA type 2: pump (o) --> idler (o) + signal (e) \n pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching.png', dpi=300, bbox_inches='tight')


Profile_THz=np.zeros(np.size(Angle))
for i, a in enumerate (Angle):
    Profile_THz[i] = np.sum(PM2_KTA[i,:])
Max_angle=np.where(Profile_THz==np.max(Profile_THz))
print('index ' + str(Max_angle[0][0]))
print('angle ' + str(Angle[Max_angle[0][0]]) + '°')


fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
ax.plot(THz_freq,PM2_KTA[Max_angle[0][0],:], color='navy', label= str(Angle[Max_angle[0][0]]) +'°')
ax.legend()
ax.set_ylabel('Relative intensity', fontsize=s) 
ax.set_xlabel('signal frequency (THz)', fontsize=s) 
ax.text(1,1.3,'KTA type 2: pump (o) --> idler (o) + signal (e) \n Spectrum of the amplified signal for the optimum angle')
secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,1.2,1.5,2,3,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_2.png', dpi=300, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(LTHz, Lidler, color='green')
ax.set_ylabel('idler wavelength', fontsize=s) 
ax.set_xlabel('signal wavelength', fontsize=s) 
ax.set_xticks((1,2,3,4,5,6))
ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_3.png', dpi=300, bbox_inches='tight')



## On cherche maintenant pour une longueur d'onde de signal précise
## Ici, on s'intéresse à 3.4µm
## Quel est l'angle optimum pour 3.4µm ?


Target = round(lam2freq(3.4),1)  #target frequency in THz
idx = np.where(np.isclose(THz_freq, Target))[0][0]
# print(idx)
# print(THz_freq[idx])
# print(Target)
idxmax = np.where(PM2_KTA[:, idx] == np.max(PM2_KTA[:, idx]) )[0][0]
fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(Angle, PM2_KTA[:, idx], color='red')
#ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
fig.text(0,1.1,'Angular profile at 3.4 µm  ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
ax.set_xlabel('Angle (deg)')
ax.set_ylabel('Output signal')
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_4.png', dpi=300, bbox_inches='tight')







fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
ax.plot(THz_freq, PM2_KTA[idxmax, :], color='red')
# ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
ax.set_xlabel('Signal frequency (THz)')
ax.set_ylabel('Output signal')
secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
secax.set_xlabel('signal wavelength (µm)', fontsize=s)
secax.set_xticks((1,2,3,4,5,6))
fig.savefig(r'D:\phase matching OPA\python\results\KTA_phase_matching_5.png', dpi=300, bbox_inches='tight')


# #FWHM function
# def FWHM (x, y):
#     limit = 1200
#     ymax = np.max(y[0:limit])
#     print('\n', ymax*0.5 , 'max*0.5')
#     for i, value in enumerate (y):
#         if value==ymax:
#             index_max=i
#             print(index_max)
#     for i, value in enumerate ( y[0:index_max] ):
#         if value >= 0.5*y[index_max]:
#             low = i
#             break
#     print(low)
#     print(y[low] , 'y[low]')
#     for i, value in enumerate (y):
#         if i > index_max+1:
#             if i<limit:
#                 if value <= 0.5*y[index_max]:
#                     high = i
#                     break
#     print(high)
#     print(y[high] , 'y[high]', '\n')
#     return (x[high] - x[low])


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


FWHM_freq = round(FWHM(THz_freq, PM2_KTA[idxmax, :]),1)
print(str(FWHM_freq) + ' THz FWHM' )
#assuming gaussian shape
print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')






#%% CASE 3 : LGS crystal LiGaS2
#biaxial crystal
#non-collinear geometry. We denote by alpha the angle between the pump and the seed.



# #%% GaSe
# # o and e
# c = 299792458
# Lis = np.arange(c/(100*1e6),c/(15*1e6),0.01) #10 #idler wavelength (µm)        ordinary
# result1 = np.zeros(np.size(Lis))
# angles=np.arange(1,90,0.1)
# I3=np.zeros((np.size(Lis), np.size(angles)) )


# for n, Li in enumerate (Lis):
#     #for m, Beta in enumerate (angles):
        
#      #print('\n GaSe crystal')
#      Lp = 1.95 #pump wavelength (µm)        extraordinary
#      Ls = (1/Lp - 1/Li)**(-1) #signal wavelength (µm)          ordinary
#      #print('Lp=' + str(Lp) +'µm, Ls=' + str(round(Ls,4)) + 'µm, Li=' + str(Li) + 'µm')
#      def no(l):
#          return np.sqrt(10.6409 + 0.3788/((l**2) - 0.1232) + 7090.7/(l**2 -2216.3))
#      def ne(l):
#          return np.sqrt(8.2477 + 0.2881/((l**2) - 0.1669) +4927.5 /((l)**2-1990.1))
#      def net(l,tet):
#          return np.sqrt(1/( np.cos(tet)**2/no(l)**2 + np.sin(tet)**2/ne(l)**2 ))
#      #Type 1 phase matching : e (pump, 0.8µm) --> o + o 
#      def y1(x) :
#          return net(Lp,x)/Lp - no(Ls)/Ls - no(Li)/Li 
#     # def y2(x) :
#         # return net(Lp,Beta)/Lp - no(x)/x - no(Li)/Li 
#      # plt.figure()
#      # plt.plot(np.rad2deg(np.arange(-3.14/2,3.14/2,0.001)), y1(np.arange(-3.14/2,3.14/2,0.001)))
#      result1[n] = scipy.optimize.fsolve(y1,0.6)
#     # Ls = scipy.optimize.fsolve(y2,0.6)
     
#      Dk = net(Lp,np.deg2rad( angles))*2*np.pi/(Lp*1e-6)  - no(Li)*2*np.pi/(Li*1e-6) - no(Ls)*2*np.pi/(Ls*1e-6)
#      I3[n,:]=(np.sinc(Dk*2e-3/2))**2
     
        
        
        
        
    
#     #print(result)
# print('Type 1 phase matching angle: ' + str(np.round(np.rad2deg(result1[0]),1)) + '°')
# plt.figure(figsize=(4,7))
# plt.title('GaSe Type 1')
# plt.plot(c/Lis*1e-6,np.rad2deg(result1))
# plt.xlabel('Frequency of the idler (THz)')
# plt.ylim(0,90)
# plt.ylabel('Phase matching angle (deg)')


# # plt.figure()
# # plt.pcolor(Lis, angles, I3)

# plt.figure(figsize=(4,6))
# # plt.pcolor(c/Lis*1e-6,angles,np.transpose(np.log(I3)))
# # plt.figure()
# plt.ylim(0,90)
# plt.pcolor(c/Lis*1e-6,angles,np.transpose((I3)))

# #I need the generation efficiency when it's not perfectly phase matched

# #sinc²(Delta k L/2)

# # #%%
# # plt.figure()
# # xx = np.arange(-10,10,0.1)
# # plt.plot(xx, np.sinc(xx)**2)


# #%% AGSe, 1.55µm

# # o and e
# c = 299792458
# Lis = np.arange(c/(50*1e6),c/(15*1e6),0.01) #10 #idler wavelength (µm)        ordinary
# result1 = np.zeros(np.size(Lis))
# angles=np.arange(1,90,0.1)
# I3=np.zeros((np.size(Lis), np.size(angles)) )


# for n, Li in enumerate (Lis):
#     #for m, Beta in enumerate (angles):
        
#      #print('\n GaSe crystal')
#      Lp = 1.55 #pump wavelength (µm)        extraordinary
#      Ls = (1/Lp - 1/Li)**(-1) #signal wavelength (µm)          ordinary
#      #print('Lp=' + str(Lp) +'µm, Ls=' + str(round(Ls,4)) + 'µm, Li=' + str(Li) + 'µm')
#      def no(l):
#          return np.sqrt(6.8507 + 0.4297/((l**2) - 0.1584) -0.00125*l**2)
#      def ne(l):
#          return np.sqrt(6.6792 + 0.4598/((l**2) - 0.2122) -0.00126*l**2)
#      def net(l,tet):
#          return np.sqrt(1/( np.cos(tet)**2/no(l)**2 + np.sin(tet)**2/ne(l)**2 ))
#      #Type 1 phase matching : e (pump, 0.8µm) --> o + o 
#      def y1(x) :
#          return net(Lp,x)/Lp - no(Ls)/Ls - no(Li)/Li 
#     # def y2(x) :
#         # return net(Lp,Beta)/Lp - no(x)/x - no(Li)/Li 
#      # plt.figure()
#      # plt.plot(np.rad2deg(np.arange(-3.14/2,3.14/2,0.001)), y1(np.arange(-3.14/2,3.14/2,0.001)))
#      result1[n] = scipy.optimize.fsolve(y1,0.6)
#     # Ls = scipy.optimize.fsolve(y2,0.6)
     
#      Dk = net(Lp,np.deg2rad( angles))*2*np.pi/(Lp*1e-6)  - no(Li)*2*np.pi/(Li*1e-6) - no(Ls)*2*np.pi/(Ls*1e-6)
#      I3[n,:]=(np.sinc(Dk*0.5e-3/(2*np.pi)))**2
     
     
    
#     #print(result)
# print('Type 1 phase matching angle: ' + str(np.round(np.rad2deg(result1[0]),1)) + '°')
# plt.figure(figsize=(4,7))
# plt.title('GaSe Type 1')
# plt.plot(c/Lis*1e-6,np.rad2deg(result1))
# plt.xlabel('Frequency of the idler (THz)')
# plt.ylim(0,90)
# plt.ylabel('Phase matching angle (deg)')

# #%%

# plt.figure(figsize=(5,5.5))
# plt.ylim(0,90)
# plt.pcolor(c/Lis*1e-6,angles,np.transpose((I3)), cmap='viridis')
# plt.colorbar()


# #%%




























#%% 
"""objective here: find the best crystal for SHG of 3400 nm
find the best thickness too
(spectral bandwidth <--> thickness)
"""

#%%  SHG for 3400 nm --> 1700 nm
##  AGS
print('AgGaS2/AGS, SHG @3400nm')
## o + o --> e
#refractive index of AGS depending on the angle (refractiveindex.info, AgGaS2, 293.15 Kelvin) OK !!
def no_AGS(l):
    return np.sqrt(5.7975 + 0.2311/((l**2)-0.0688)-0.00257*(l)**2)
def ne_AGS(l):
    return np.sqrt(5.5436 + 0.2230/((l**2)-0.0946)-0.00261*(l)**2)
def net_AGS(l,tet): 
    return np.sqrt(1/(np.cos(tet)**2/no_AGS(l)**2 + np.sin(tet)**2/ne_AGS(l)**2 ))


Angle_FWHM=[]
Duration_FWHM=[]

C_thicks = [500,1000,2000,3000,4000,5000,6000,7000,8000]
for it, C_thick in enumerate (C_thicks):
    
    Angle = np.arange(20,50,0.01)
    c=0.299792458   #µm/fs
    #C_thick=8000 #µm
    freq_min=55#freq2lam(3.4+1)
    freq_max=120#freq2lam(3.4-1)
    THz_freq=np.arange(freq_min,freq_max,0.01)      #
    Lpump=3.4 #µm   
    LTHz = 1/(THz_freq/c)*1000
    
    ## frequency (THz) range for the FWHM calculation:
    r1 =56
    r2=110
    
    # angle (deg) range for the FWHM calculation:
    ra1= 25
    ra2= 40    
    
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
    
    
    
    
    #%Type 1 phase matching for LN
    print('o + o --> e')
    
    Pump_freq=c/Lpump*1000
    idler_freq=Pump_freq+THz_freq
    Lidler=1/(idler_freq/c)*1000
    PM_AGS = np.zeros((np.size(Angle), np.size(LTHz)))
    
    for i, b in enumerate (Angle):
          CA = b*np.pi/180
          Dk = 2*np.pi*no_AGS(Lpump)/Lpump - 2*np.pi*net_AGS(Lidler,CA)/Lidler + 2*np.pi*no_AGS(LTHz)/LTHz
          #Dk = 2*np.pi*net_GS(Lpump,CA)/Lpump - 2*np.pi*ne_GS(Lidler)/Lidler - 2*np.pi*net_GS(LTHz,CA)/LTHz
          PM_AGS[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
          for j, f in enumerate (THz_freq):
              if math.isnan(PM_AGS[i,j]):
                  PM_AGS[i,j]=0
         
    ##### COLOR PLOT (optional)
    # fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
    # im = ax.pcolor(THz_freq, Angle, PM_AGS, cmap='jet')
    # ax.set_ylabel('angle (deg)', fontsize=s) 
    # ax.set_xlabel('signal frequency (THz)', fontsize=s)   
    # ax.tick_params(labelsize=s)
    # secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    # secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    # secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    # secax.tick_params(labelsize=s)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # fig.colorbar(im, cax=cax, orientation='vertical')
    # fig.suptitle('AGS, SHG, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    # fig.savefig(r'D:\phase matching OPA\python\results\AGS_SHG_3400nm_'+ str(C_thick)+'_plot1.png', dpi=300, bbox_inches='tight')
    
    
    # spectrum at the angle which maximizes the total signal (over the full frequency range)
    Profile_THz=np.zeros(np.size(Angle))
    for i, a in enumerate (Angle):
        Profile_THz[i] = np.sum(PM_AGS[i,:])
    Max_angle=np.where(Profile_THz==np.max(Profile_THz))
    print('index ' + str(Max_angle[0][0]))
    print('angle ' + str(Angle[Max_angle[0][0]]) + '°')
    
    
    fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
    ax.plot(THz_freq,PM_AGS[Max_angle[0][0],:], color='navy', label= str(round(Angle[Max_angle[0][0]],2)) +'°')
    ax.legend()
    ax.set_ylabel('Relative intensity', fontsize=s) 
    ax.set_xlabel('signal frequency (THz)', fontsize=s) 
    fig.text(0,1.1,'AGS SHG: pump 3.4µm (o) + signal (o) --> idler (e) \n Spectrum of the amplified signal for the optimum angle \n which gives the highest integrated signal')
    secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    fig.savefig(r'D:\phase matching OPA\python\results\AGS_SHG_3400nm_'+ str(C_thick)+'_plot2.png', dpi=300, bbox_inches='tight')
    
    
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
    #ax.set_xticks((1,2,3,4,5,6))
    #ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
    fig.text(0,1.1,'AGS SHG: pump 3.4µm (o) + signal (o) --> idler (e) ')
    fig.savefig(r'D:\phase matching OPA\python\results\AGS_SHG_3400nm_'+ str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    ## Ici, on s'intéresse à 3.4µm
    ## Quel est l'angle optimum pour 3.4µm ?
    
    Target = round(lam2freq(3.4),1)  #target frequency in THz
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
    idxmax = np.where(PM_AGS[:, idx] == np.max(PM_AGS[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM_AGS[:, idx], color='red')
    #ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
    fig.text(0,1.2,'AGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 3.4µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM_AGS[:, idx], ra1,ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Output signal')
    fig.savefig(r'D:\phase matching OPA\python\results\AGS_SHG_3400nm_'+ str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
    
    
    
    
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM_AGS[idxmax, :], color='red')
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    # ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
    fig.text(0,1.2,'AGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
    ax.set_xlabel('Signal frequency (THz)')
    ax.set_ylabel('Output signal')
    secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    
    
    
    FWHM_freq = round(FWHM(THz_freq, PM_AGS[idxmax, :],r1,r2),4)
    print(str(FWHM_freq) + ' THz FWHM' )
    fig.text(0,1.4,str(FWHM_freq) + ' THz FWHM', fontsize=s)
    
    #assuming gaussian shape
    print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')
    fig.text(0,1.3,'Assuming Gaussian shape: ' + str("{:10.4}".format(2*np.log(2)/(np.pi*FWHM_freq*1e12))) + ' s FWHM', fontsize=s)
    
    fig.text(0,1.5,'Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    
    fig.savefig(r'D:\phase matching OPA\python\results\AGS_SHG_3400nm_'+ str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
    
    
    
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM(Angle, PM_AGS[:, idx], ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    
        
        
    
    
    
fig_F=plt.figure()
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
fig_F.savefig(r'D:\phase matching OPA\python\results\AGS_SHG_3400nm_plotFINAL.png', dpi=300, bbox_inches='tight')
  
    
#%%
    ###### what's the spatial walkoff?   tan rho = -1/ne dne/dtheta
def dnedt(l,t):
    return -1/2*(   (np.cos(t))**2/(no_AGS(l))**2   +    (np.sin(t))**2/(ne_AGS(l))**2    )**(-3/2) * (2*np.sin(t)*np.cos(t)*(1/(ne_AGS(l))**2 - 1/(no_AGS(l))**2 )       )
def SpWO(l,t):
    return np.arctan(-1/net_AGS(l,t)*dnedt(l,t))
print( '\n Spatial Walkoff')
print( str(SpWO(1.7 , np.deg2rad(34.4)) )  + 'rad')
print( str(np.rad2deg(SpWO(1.7 , np.deg2rad(34.4)))) + 'deg' )    
th=0.4 #mm
print ('If thickness =' + str(th) + 'mm:')
print( 'vertical shift will be (in mm): ' + str(th*np.tan(SpWO(3.4 , np.deg2rad(34.4)))  ) )





#%% SHG for 3400 nm --> 1700nm
## GaSe (called GS in SNLO) just for ref, bc absorption is too high <-- YUYA
# refractiveindex.info:
print ('GaSe/GS, SHG @3400nm')
    
Angle = np.arange(10,16,0.01)
c=0.299792458   #µm/fs
C_thick=4000 #µm
freq_min=50#freq2lam(3.4+1)
freq_max=120#freq2lam(3.4-1)
THz_freq=np.arange(freq_min,freq_max,0.01)      #
Lpump=3.4 #µm   
LTHz = 1/(THz_freq/c)*1000

## frequency (THz) range for the FWHM calculation:
r1 =60
r2=110

# angle (deg) range for the FWHM calculation:
ra1= 11
ra2= 15


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

C_thicks = [500,1000,2000,3000,4000,5000,6000,7000,8000]
for it, C_thick in enumerate (C_thicks):

    
    for i, b in enumerate (Angle):
          CA = b*np.pi/180
          Dk = 2*np.pi*no_GS(Lpump)/Lpump - 2*np.pi*net_GS(Lidler,CA)/Lidler + 2*np.pi*no_GS(LTHz)/LTHz
          #Dk = 2*np.pi*net_GS(Lpump,CA)/Lpump - 2*np.pi*ne_GS(Lidler)/Lidler - 2*np.pi*net_GS(LTHz,CA)/LTHz
          PM_GS[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
          for j, f in enumerate (THz_freq):
              if math.isnan(PM_GS[i,j]):
                  PM_GS[i,j]=0
         
    ###" OPTIONAL COLOR PLOT
    # fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
    # im = ax.pcolor(THz_freq, Angle, PM_GS, cmap='jet')
    # ax.set_ylabel('angle (deg)', fontsize=s) 
    # ax.set_xlabel('signal frequency (THz)', fontsize=s)   
    # ax.tick_params(labelsize=s)
    # secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    # secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    # secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    # secax.tick_params(labelsize=s)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # fig.colorbar(im, cax=cax, orientation='vertical')
    # fig.suptitle('GaSe (GS), SHG, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    # fig.savefig(r'D:\phase matching OPA\python\results\GaSe_SHG_3400nm_'+str(C_thick)+'_plot1.png', dpi=300, bbox_inches='tight')
    
    
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
    fig.text(0,1.1,'GaSe (GS) SHG: pump 3.4µm (o) + signal (o) --> idler (e) \n Spectrum of the amplified signal for the optimum angle \n which gives the highest integrated signal')
    secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    fig.savefig(r'D:\phase matching OPA\python\results\GaSe_SHG_3400nm_'+str(C_thick)+'_plot2.png', dpi=300, bbox_inches='tight')
    
    
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
    #ax.set_xticks((1,2,3,4,5,6))
    #ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
    fig.text(0,1.1,'GaSe (GS) SHG: pump 3.4µm (o) + signal (o) --> idler (e) ')
    fig.savefig(r'D:\phase matching OPA\python\results\GaSe_SHG_3400nm_'+str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    ## Ici, on s'intéresse à 3.4µm
    ## Quel est l'angle optimum pour 3.4µm ?
    
    Target = round(lam2freq(3.4),1)  #target frequency in THz
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
    idxmax = np.where(PM_GS[:, idx] == np.max(PM_GS[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM_GS[:, idx], color='red')
    #ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
    fig.text(0,1.2,'GaSe (GS) SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 3.4µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM_GS[:, idx], ra1,ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Output signal')
    fig.savefig(r'D:\phase matching OPA\python\results\GaSe_SHG_3400nm_'+str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
    
    
    
    
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM_GS[idxmax, :], color='red')
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    # ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
    fig.text(0,1.2,'GaSe (GS) SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
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
    
    fig.savefig(r'D:\phase matching OPA\python\results\GaSe_SHG_3400nm_'+str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
    
    
    
    
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM(Angle, PM_GS[:, idx], ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    
        
        
    
    
    
fig_F=plt.figure()
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
plt.title('GaSe for SHG 3.4µm (o) + 3.4µm (o) --> 1.7µm (e)')
fig_F.savefig(r'D:\phase matching OPA\python\results\GS_SHG_3400nm_plotFINAL.png', dpi=300, bbox_inches='tight')
  
    
    
#%%  SHG of 3400nm using LGS?
######## 3400 (o) + 3400 (o) --> 1700 (e)
############

print ('LGS, SHG @3400nm')
######## ça ne MARCHE PAS!!! LES ANGLES SONT PAS BONS
##########TO BE CHECKED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
Angle = np.arange(53,56,0.01)
c=0.299792458   #µm/fs
C_thick=4000 #µm
freq_min=20#freq2lam(3.4+1)
freq_max=160#freq2lam(3.4-1)
THz_freq=np.arange(freq_min,freq_max,0.01)      #
Lpump=3.4 #µm   
LTHz = 1/(THz_freq/c)*1000

## frequency (THz) range for the FWHM calculation:
r1 =80
r2=150

# angle (deg) range for the FWHM calculation:
ra1= 53.5
ra2= 55.5


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

#%Type 1 phase matching for LN
print('o + o --> e')

Pump_freq=c/Lpump*1000
idler_freq=Pump_freq+THz_freq
Lidler=1/(idler_freq/c)*1000
PM_LGS = np.zeros((np.size(Angle), np.size(LTHz)))




Angle_FWHM=[]
Duration_FWHM=[]

C_thicks = [10000,15000,20000,25000,30000,35000,40000]
for it, C_thick in enumerate (C_thicks):

    
    for i, b in enumerate (Angle):
          CA = b*np.pi/180
          Dk = 2*np.pi*noLGS(Lpump)/Lpump - 2*np.pi*neLGS(Lidler,CA)/Lidler + 2*np.pi*noLGS(LTHz)/LTHz
          PM_LGS[i,:]=np.sinc(Dk*C_thick/2/np.pi)**2
          for j, f in enumerate (THz_freq):
              if math.isnan(PM_LGS[i,j]):
                  PM_LGS[i,j]=0
         
    # ##" OPTIONAL COLOR PLOT
    # fig, ax = plt.subplots(figsize=(6,5),constrained_layout=True)
    # im = ax.pcolor(THz_freq, Angle, PM_LGS, cmap='jet')
    # ax.set_ylabel('angle (deg)', fontsize=s) 
    # ax.set_xlabel('signal frequency (THz)', fontsize=s)   
    # ax.tick_params(labelsize=s)
    # secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    # secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    # secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    # secax.tick_params(labelsize=s)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # fig.colorbar(im, cax=cax, orientation='vertical')
    # fig.suptitle('LGS, SHG, pump ='+ str(Lpump) + 'µm, Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    # fig.savefig(r'D:\phase matching OPA\python\results\LGS_SHG_3400nm_'+str(C_thick)+'_plot1.png', dpi=300, bbox_inches='tight')
    
    
    # spectrum at the angle which maximizes the total signal (over the full frequency range)
    Profile_THz=np.zeros(np.size(Angle))
    for i, a in enumerate (Angle):
        Profile_THz[i] = np.sum(PM_LGS[i,:])
    Max_angle=np.where(Profile_THz==np.max(Profile_THz))
    print('index ' + str(Max_angle[0][0]))
    print('angle ' + str(Angle[Max_angle[0][0]]) + '°')
    
    
    fig, ax = plt.subplots(figsize=(6,4),constrained_layout=True)
    ax.plot(THz_freq,PM_LGS[Max_angle[0][0],:], color='navy', label= str(round(Angle[Max_angle[0][0]],2)) +'°')
    ax.legend()
    ax.set_ylabel('Relative intensity', fontsize=s) 
    ax.set_xlabel('signal frequency (THz)', fontsize=s) 
    fig.text(0,1.1,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e) \n Spectrum of the amplified signal for the optimum angle \n which gives the highest integrated signal')
    secax=ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    fig.savefig(r'D:\phase matching OPA\python\results\LGS_SHG_3400nm_'+str(C_thick)+'_plot2.png', dpi=300, bbox_inches='tight')
    
    
    
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(LTHz, Lidler, color='green')
    ax.set_ylabel('idler wavelength', fontsize=s) 
    ax.set_xlabel('signal wavelength', fontsize=s) 
    #ax.set_xticks((1,2,3,4,5,6))
    #ax.set_yticks((1,2,3,4,5,6,7,8,9,10))
    fig.text(0,1.1,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e) ')
    fig.savefig(r'D:\phase matching OPA\python\results\LGS_SHG_3400nm_'+str(C_thick)+'_plot3.png', dpi=300, bbox_inches='tight')
    
    
    
    ## On cherche maintenant pour une longueur d'onde de signal précise
    ## Ici, on s'intéresse à 3.4µm
    ## Quel est l'angle optimum pour 3.4µm ?
    
    Target = round(lam2freq(3.4),1)  #target frequency in THz
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
    idxmax = np.where(PM_LGS[:, idx] == np.max(PM_LGS[:, idx]) )[0][0]
    fig, ax = plt.subplots(figsize=(3,3),constrained_layout=True)
    ax.plot(Angle, PM_LGS[:, idx], color='red')
    #ax.set_xlim(Angle[idxmax-50], Angle[idxmax+50])
    fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Angular profile at 3.4µm signal ' + r'($\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°)', fontsize=s)
    fig.text(0,1.3, str(round(np.deg2rad(FWHM(Angle, PM_LGS[:, idx], ra1,ra2))*1000,3)) + 'mrad FWHM', fontsize=s)
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Output signal')
    fig.savefig(r'D:\phase matching OPA\python\results\LGS_SHG_3400nm_'+str(C_thick)+'_plot4.png', dpi=300, bbox_inches='tight')
    
    
    
    
    
    
    ##### CHECK THE FWHM
    fig, ax = plt.subplots(figsize=(4,3),constrained_layout=True)
    ax.plot(THz_freq, PM_LGS[idxmax, :], color='red')
    ax.fill_between( [r1,r2] ,1, color='yellow', alpha=0.3)
    # ax.set_xlim(THz_freq[idx-2000], THz_freq[idx+2000])
    fig.text(0,1.2,'LGS SHG: pump 3.4µm (o) + signal (o) --> idler (e)', fontsize=s)
    fig.text(0,1.1,'Spectrum at ' + r'$\theta_{max}$='+ str(round(Angle[idxmax],1) )+'°', fontsize=s)
    ax.set_xlabel('Signal frequency (THz)')
    ax.set_ylabel('Output signal')
    secax = ax.secondary_xaxis('top', functions=(lam2freq,freq2lam)   )
    secax.set_xlabel('signal wavelength (µm)', fontsize=s)
    secax.set_xticks((2,3,4,5,6,7,8,9,10,50))
    
    
    
    FWHM_freq = round(FWHM(THz_freq, PM_LGS[idxmax, :],r1,r2),4)
    print(str(FWHM_freq) + ' THz FWHM' )
    fig.text(0,1.4,str(FWHM_freq) + ' THz FWHM', fontsize=s)
    
    #assuming gaussian shape
    print('Assuming Gaussian shape: ' + str(2*np.log(2)/(np.pi*FWHM_freq*1e12)) + ' s FWHM')
    fig.text(0,1.3,'Assuming Gaussian shape: ' + str("{:10.4}".format(2*np.log(2)/(np.pi*FWHM_freq*1e12))) + ' s FWHM', fontsize=s)
    
    fig.text(0,1.5,'Effective thickness =' + str(C_thick) + 'µm', fontsize=s)
    
    fig.savefig(r'D:\phase matching OPA\python\results\LGS_SHG_3400nm_'+str(C_thick)+'_plot5.png', dpi=300, bbox_inches='tight')
    
    
    
    
        
    Angle_FWHM=np.append(Angle_FWHM,  (round(np.deg2rad(FWHM(Angle, PM_LGS[:, idx], ra1,ra2))*1000,3)))
    Duration_FWHM=np.append(Duration_FWHM,  (2*np.log(2)/(np.pi*FWHM_freq*1e12))  )
    
        
        
    
    
    
fig_F=plt.figure()
plt.scatter(C_thicks,Angle_FWHM, color='red')
plt.ylabel('FWHM angular profile (mrad)', color='red')
plt.xlabel('crystal effective thickness (µm)')
plt.twinx()
plt.scatter(C_thicks, Duration_FWHM*1e15, marker='+', s=100, color='green')    
plt.ylabel('FWHM duration (fs)', color='green')
fig_F.savefig(r'D:\phase matching OPA\python\results\LGS_SHG_3400nm_plotFINAL.png', dpi=300, bbox_inches='tight')
  
