"""
 STRUCTURE OF THE HDF FILES :    
 The hdf5 file is divided into 3 parts : CEP values ; time stamps ; configuration
 Each carrier-enveloppe phase value is recorded as an integer : you must divide by 10 000 to get it in radians
 phase values are recorded in an array
 Each line of this array corresponds to a 20second-long measurement (parameter 'phase points' defined through the Fringeezz software : if phase points = 20 000, each line corresponds to 20sec)
 timestamps : list of times (sec) when each measurement line ended 
 printing timestamps[1] - timestamps[0] should therefore return 20 as there are 20 seconds between two lines of measurements

"""





# libraries
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import welch # PSD using Welch method



## FIRST : EXTRACT DATA 
# option 1 = ask for file location and filename
#print ('Where is the file located?')
#filepath = input('The path of the file is :')  # e.g. : C:\Users\ouille\Desktop\GitHub_users\GitHub_MarieOuille\Python-scripts\Fringeezz
#print ('What is the name of the hdf5 file?')
#filename = input('The name of the file is :')   # e.g. :  'phase_log_20170512_110754' 

#option 2 : enter it manually here :
filepath = r'Z:\Laser\CEP\20191120\fringeezz'
filename = 'phase_log_20191120_121032'



file= str(filename) + '.hdf5' 

#output directory = new directory created where "fringeezz.py" is located
path = os.getcwd()
if not os.path.exists('analyzed_data'):
    os.mkdir('analyzed_data')
output_dir=path+'\\analyzed_data'

def fringeezz_load_h5_data_log(file, filepath, print_structure=True):
    """ load data log file in HDF5 format from the fringeezz from default directory.
    input: file name, filepath if not default, print_structure=True if you want a print of the hdf5 file structure, else = False
    return: fringeezz configuration as string, timestamps array, phase matrix (rad)
    """
    with h5py.File(os.path.join(filepath,file)) as h5f:
        if print_structure:
            for name in h5f:
                print(h5f[name])
        config = np.array(h5f['configuration'], dtype = '|S1545')
        if 'phases' in h5f:
            phases = np.array(h5f['phases'], dtype = int)
        else:
            phases = float('nan')
        timestamps = np.array(h5f['timestamps'], dtype = int).ravel()
    return config[0].decode('UTF-8'), timestamps, phases*1e-4
config, timestamps, phases = fringeezz_load_h5_data_log(file,filepath) #calls the function defined above to extract data from your hdf5 file 
print(config) #To get info on the configuration of the Fringeezz at the time of the measurement :

# so now, let's deal with these raw data and turn it into nice basic lists :
    
list_phases = phases[0]
for i in np.arange(1,np.size(phases[:,1])):
    list_phases = np.append(list_phases, phases[i])
ecart_type = int( np.sqrt(np.mean(list_phases**2)-(np.mean(list_phases))**2) * 1000)
#getting rid of the 'nan' :
for i in np.arange(0,np.size(list_phases)-1):
    if list_phases[i] != 0:
        print(i) 
        start = i
        break
list_phases = list_phases[start:np.size(list_phases)] - np.mean(list_phases[start:np.size(list_phases)])
#we don't actually use the time stamps as we know there is 1ms between each phase measurement
time = np.arange(0,np.size(phases[0])*np.size(phases[:,1]))
time = (time[start:np.size(time)]-start)/60000

       
       

##m-shots averages 
#m = 3
#list_phasesm = [] ; timem =[] ; ecart_typem = [];
#for k, value in enumerate (np.arange(0,np.size(list_phases)/m)):
#    list_phasesm = np.append ( list_phasesm, np.mean (list_phases[k*m:k*m+m]) )
#    timem =  np.append ( timem ,  np.mean(time[k*m:k*m+m]) )
#    print (str(k) + '/' + str (np.size(list_phases)/m))
#ecart_typem = int( np.sqrt(np.mean(list_phasesm**2)-(np.mean(list_phasesm))**2) * 1000)
#plt.figure(figsize=(7,3))
#plt.plot(timem, list_phasesm, '.', color = 'navy')   
#plt.ylim([-np.pi,np.pi])
#plt.yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ['$-\pi$', '-$\pi$/2', '0', '$\pi$/2' , '$\pi$']) 
#plt.xlabel('time (min)')
#plt.ylabel('CEP (rad)')
#plt.title(str(ecart_typem)+ ' mrad rms (averages over ' + str(m)+ ' shots)')
#plt.savefig(str(output_dir) + '/' +  str(file[0:len(file)-5]) + '_averaged'+str(m)+ 'shots' + '.png', dpi = 300, bbox_inches='tight')
       
       




##to save CEP values as a simple basic txt file :
#newfile=open(str(output_dir) + '/' + str(filename) + '_txt_data.txt','w')  
#newfile.write('CEP(rad)' +'  ' + 'time(min)' +"\n")
#for i in np.arange(0,np.size(list_phases)):
#    newfile.write(str(list_phases[i])+ '  ' +  str(time[i]) + "\n")
#newfile.close()       








      
# SECOND : PLOT AND SAVE PHASE/TIME + HISTOGRAM
s = 12  #size of the text
c = 'navy'  #color
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(9,3), gridspec_kw = {'width_ratios':[2.3, 0.7]})
plt.subplots_adjust(left=0.05, right=0.97, wspace=0.01)
ax1.plot(time, list_phases, '.', markersize = 1, color = c)
ax1.set_ylim([-np.pi,np.pi])
ax1.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
ax1.set_yticklabels(['$-\pi$', '-$\pi$/2', '0', '$\pi$/2' , '$\pi$'])
ax1.set_xlim([0,max(time)])
ax1.tick_params(labelsize=s)
ax1.set_xlabel('time (min)', fontsize=s)
ax1.set_ylabel('CEP (rad)', fontsize=s)
ax2.hist(list_phases, bins=25, orientation='horizontal',color = c)
ax2.set_xlabel('occurrence', fontsize=s)
ax2.tick_params(labelsize=s)
ax2.annotate('$\sigma_{rms} = $' + str(ecart_type) + ' mrad', xy=(100,2), size = s +2)
ax2.set_axis_off()
f.savefig(str(output_dir) + '/' + str(file[0:len(file)-5]) + '_Phase_and_Histo' + '.png', dpi=1000,bbox_inches='tight')




#THIRD : PLOT PSD_IPN
reprate = 1000.  #repetition rate of the laser = 1kHz
def intg_phase_noise(df, PSD):
    """ Calculate the inegrated phase noise from the PSD
    input: df: frequency step of the PSD in Hz
    input: PSD: the PSD in rad^2/Hz
    return: ipn : integrated phase noise in rad"""
    ipn = np.sqrt((df*PSD[::-1]).cumsum())[::-1]
    return ipn

def PSD_IPN(phases, reprate, precision = 1.0):
    """ Calculate PSD and IPN from fringeezz measurement (data_log)
    input : phases (direct from HDF5 data_log file 'phases'), reprate of the measure, and precision of the PSD calculation in Hz (default = 1Hz)
    Print stats about the measurement: mean and std
    return freq (frequencies axes in Hz), PSD (rad^2/Hz) and IPN (mrad)"""
    N = np.ceil(np.log2(reprate/precision)) # auto detect 2**N nperseg to reach 'precision'Hz with PSD
    phi = phases.ravel()
    print('mean: {:.3f} rad, std: {:.1f} mrad RMS'.format(np.mean(phi), np.std(phi)*1e3)) # some stats
    fw, psd = welch(phi, reprate, nperseg=2**N, scaling='density') # PSD
    ipn = intg_phase_noise(fw[1]-fw[0], psd) # integrated phase noise
    return fw, psd, ipn

freq, psd, ipn = PSD_IPN(phases,reprate, precision = 0.1) 
                        
def plot_save_PSD(freq, psd, ipn):
    """plot and save PSD and IPN
    input: freq, psd and ipn direct from PSD_IPN function
    return: nothing
    plot the graph
    save the graph """
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].set_title('Power spectral density')
    axarr[0].set_ylabel('PSD [rad$^2$/Hz]')
    axarr[1].set_title('Integrated phase noise')
    axarr[1].set_xlabel('Frequency [Hz]')
    axarr[1].set_ylabel('IPN [mrad]')

    axarr[0].loglog(freq, psd, label = 'psd', color = c)
    axarr[1].semilogx(freq, ipn*1e3, label = 'ipn', color = c)

    axarr[0].legend(loc='best')
    axarr[1].legend(loc='best')
    plt.savefig(str(output_dir) + '/' +  str(file[0:len(file)-5]) + '_PSD_IPN' + '.png', dpi = 300, bbox_inches='tight')

plot_save_PSD(freq, psd, ipn) 

print('2 plots have been saved in the same folder as the file : PSD and IPN + CEP vs time')
print('A text file with CEP values (centered around 0) has also been saved')