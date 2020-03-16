import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(-30*10**(-15),30*10**(-15),2*10000, endpoint=False)
f=3/8*10**15   #800nm
a=(3.5/2.35)*10**(-15)    #durée en fs fwhm = 2.35 a

def A(x):
    return np.exp(-x**2/(2*a**2))

def E(x, cep):
    return A(x)*np.exp(-1j*(2*np.pi*f*x+cep))


CEPin = input('What CEP value do you want (in multiple of Pi, and no fraction please) ?')

plt.plot(t, A(t), 'b--')
plt.plot(t, -1*A(t), 'b--')
plt.plot(t, (E(t, float(CEPin) * np.pi)).real, color='red')
plt.xlabel('Time (fs)', fontsize='16')
plt.ylabel('Electric field (a.u.)', fontsize='16')
plt.yticks((-1,-0.5,0,0.5,1))
plt.xlim(-6*10**(-15),6*10**(-15))
plt.xticks((-5*10**(-15),-2.5*10**(-15), 0, 2.5*10**(-15), 5*10**(-15)), ('-5', '-2.5', '0', '2.5', '5'))
plt.axis('Off')


save = input('Do you want to save the file ? (yes or no) ')
if save == 'yes':
    plt.savefig('carrier_enveloppe.png')
else :
    pass
