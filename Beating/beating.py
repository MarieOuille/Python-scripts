import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import math
from scipy import signal

t = np.linspace(-1,1,2*10000, endpoint=False)


E=1


def A(x,w):
    return E*np.cos(w*x)

plt.plot(t, A(t,100))
plt.plot(t,A(t,110))
plt.ylim(-2,2)
plt.xlabel('time (u.a.)')
plt.ylabel('amplitude (u.a.)')
plt.xlim(-1,0.3)
plt.savefig('close_frequencies.png', dpi=500)
plt.close()
def Somme(p,k):
    return p+k
plt.plot(t, (A(t,110)))
#plt.plot(t, (A(t,100)))
plt.xlim(-1,0.3)
plt.ylim(-2,2)
plt.xlabel('time (u.a.)')
plt.ylabel('amplitude (u.a.)')
plt.savefig('freq1.png', dpi=500)
plt.close()

plt.plot(t, Somme(A(t,100),A(t,110)))
plt.xlabel('time (u.a.)')
plt.ylabel('amplitude (u.a.)')
plt.savefig('freq1plus2.png', dpi=500)
plt.close()
