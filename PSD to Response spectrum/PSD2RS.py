# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:13:08 2018

General approach

PSD1: MIL-STD-1540C Acceptance Level
PSD2: NASA random vibration for mass less than 22.7 kg (Acceptable condition in CubeSat)
Conservator Damping ratio = 0.05

@author: Elias Obreque
"""
import matplotlib.pyplot as plt
import numpy as np


# PSD input: MIL-STD-1540C

#HZ      = [20, 150, 600, 2000]
#GRMS2   = [0.0053, 0.04, 0.04, 0.0036]


# PSD input: NASA
HZ      = [20, 50, 800, 2000]
GRMS2   = [0.013, 0.08, 0.08, 0.013]

#==========================================================================

# Natural frequency 

FN = np.arange(20, 2010, 10)
xi = 0.05

#============================================================================
# Slope of the curve
m = []
for i in range(1, len(HZ)):
    base = np.log10(HZ[i]) - np.log10(HZ[i - 1])
    alt = np.log10(GRMS2[i]) - np.log10(GRMS2[i - 1])
    m.append(alt/base)

# GRMS2(f) 
def G_f(f):
    Gn = 0
    for i in range(1, len(HZ)):
        if f >= HZ[i - 1] and f <= HZ[i]:
            Const = GRMS2[i - 1]/HZ[i - 1]**m[i-1]
            Gn = Const*f**m[i-1]  
        elif f > max(HZ):
            Const = GRMS2[-1]/HZ[-1]**m[-1]
            Gn = Const*f**m[-1] 
    return Gn

def AreaRA(fx, x):
    areaRA = 0
    dx = x[1] - x[0]
    for i in range(len(x)):
        areaRA = areaRA + fx[i]*dx
    return areaRA**0.5


#============================================================================
areaPSD = 0 
for i in range(1, len(HZ)):
    base = HZ[i] - HZ[i -1]
    alt = np.abs(GRMS2[i] - GRMS2[i - 1])
    slope = (np.log(GRMS2[i]) - np.log(GRMS2[i - 1]))/(np.log(HZ[i]) - np.log(HZ[i - 1]))
    offset = GRMS2[i -1]/HZ[i - 1]**slope
    
    if slope != -1:
        areaPSD = areaPSD + (offset/(slope + 1))*(HZ[i]**(slope+1) - HZ[i - 1]**(slope+1))
    else:
        areaPSD = areaPSD + offset*(np.log(HZ[i]) - np.log(HZ[i - 1]))
        
Grms = np.sqrt(areaPSD)
GPeak = np.sqrt(2)*Grms
print('\nValor Grms: ', Grms)
print('Valor Gpeak: ', GPeak,"\n")

#==========================================================================

Acc = []

F = np.linspace(min(HZ), 2000, 10000)
df = F[1] - F[0]

k = 0
for fn in FN:
    Acc.append([])
    for i in range(len(F)):
        p = F[i]/fn
        C = (1 + (2*xi*p)**2)/((1 - p**2)**2 + (2*xi*p)**2)
        Acc[k].append(C*G_f(F[i]))
    k = k + 1


AreaGRMS = []
for m in range(len(FN)):
    AreaGRMS.append(AreaRA(Acc[m], F))
    #print("Response Accel (GRMS) [",FN[m], "Hz] =",AreaGRMS[m])


#==========================================================================

maxAccGRMS = max(AreaGRMS)
k = list(AreaGRMS).index(maxAccGRMS)
maxFn = FN[k]

print('Worst-Case point is:', maxAccGRMS,"g at", maxFn,"Hz")
    


#==========================================================================



# PLOT
#%%

textlegeng = []
plt.figure()
plt.title('Responce Power Spectral Density Curves')
plt.ylabel('PSD [$G^2 / Hz$]')
plt.xlabel('Frequency [$Hz$]')
plt.yscale('log')
plt.xscale('log')
for i in np.arange(0, int(100/20), 3):
    plt.plot(F, Acc[i], '--')
    textlegeng.append(str(FN[i]) +" Hz")
for i in np.arange(int(100/20) + 3, len(Acc), 50):
    plt.plot(F, Acc[i], '--')
    textlegeng.append(str(FN[i]) +" Hz")

plt.plot(HZ, GRMS2, 'k')
textlegeng.append("PSD")
plt.legend(textlegeng)
plt.ylim(0.001, 15)
plt.xlim(10, 10000)
plt.grid(which='both', axis='both')


plt.figure()
plt.title('Vibration Response Spectrum')
plt.ylabel('Accel [$G_{RMS}$]')
plt.xlabel('Frequency [$Hz$]')
plt.yscale('log')
plt.xscale('log')
plt.plot(FN, AreaGRMS)
plt.grid(which='both', axis='both')
