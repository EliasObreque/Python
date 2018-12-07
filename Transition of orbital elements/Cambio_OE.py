# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 17:25:33 2018

@author: Elias Obreque
"""

import matplotlib.pyplot as plt
import math
import numpy as np
from sympy import Integral, integrate


NombreTle = 'sat42788'



m_earth     = 5.9722e24
G           = 6.674e-11
mu          = 3.986044418e5
R_earth     = 6378.0000

def JulianEpoch(y):
    m   = 1
    J0  = 367.0*y - math.floor(7.0*(y + math.floor((m + 9.0)/12.0))/4.0) + math.floor(275.0*m/9.0) + 1721013.5
    return J0

def getOE(Data):
    leng = len(Data)
    i    = [ ]
    e    = [ ]
    n    = [ ]
    a    = [ ]
    t    = [ ]
    for k in range(leng-1):
        if Data[k][0] == '2':
            i.append(float(Data[k][9:16]))
            e.append(float(Data[k][26:33])/10000000.0)
            n.append(float(Data[k][52:63]))
            nsat = 2.0*math.pi*float(Data[k][52:63])/86400.0
            a.append(((mu**(1.0/3.0))/(nsat**(2.0/3.0))) - R_earth)
            
        elif Data[k][0] == '1':
            if float(Data[k][18:20])<20:
                year = float('20'+Data[k][18:20])
            else:
                year = float('19'+Data[k][18:20])    
            J0epoch = JulianEpoch(year)
            epoch   = float(Data[k][20:32])
            Jepoch  = J0epoch + epoch
            t.append(Jepoch)    
    return i, a, n, e, t, 


def Trap(fx, x):
    area = 0
    for i in range(len(t) - 1):
        area = area + 0.5*(x[i + 1] - x[i])*(fx[i + 1] + fx[i])
    return area

TLE_open    = open(NombreTle+'.txt','r')
TLE_read    = TLE_open.read()
Data        = TLE_read.split('\n')
TLE_open.close()

inc, a, n, e, tJD = getOE(Data)

t = (np.array(tJD) - tJD[0])

NumRev = Trap(n,t)

KilTrav = 2*math.pi*(np.mean(a) + R_earth)*NumRev

print("Number of revolutions:", NumRev)
print("Kilometers traveled:", KilTrav)

plt.figure()
plt.title('Inclination')
plt.ylabel('$i$ [°]')
plt.xlabel('Time [day]')
plt.plot(t, inc)
plt.grid()
plt.show()

plt.figure()
plt.title('Apogee')
plt.ylabel('$Apogee$ [km]')
plt.xlabel('Time [day]')
plt.plot(t, a)
plt.grid()
plt.show()

plt.figure()
plt.title('Revolutions by day')
plt.ylabel('$n$ [rev/día]')
plt.xlabel('Time [day]')
plt.plot(t, n)
plt.grid()
plt.show()

plt.figure()
plt.title('Eccentricity')
plt.ylabel('$e$ [-]')
plt.xlabel('Time [day]')
plt.plot(t, e)
plt.grid()
plt.show()





