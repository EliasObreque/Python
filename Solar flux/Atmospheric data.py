# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 23:38:26 2018

@author: EO
"""


import matplotlib.pyplot as plt


SF_open    = open('Solarflux.txt','r')
SF_read    = SF_open.read()
SF_Data    = SF_read.split('\n')
SF_open.close()

def getSF(Data):
    leng = len(Data)
    sf    = [ ]
    jd    = [ ]
    for k in range(leng-2):
        
        if SF_Data[k].split(' ')[3] != '':
            
            sf.append(float(Data[k].split(' ')[3]))
            jd.append(float(Data[k].split(' ')[0]))      
        else:
            sf.append(float(Data[k].split(' ')[4]))
            jd.append(float(Data[k].split(' ')[0]))
           
    return sf, jd

SF, JD = getSF(SF_Data)


plt.plot(JD, SF)
plt.grid()
