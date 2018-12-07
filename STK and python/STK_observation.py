# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 17:11:39 2018

Análsis de observación y pasados según ángulo de cámara [ConeHalfAng]
y número de satélites [Nsat].

@author: Elias Obreque


"""
import math
import numpy as np
import matplotlib.pyplot as plt
from win32api import GetSystemMetrics
from comtypes.client import CreateObject, GetActiveObject

def JulianDate(t):
        Date = t[0:11].split(' ')
        Month_num = {'Jan': 1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7,
                     'Aug':8, 'Oct':9, 'Sep':10,'Nov':11, 'Dec':12}
        
        Time = t[12:20].split(':')
        
        m   = Month_num[Date[1]]
        d   = float(Date[0])
        y   = float(Date[2])
        
        h    = float(Time[0])
        minu = float(Time[1])
        seg  = float(Time[2]) 
        
        UT  = h + minu/60 + seg/3600
        
        J0  = 367.0*y - math.floor(7.0*(y + math.floor((m + 9.0)/12.0))/4.0) + math.floor(275.0*m/9.0) + d + 1721013.5
        JD  = J0 + UT/24.0
        return JD   

#Crear objeto
uiApplication = CreateObject('STK11.Application')
uiApplication.Visible = True
uiApplication.UserControl=True

root = uiApplication.Personality2

root.NewScenario("Constellation_Analysis_5") 

   
from comtypes.gen import STKUtil
from comtypes.gen import STKObjects

scenario = root.CurrentScenario
#Create Static STK Object
scenario2 = scenario.QueryInterface(STKObjects.IAgScenario)
   

scenario2.SetTimePeriod("01 Jul 2018 00:00:00", "31 Dec 2018 00:00:00")

#1. Add a target object to the scenario.
if scenario.Children.Contains(STKObjects.eTarget,"Bustamante"):
    print('Target creado')
else:
    target = scenario.Children.New(STKObjects.eTarget,"Bustamante");
    target2 = target.QueryInterface(STKObjects.IAgTarget)

#2. Move the Target object to a desired location.

#Chillan
target2.Position.AssignGeodetic(-36.5923,-71.7643,0)



#%%
#Add and set a Satellite object to the scenario


Nsat    = 20 #numero de satelites

# Orbital elements
mu      = 3.986044418e14
R0      = 6371000
a       = 500000 #m
n       = 2*math.pi*math.sqrt((a + R0)**3/mu)  #MeanMotion
Rev     = 24*60*60/n
e       = 0   #eccentricity
RAAN    = 46
f       = 0   #True anomaly
incl    = 97  #Inclination
w       = 0   #ArgOfPerige

ConeHalfAng = 10


N = 1
vm = 0
ValorMedio = []
delay = []
Start       = []
Stop        = []

while N <= Nsat:
    delay.append([])
    f           = []
    satellite   = []
    satellite2  = []
    satProp     = []
    keplerian   = []
    keplerian2  = []
    
    sensor      = []
    sensor2     = []    
    
    # Para analisis equidistante entre satelite usar:    
    #phase   = 360/N
    
    # Para analisis de distancia fija
    phase   = 10
    
    for i in range(N):
        f.append(i*phase)
    
    #-------------------------------------------------------
        
    for i in range(N):
        if (scenario.Children.Contains(STKObjects.eSatellite, "SUCHAI"+str(i))):
            satellite.append(root.GetObjectFromPath('Satellite/SUCHAI'+str(i)))
        else:
            satellite.append(scenario.Children.New(STKObjects.eSatellite, "SUCHAI"+str(i)))
            
        satellite2.append(satellite[i].QueryInterface(STKObjects.IAgSatellite))
        satellite2[i].SetPropagatorType(STKObjects.ePropagatorTwoBody)
        satProp.append(satellite2[i].Propagator.QueryInterface(STKObjects.IAgVePropagatorTwoBody))
    
        # Use the Classical Element interface
        keplerian.append(satProp[i].InitialState.Representation.ConvertTo(STKUtil.eOrbitStateClassical))
        # Changes from SMA/Ecc to MeanMotion/Ecc
        keplerian2.append(keplerian[i].QueryInterface(STKObjects.IAgOrbitStateClassical))
        kp2 = keplerian2[i]
        kp2.SizeShapeType = (STKObjects.eSizeShapeMeanMotion)
        # Makes sure Mean Anomaly is being used
        kp2.LocationType = (STKObjects.eLocationMeanAnomaly)
        # Use RAAN instead of LAN for data entry
        kp2.Orientation.AscNodeType = (STKObjects.eAscNodeRAAN)
        
        # Set unit preferences for revs/day
        root.UnitPreferences.Item('AngleUnit').SetCurrentUnit('revs')
        root.UnitPreferences.Item('TimeUnit').SetCurrentUnit('day')
    
        # Assign the perigee and apogee altitude values:
        keplerian2[i].SizeShape.QueryInterface(STKObjects.IAgClassicalSizeShapeMeanMotion).MeanMotion = Rev
        keplerian2[i].SizeShape.QueryInterface(STKObjects.IAgClassicalSizeShapeMeanMotion).Eccentricity = e
    
    
        # Return unit preferences for degrees and seconds
        root.UnitPreferences.Item('AngleUnit').SetCurrentUnit('deg')
        root.UnitPreferences.Item('TimeUnit').SetCurrentUnit('sec')
    
        # Assign the other desired orbital parameters:
        kp2.Orientation.Inclination = incl
        kp2.Orientation.ArgOfPerigee = w
        kp2.Orientation.AscNode.QueryInterface(STKObjects.IAgOrientationAscNodeRAAN).Value = RAAN
        kp2.Location.QueryInterface(STKObjects.IAgClassicalLocationMeanAnomaly).Value = f[i]
    
        # Apply the changes made to the satellite's state and propagate:
        satProp[i].InitialState.Representation.Assign(keplerian[i])
        satProp[i].Propagate()
    
    
    for s in range(N):
        #Add and set a sensor
        if(satellite[s].Children.Contains(STKObjects.eSensor, "Camara")):
            sensor.append(root.GetObjectFromPath('Satellite/SUCHAI'+str(s)+'/Sensor/Camara'))
        else:
            sensor.append(satellite[s].Children.New(STKObjects.eSensor, "Camara"))
            
        sensor2 = sensor[s].QueryInterface(STKObjects.IAgSensor)
        sensor2.CommonTasks.SetPatternSimpleConic(ConeHalfAng, 0.1)
      
    access = []
    accessDP = []
    Results = []
    r = 0
    
    try:
        for acc in range(N):
            Results.append([])
            access.append(sensor[acc].GetAccessToObject(target))
            access[acc].ComputeAccess()
        
            accessDP.append(access[acc].DataProviders.Item('Access Data'))
            accessDP2 = accessDP[acc].QueryInterface(STKObjects.IAgDataPrvInterval)
        
            results = accessDP2.Exec(scenario2.StartTime, scenario2.StopTime)
        
            
            accessStartTimes = results.DataSets.GetDataSetByName('Start Time').GetValues()
            accessStopTimes = results.DataSets.GetDataSetByName('Stop Time').GetValues()
            
            Results[r].append(accessStartTimes)
            Results[r].append(accessStopTimes)
            
            r = r + 1
            
            
                    
    except:
        print('No existe acceso, aumentar ángulo de cono en la cámara')
        
    Startime    = []
    Stoptime    = []
    StartimeJD  = []
    StoptimeJD  = []
    estadoSat   = []
    TimeJD      = []

    k           = 0
    
    for i in range(N):
        StartimeJD.append([])
        StoptimeJD.append([])
        estadoSat.append([])
        TimeJD.append([])
    
        Startime.append(Results[i][0])
        Stoptime.append(Results[i][1])
    
        for pas in range(len(Results[i][0])):
            StartimeJD[k].append(JulianDate(Results[i][0][pas]))
            TimeJD[k].append(JulianDate(Results[i][0][pas]))
            Start.append(JulianDate(Results[i][0][pas]))
            estadoSat[k].append(1)
            StoptimeJD[k].append(JulianDate(Results[i][1][pas]))
            TimeJD[k].append(JulianDate(Results[i][1][pas]))
            Stop.append(JulianDate(Results[i][0][pas]))
            estadoSat[k].append(0)
        k = k + 1
    
    Pasada  = []
    Start.sort()
    Stop.sort()
    
    for i in range(1, int(len(Start)/N)):
        Pasada.append(i + 2)
        delay[N - 1].append(Start[i] - Stop[i - 1])
    
    ValorMedio.append(np.mean(delay[N - 1]))
    print('Valor Medio con', N,'Sat: ',ValorMedio[vm], 'dias')     
    N = N + 1
    vm = vm + 1 
    
    
    
#%%



CantSat = np.linspace(1, N-1, N-1)

plt.figure()
plt.grid()
plt.xlim([0, N-1])
plt.plot(CantSat, ValorMedio)
plt.xlabel('Número de satélites')
plt.ylabel('Intervalo respecto al ultimo ingreso [dias]')

    

for i in range(N-1):
    plt.figure()
    plt.grid()
    plt.title('Numero de satelites: ' + str(i + 1))
    plt.plot(np.array(delay[i]))
    plt.xlabel('Número de pasada')
    plt.ylabel('Intervalo respecto al ultimo ingreso [dias]')
plt.show()
