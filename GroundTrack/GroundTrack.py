"""
Created on Wed Oct 31 11:23:34 2018

@author: Elias Obreque
"""

#Satellite ground track with TLE
#================================
#TLE must be. txt file
#SUCHAI
#1 42788U 17036Z   18203.83802595  .00001168  00000-0  54122-4 0  9994
#2 42788  97.3987 263.1397 0010661 317.6944  42.3469 15.22024574 60010
#================================
# librery: 
# - matplotlib.ticker, cartopy.crs, cartopy.mpl.gridliner, TLESatellite
#================================

#Import
import matplotlib.ticker as mticker
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import math
import time
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from TLESatellite import Sat


#Santiago, Chile
SanLat = -33.4541
SanLon = -70.6502
SanAlt =  570

#=============================================================================

NombreTle = 'suchai'
Area   = 0.01
Masa = 1.13
CD = 2.2

#Julian day initial
t1      = time.strftime('%X %x %Z')
       
 #=============================================================================
#     Function      
 
def JulianDate(t):
    m   = float(t[9] + t[10])
    if float(t[0] + t[1])>= 20:
        d = float(t[12] + t[13])+1
    else:
        d = float(t[12] + t[13])
        
    y   = float('20' + t[15] + t[16])
    UT  = (3 + float(t[0] + t[1]) + float(t[3] + t[4])/60.0 + float(t[6] + t[7])/3600.0) % 24.0
    J0  = 367.0*y - math.floor(7.0*(y + math.floor((m + 9.0)/12.0))/4.0) + math.floor(275.0*m/9.0) + d + 1721013.5
    JD  = J0 + UT/24.0
    return JD

def JulianEpoch(y):
    m   = 1
    J0  = 367.0*y - math.floor(7.0*(y + math.floor((m + 9.0)/12.0))/4.0) + math.floor(275.0*m/9.0) + 1721013.5
    return J0

def Kepler_ite(M, e):
    err     = 1.0e-10
    E0      = M
    t       = 1
    itt     = 0
    while(t):
          E     =  M + e*math.sin(E0);
          if ( abs(E - E0) < err):
              t = 0
          E0    = E
          itt   = itt + 1
          E     = E0
    return E

def Plane(M0, n, t, ec, a):
    E   = []
    r   = []
    f   = []
    x   = []
    y   = []
    vr  = []
    vf  = []
    p   = a*math.sqrt(1 - ec**2)
    h   = math.sqrt(p*mu)
    M   = (M0 + n*t)%(2*math.pi)
    #Iter
    k   = 0
    for ma in M:
        E.append(Kepler_ite(ma, ec))
        
    for ee in E:
        r.append(a*(1 - ec*math.cos(ee)))
        f.append(2*math.atan(math.sqrt((1.0 + ec)/(1.0 - ec))*math.tan(0.5*ee)))
        x.append(r[k]*math.cos(f[k]))
        y.append(r[k]*math.sin(f[k]))
        vr.append((mu/h)*math.sin(f[k]))
        vf.append((mu/h)*(ec + math.cos(f[k])))
        k = k + 1 
    return [r, f, x, y, M, vr, vf]

def VelPlane(a, e, f, mu):
    vx  = []
    vy  = []
    p   = a*math.sqrt(1 - e**2)
    h   = math.sqrt(p*mu)
    for ang in f:
        vx.append((mu/h)*math.sin(ang))
        vy.append((mu/h)*(e + math.cos(ang)))
    return [h, vx, vy]

def VelViento(x, y, z):
    omegav  = np.array([0, 0, rot_earth])
    pos     = np.array([x, y, z])
    velv    = np.cross(omegav, pos)*R_earth/np.linalg.norm(pos)
    V       = np.linalg.norm(velv)
    return velv[0], velv[1], velv[2], V

def VelRelativa(vx, vy, vz, Vx, Vy, Vz):
	vxr = Vx - vx
	vyr = Vy - vy
	vzr = Vz - vz
	vr 	= math.sqrt(vxr**2 + vyr**2 + vzr**2) 
	return vxr, vyr, vzr, vr

def ECI(t, x, y, W, W_dot, RAAN, RAAN_dot, i):
	X   	= [ ]
	Y   	= [ ]
	Z   	= [ ]
	k   	= 0
	RaAn 	= [ ]
	WW 		= [ ]
	for DT in t:
        #Perifocal frame to geocentric equatorial frame
		 w  = W + W_dot*DT
		 raan  = RAAN + RAAN_dot*DT
		 RaAn.append(math.degrees(raan))
		 WW.append(math.degrees(w))
		 PO = [[- math.sin(raan)*math.sin(w)*math.cos(i) + math.cos(raan)*math.cos(w), - math.sin(raan)*math.cos(w)*math.cos(i) - math.cos(raan)*math.sin(w), 0],
               [math.cos(raan)*math.sin(w)*math.cos(i) + math.sin(raan)*math.cos(w), math.cos(raan)*math.cos(w)*math.cos(i) - math.sin(raan)*math.sin(w), 0], 
               [math.sin(w)*math.sin(i), math.cos(w)*math.sin(i), 0]]
		 X.append(PO[0][0]*x[k] + PO[0][1]*y[k])
		 Y.append(PO[1][0]*x[k] + PO[1][1]*y[k])
		 Z.append(PO[2][0]*x[k] + PO[2][1]*y[k])
		 k = k + 1     
	return [X, Y, Z, RaAn, WW]

def ECIV(t, x, y, W, W_dot, RAAN, RAAN_dot, i):
	X   = [ ]
	Y   = [ ]
	Z   = [ ]
	k   = 0
	for DT in t:
		w  = W + W_dot*DT
		raan  = RAAN + RAAN_dot*DT
		PO = [[- math.sin(raan)*math.sin(w)*math.cos(i) + math.cos(raan)*math.cos(w), - math.sin(raan)*math.cos(w)*math.cos(i) - math.cos(raan)*math.sin(w), 0],
               [math.cos(raan)*math.sin(w)*math.cos(i) + math.sin(raan)*math.cos(w), math.cos(raan)*math.cos(w)*math.cos(i) - math.sin(raan)*math.sin(w), 0], 
               [math.sin(w)*math.sin(i), math.cos(w)*math.sin(i), 0]]
		Z.append((PO[2][0]*x[k] + PO[2][1]*y[k]))
		X.append(-(PO[0][0]*x[k] + PO[0][1]*y[k]))
		Y.append((PO[1][0]*x[k] + PO[1][1]*y[k]))
		k = k + 1
	return [X, Y, Z]

def RotationEarth(X, Y, Z, GM, rot_earth):
    Xe  = [ ]
    Ye  = [ ]
    Ze  = [ ]
    ke  = 0
    for DT in t:
        theta = GM*math.pi/180.0 + rot_earth*DT
        RO = [[math.cos(theta), math.sin(theta), 0],
               [-math.sin(theta), math.cos(theta), 0],
               [0, 0, 1]]
        Xe.append(RO[0][0]*X[ke] + RO[0][1]*Y[ke])
        Ye.append(RO[1][0]*X[ke] + RO[1][1]*Y[ke])
        Ze.append(Z[ke])
        ke = ke + 1
    return [Xe, Ye, Ze]

def GroundTrack(Xe, Ye, Ze, r):
    lat 	= [ ]
    lon 	= [ ]
    alt   = [ ]
    #Latitude and longitude (geocentric equatorial frame)
    for k in range(len(Xe)):
        alt.append(r[k] - R_earth)
        lon.append(math.degrees(math.atan2(Ye[k],Xe[k])))
        lat.append(math.degrees(math.asin(Ze[k]/r[k])))
    return [lon, lat, alt]

def PerturbationJ2(a, ec, mu, R_earth, i):
    Const    = -(1.5)*J2*(math.sqrt(mu))*((R_earth)**2)/(((1-ec**2)**2)*((a)**(3.5)))
    RAAN_dot = Const*(math.cos(i)) 
    W_dot    = Const*(2.5*(math.sin(i))**2 - 2)
    M_dot    = 0.5*Const*math.sqrt(1 - ec**2)*(3*(math.cos(i))**2 - 1 )   
    return RAAN_dot, W_dot, M_dot

def JD2GMST(JD):
    JD0     = np.round(JD) - 0.5
    UT      = (JD - JD0)*24.0      #%Time in hours past previous midnight
    D0      = JD0 - 2451545.0   #%Compute the number of days since J2000
    T       = D0/36525          #%Compute the number of centuries since J2000
    GMST0   = 100.4606184 + 36000.77004*T + 0.000387933*T**2.0 - (2.583e-8)*T**3.0
    GMST0   = GMST0 - math.floor(GMST0/360.0)*360.0
    GMST    = GMST0 + 360.98564724*UT/24.0
    return GMST

def animate(i):
    GT.set_data(lonsat[i], latsat[i])  # update the data
    return GT,

def init():
    GT.set_data([],[])
    return GT, 

#=============================================================================
    
InfoTLE = Sat.LeerTLE(NombreTle + '.txt')
sat = Sat.setTLE(InfoTLE)

#LINE 1 TLE

#Time
JD1     	= JulianDate(t1)
JD2     	= JD1 + 12/24.0
GM      	= JD2GMST(JD1)%360

Yepoch   	= sat.EY
J0epoch  	= JulianEpoch(Yepoch)
Jepoch   	= J0epoch + sat.DEY

dt       	= 5 #s
tset     	= (JD1 - Jepoch)*24*3600
t        	= np.linspace(0.0, (JD2 - JD1)*24.0*3600.0, int((JD2 - JD1)*24.0*3600.0/dt))

#Earth grav. parameter
#=============================================================================
rot_earth   = 7.292115854670501e-5
R_earth     = 6378.1350
m_earth     = 5.9722e24
G           = 6.674e-11
mu          = 3.986044418e5
J2          = 1.08263e-3
k_2         = (5.413080e-4)*R_earth**2

#ORBITAL ELEMENTS
isat        = sat.incl*math.pi/180.0
ecsat       = sat.ec/10000000.0
nsat        = (2.0*math.pi*sat.nm/86400.0)
asat        = (mu**(1.0/3.0))/(nsat**(2.0/3.0))
[RAAN_dotsat, W_dotsat, M_dotsat] = PerturbationJ2(asat, ecsat, mu, R_earth, isat)
RAANsat     = sat.RAAN*math.pi/180.0 + RAAN_dotsat*tset   
Wsat        = sat.ap*math.pi/180.0 + W_dotsat*tset
M0sat       = (sat.MA*math.pi/180.0 + nsat*tset + M_dotsat)%(2*math.pi)

#=============================================================================

[rsat, fsat, xsat, ysat, Msat, VR, VF]  = Plane(M0sat, nsat, t, ecsat, asat)
[Hsat, v_x, v_y]          				    = VelPlane(asat, ecsat, fsat, mu)
[Xsat, Ysat, Zsat, RAan, W]          	 = ECI(t, xsat, ysat, Wsat, W_dotsat, RAANsat, RAAN_dotsat, isat)
[Vxsat, Vysat, Vzsat] 					    = ECIV(t, v_x, v_y, Wsat, W_dotsat, RAANsat, RAAN_dotsat, isat)
[Xesat, Yesat, Zesat]           		    = RotationEarth(Xsat, Ysat, Zsat, GM, rot_earth)
[Vxesat, Vyesat, Vzesat]                = RotationEarth(Vxsat, Vysat, Vzsat, GM, rot_earth)
[lonsat, latsat, Alt]             		 = GroundTrack(Xesat, Yesat, Zesat, rsat)


Vxv 	= []
Vyv 	= []
Vzv 	= []
Vv    = []
Vxr	= []
Vyr 	= []   
Vzr 	= []
Vr 	= []	
Vsat 	= []

for k in range(len(lonsat)):
	[vxv, vyv, vzv, vv] = VelViento(Xesat[k], Yesat[k], Zesat[k])
	Vxv.append(vxv)
	Vyv.append(vyv)
	Vzv.append(vzv)
	Vv.append(vv)
	[vxr, vyr, vzr, vr] = VelRelativa(vxv, vyv, vzv, Vxesat[k], Vyesat[k], Vzesat[k])
	Vxr.append(vxr)
	Vyr.append(vyr)
	Vzr.append(vzr)	
	Vr.append(vr)
	Vsat.append(math.sqrt(Vxesat[k]**2 + Vyesat[k]**2 + Vzesat[k]**2))
	  
#=============================================================================

#PLOT GroundTrack
print(lonsat[0])
print(latsat[0])
plt.figure(1)
fig = plt.gcf()
fig.set_size_inches(9, 5)
ax 	= plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
gl 	= ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_left = False
gl.xlines = True
gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER


plt.ylim([-90, 90])
plt.xlim([-180, 180])
plt.plot(SanLon, SanLat,'+k')
plt.plot(lonsat,latsat,'.', markersize = 1)


GT,     = ax.plot([], [], '*r')



ani = animation.FuncAnimation(fig, animate, np.arange(1, 2500), init_func=init,
                              interval = dt*1000, blit=True)
plt.show()


# Plot Altitude, latitude and longitude
plt.figure()
plt.plot(t, Alt)
plt.xticks()
plt.plot(t, lonsat)
plt.plot(t, latsat)
plt.grid()	


