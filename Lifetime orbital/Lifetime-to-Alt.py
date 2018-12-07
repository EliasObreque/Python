# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 12:31:52 2018

@author: ELias Obreque.

REf: DAS (NASA) and STK simulators were used
"""

import matplotlib.pyplot as plt

masa3U  = 1.13*3
masa6U  = 1.13*6
masa12U = 1.13*12

Area3U  = 0.055
Area6U  = 0.055
Area12U = 0.08


Area2masa3U = Area3U/masa3U
Area2masa6U = Area6U/masa6U
Area2masa12U = Area12U/masa12U
Alt = [350, 360, 370, 380, 420, 430, 440, 450]


print('\nRelacion de Area-to-masa')
print('3U: ',Area2masa3U, '|6U: ', Area2masa6U,'|12U: ', Area2masa12U)

info3U_0  = '2020, Inclinacion: 98'
info3U_1  = '2024, Inclinacion: 98'


Olt3U_0 = [0.531, 0.635, 0.750, 0.860, 1.325, 1.446, 1.566, 1.692] 
Olt6U_0 = [0.838, 0.969, 1.101,1.232, 1.785, 1.938, 2.097,  2.267] 
Olt12U_0 = [1.008, 1.144, 1.287, 1.429, 2.048, 2.229, 2.415,2.623] 

Olt3U_1 = [0.093 ,0.110 ,0.131 ,0.153 , 0.307, 0.361, 0.422, 0.498]
Olt6U_1 = [0.175,0.214 ,0.252 , 0.301, 0.613,0.739 ,0.887 ,1.073]
Olt12U_1 = [0.241,0.290 , 0.345,0.416 , 0.871, 1.057, 1.287, 1.593]


plt.figure()
plt.title('Año lanzamiento: '+ info3U_0)
plt.plot(Alt, Olt3U_0)
plt.plot(Alt, Olt6U_0)
plt.plot(Alt, Olt12U_0)
plt.grid()
plt.xlabel('Altura [km]')
plt.ylabel('Lifetime [año]')
plt.legend(['3U', '6U', '12U'])

plt.figure()
plt.title('Año lanzamiento: '+ info3U_1)
plt.plot(Alt, Olt3U_1)
plt.plot(Alt, Olt6U_1)
plt.plot(Alt, Olt12U_1)
plt.grid()
plt.xlabel('Altura [km]')
plt.ylabel('Lifetime [año]')
plt.legend(['3U', '6U', '12U'])

plt.figure()
plt.title('Comparación de vida orbital')
plt.plot(Alt, Olt3U_0, 'b--')
plt.plot(Alt, Olt6U_0, 'r--')
plt.plot(Alt, Olt12U_0, 'g--')
plt.plot(Alt, Olt3U_1, 'b')
plt.plot(Alt, Olt6U_1, 'r')
plt.plot(Alt, Olt12U_1, 'g')
plt.grid()
plt.xlabel('Altura [km]')
plt.ylabel('Vida orbital [año]')
plt.legend(['3U_2020', '6U_2020', '12U_2020','3U_2024', '6U_2024', '12U_2024'])
