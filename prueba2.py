# -*- coding: utf-8 -*-
"""
Created on Mon May  4 03:54:48 2020

@author: maria
"""

import numpy as np
import math 
from math import sin #exp #
from matplotlib import pyplot as plt

# Constantes fisicas
epsilon0 = 8.8542*10**-12   # Valor de la constante epsilon_o
miu0 = 1.2566*10**-6   # Valor de la constante mu_o
epsilonr1 = 12  # Constante para el plastico
c = 1/ math.sqrt(miu0*epsilon0)  ## Valor velocidad de la luz en el vacio m/s^2

# Declaracion de variables
Pt = 250     # Definición de mi mundo
El = 0.0254/2  # 1 pulgada = 0.0254 m  1/2''
Pc = int(Pt/2) # Mitad de mi mundo
Hy = np.zeros(Pt)
Ex = np.zeros(Pt)
Exr_t = np.zeros(Pt)
Ext_t = np.zeros(Pt)
cons = np.zeros((1,Pt) , dtype = float) 
dx = 0.01   # Discretización de mi mundo
dt = dx/c    # Calcular el dt segun el dx que tengamos y c
E0 = 1.0
Total = 1000
pi = 3.1416
f = 10.5*10**9
d2 = ((c/f)/4*1.49)  # Espesor de las capas antireflejo con un indice de refraccion entre el del aire y el del vidrio
epsilonr0 = 50 

# Coondiciones de frontera
Fron_U = [0, 0]
Fron_D = [0, 0]
cb = np.ones(Pt)
cb = 0.5 * cb
cb_start = 100
cb[cb_start:] = 0.5 / epsilonr1

med1i = int((Pt*d2)/0.0254)  # Espacio de las capas antireflejo
med2i = int((Pt*((El)-1/4))/0.0254)  # Espacio del aire
med3i = int(med2i + (((El)*Pt))/0.0254)  # Espacio del plastico
med4i = med1i

cons[:,1]= miu0 
cons[0:med1i,0] = epsilon0
cons[med1i:med2i,0] = epsilonr0*epsilon0
cons[med2i:med3i,0] = epsilonr1*epsilon0
cons[med3i:med4i,0] = epsilonr0*epsilon0
cons[med4i:Pt,0] = epsilon0

# Coeficioentes de transmisión y reflexión
n0 = math.sqrt( miu0/epsilon0)  # Impedancia intrinseca del aire
n1 = math.sqrt( miu0/(epsilon0*epsilonr1))  # Impedancia intrinseca del plastico
#n2 = 
R1 = (n1 - n0)/(n1 + n0)  # Coeficiente de reflexión cuando va del aire al plastico
T1 = (2*n1)/(n1 + n0)    # Coeficiente de transmición cuando va del aire al plastico
R2 = (n0 - n1)/(n0 + n1)  # Coeficiente de reflexión cuando va del plastico al aire
T2 = (2*n0)/(n0 + n1)   # Coeficiente de transmición cuando va del plastico al aire

# Method FDTD in 2D
for n in range(1, Total+1):
    
    #Se calcula Ex       
    for i in range(1, Pt):
        Ex[i] = Ex[i] + cb[i] * (Hy[i - 1] - Hy[i])
        Exr_t[i] = R1*Ex[i] + cb[i] * (Hy[i - 1] - Hy[i])
        Ext_t[i] = T1*Ex[i] + cb[i] * (Hy[i - 1] - Hy[i])
        
    # Absorvente
    Ex[0] = Fron_D.pop(0)
    Fron_D.append(Ex[1])
    Ex[Pt - 1] = Fron_U.pop(0)
    Fron_U.append(Ex[Pt - 2])
     #Ez[n+1,imax] = Ez[n,imax-1]+((c*(dt-dx))/(c*(dt+dx)))*(Ez[n+1,imax-1]-Ez[n,imax])
    
    #Se calcula Hy 
    for i in range(Pt - 1):
        Hy[i] = Hy[i] + 0.5 * (Ex[i] - Ex[i + 1])
    
    # Se pone un pulso    
     #Ez[Pc,Pc]= E0*exp((-(n-8.0)**2)/16.0)
    pulse = sin(2 * pi * f * dt * Total) 
    Ex[5] = pulse + Ex[5]
#    pulse = exp(-0.5 * ((t0 - time_step) / spread) ** 2)
#    Ez[Pci, Pcj] = pulse

# Grafica
X0 = (np.linspace(-50,50,num = Pt))*dx
plt.plot(X0,Hy[0,:]) #plot pulse
fig1 = plt.figure()
f1 = plt.contourf(Hy, 15, cmap = 'Set1')
plt.colorbar(f1)
plt.ylabel('Y(m)')
plt.xlabel('X(m)')
plt.show()
        
