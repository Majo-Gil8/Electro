# -*- coding: utf-8 -*-
"""
Created on Mon May  4 03:54:48 2020

@author: maria
"""

import numpy as np
import math 
from math import   exp #sin
from matplotlib import pyplot as plt

# Constantes fisicas
epsilon0 = 8.8542*10**-12   # Valor de la constante epsilon_o
miu0 = 1.2566*10**-6   # Valor de la constante mu_o
epsilonr1 = 12  # Constante para el plastico
c = 1/ math.sqrt(miu0*epsilon0)  ## Valor velocidad de la luz en el vacio m/s^2

# Declaracion de variables
Pt = 100     # Definición de mi mundo
El = 0.0254/2  # 1 pulgada = 0.0254 m  1/2''
Pc = int(Pt/2) # Mitad de mi mundo
Total = 100
Hy = np.zeros((Pt,Pt), dtype = "float64") 
Ez = np.zeros((Pt,Pt), dtype = "float64") 
cons1 = np.zeros(Pt)
cons2 = np.zeros(Pt) 
dx = 0.001   # Discretización de mi mundo
dt = dx/c    # Calcular el dt segun el dx que tengamos y c
E0 = 500
pi = 3.1416
t0 = 40
spread = 12
f = 10.5*10**9   # Frecuencia utilizada 
n0 = math.sqrt( miu0/epsilon0)  # Impedancia intrinseca del aire
n1 = math.sqrt( miu0/(epsilon0*epsilonr1))  # Impedancia intrinseca del plastico
factor_norm = math.sqrt(epsilon0/miu0)

# Caracteristicas del material
Ir = 1.25  # Indice de refracción
lamda = c/f   # Longitud de onda 
m = 2   # Cualquier entero
d2 = ((m-1/2)*lamda)/2*Ir  # Espesor de las capas anti-reflejo con un indice de refraccion entre el del aire y el plastico
epsilonr2 = (miu0*Ir**2)/(epsilon0*n0)   # Constante para las capas
n2 = math.sqrt(miu0/(epsilon0*epsilonr2))  # Impedancia intrinseca de las capas anti-reflejo

# Coondiciones de frontera
#Fron_U = [0, 0]
#Fron_D = [0, 0]
#cb = np.ones(Pt)
#cb = 0.5 * cb
#cb_start = 100
#cb[cb_start:] = 0.5 / epsilonr1

med1i = int(Pt*d2)  # Espacio capa de anti-reflejo
med2i = int(Pt-(med1i+El))  # Espacio del aire
med3i = int(Pt-(med1i+med2i))  # Espacio del plastico

cons1[:]= miu0 # Miu en todo el mundo
cons2[:med2i] = epsilon0  # Epsilon en el aire de mi mundo
cons2[med1i:med2i] = epsilonr2*epsilon0  # Epsilon en las capas reflectantes de mi mundo
cons2[med2i:] = epsilonr1*epsilon0  # Epsilon en el plastico de mi mundo

# Coeficioentes de transmisión y reflexión

R1 = (n1 - n0)/(n1 + n0)  # Coeficiente de reflexión cuando va del aire al plastico
T1 = (2*n1)/(n1 + n0)    # Coeficiente de transmición cuando va del aire al plastico
R2 = (n0 - n1)/(n0 + n1)  # Coeficiente de reflexión cuando va del plastico al aire
T2 = (2*n0)/(n0 + n1)   # Coeficiente de transmición cuando va del plastico al aire

# Method FDTD in 2D
for n in range(1, Total + 1):
        
    # Absorvente
     Fron_D = [0, 0]
     Fron_U = [0, 0]
     Ez[0] = Fron_D.pop(0)
     Fron_D.append(Ez[1])
     Ez[Pt - 1] = Fron_U.pop(0)
     Fron_U.append(Ez[Pt - 2])
#     Ez[Pt] = Ez[Pt-1]+((c*(dt-dx))/(c*(dt+dx)))*(Ez[Pt-1]-Ez[Pt])
     
     # Se pone un pulso    
     #pulse = sin(2 * pi * f * dt * Total)  
     #Ez[Pc]= E0*exp((-(n-8.0)**2)/16.0)
     pulse = exp(-0.5*((t0-n)/spread)**2)
     Ez[Pc] = pulse
    
     plt.clf()
    #Se calcula Hy 
     for i in range(0,(Pt - 1)):
        Hy[i] = Hy[i] + (dt/(cons1[i])*dx) * (Ez[i + 1] - Ez[i])
       
     #Se calcula Ex       
     for i in range(1, Pt):
        Ez[i] = Ez[i] + (dt/(cons2[i])*dx) * (Hy[i] - Hy[i - 1])        

# Grafica
     plt.subplot(211)
     plt.plot(Ez, color='r', linewidth=1)
     plt.ylabel('Campo Electrico [V/m]', fontsize='8')
     plt.xlim(0, Pt)
     plt.ylim(-(E0+0.2)*factor_norm, (E0+0.2)*factor_norm)
    
     plt.pause(0.01)
     plt.show()       
