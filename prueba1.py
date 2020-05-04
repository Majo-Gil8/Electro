# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 21:55:51 2020

@author: maria
"""
#Main libraries    
import numpy as np
import math 
from math import exp #sin
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Constantes fisicas
epsilon0 = 8.8542*10**-12   # Valor de la constante epsilon_o
miu0 = 1.2566*10**-6   # Valor de la constante mu_o
epsilonr1 = 12  # Constante para el plastico
c = 1/ math.sqrt(miu0*epsilon0)  ## Valor velocidad de la luz en el vacio m/s^2

# Declaracion de variables
Pti = 50     # Definición de mi mundo
Ptj = 50     # Definición de mi mundo
El = 0.0254  # Grosor del material en m = 1/2''
M1 = [Pti,Ptj]  #Matriz para definir mi mundo
M2 = [El,El]   # Matriz para definir el plastico
Pci = int(Pti/2)
Pcj = int(Ptj/2)
Hx = np.zeros((Pti,Ptj) , dtype = "float64")
Hy = np.zeros((Pti,Ptj) , dtype = "float64")
Hxt = np.zeros((Pti,Ptj) , dtype = "float64")
Hyt = np.zeros((Pti,Ptj) , dtype = "float64")
Ez = np.zeros((Pti,Ptj) , dtype = "float64")
Ezt = np.zeros((Pti,Ptj) , dtype = "float64")
cons = np.zeros((Pti,Ptj) , dtype = "float64") 
dx = 0.001   # Discretización de mi mundo
dy = dx
dt = dx/c    # Calcular el dt segun el dx que tengamos y c
E0 = 1.0
Total = 100
pi = 3.1416
imax = 0.0254
jmax = imax
# Coondiciones de frontera
medio1i = int(Pti*(El-(El/2)))
medio1j = int(Ptj*(El-(El/2)))
medio2i = int(medio1i + (El*Pti))
medio2j = int(medio1i + (El*Pti))

cons[:,0]= miu0 
cons[0 : medio1i , 0] = epsilon0
cons[medio1i : medio2i , 2] = epsilonr1*epsilon0
cons[medio2i : Pti,1] = epsilon0
#print(cons)
# Absorvente
    

# Method FDTD in 2D
for n in range(1, Total + 1):
 
        # Se pone un pulso en el medio    
    Ez[Pci,Pcj]= E0*exp((-(n-8.0)**2)/16.0)
    #pulse = sin(2 * pi * 1500 * 1e6 * dt * time_step) 
#    pulse = exp(-0.5 * ((t0 - time_step) / spread) ** 2)
#    Ez[Pci, Pcj] = pulse
    
    # Se calculan Hy y Hx
    for i in range(0,(Pti - 1)):
        for j in range(0,(Ptj - 1)):
            Hy[i,j] = Hy[i,j] + (dt/(miu0)*dy)*((Ez[i + 1, j] - Ez[i , j]))
            Hx[i,j] = Hx[i,j] + (dt/(miu0)*dx)*((Ez[i,j] - Ez[i,j + 1]))
    # Se calcula Ez       
    for i in range(1, Pti):
        for j in range(1,Ptj):
            Ez[i,j] = Ez[i,j] + (dt/(epsilonr1)*dx)*(((Hy[i,j] - Hy[i - 1,j])/dy) - (Hx[i,j] + Hx[i,j - 1]))

# Matrices que guardan los datos para graf
Hyt = Hy
Hxt = Hx
Ezt = Ez
X0 = (np.linspace(-50,50,num = Pti))*dx
Y0 = (np.linspace(-50,50,num = Ptj))*dy


# Coeficioentes de transmisión y reflexión

n0 = math.sqrt( miu0/epsilon0)  # Impedancia intrinseca del aire

n1 = math.sqrt( miu0/(epsilon0*epsilonr1))  # Impedancia intrinseca del plastico

R1 = (n1 - n0)/(n1 + n0)  # Coeficiente de reflexión cuando va del aire al plastico

T1 = (2*n1)/(n1 + n0)    # Coeficiente de transmición cuando va del aire al plastico

R2 = (n0 - n1)/(n0 + n1)  # Coeficiente de reflexión cuando va del plastico al aire

T2 = (2*n0)/(n0 + n1)   # Coeficiente de transmición cuando va del plastico al aire

# Plot (cambiar)

# Animacion 
#fig = plt.figure(1)
#ax = fig.gca()

#def animation(n):
#    ax.clear()
#    ax.plot(X0[:n],Y0[:n])
#    plt.tittle(str(X0[n]))
  
#ani = animation.FuncAnimation(fig, animation, range (len(X0))

#plt.show()

# Grafica
X0,Y0 = np.meshgrid(X0, Y0)
#plt.plot(X0,Ez[0,:]) #plot pulse
fig1 = plt.figure()
f1 = plt.contourf(Ezt, 15, cmap = 'Set1')
plt.colorbar(f1)
plt.ylabel('Y(m)')
plt.xlabel('X(m)')
plt.show()