from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import (leastsq, curve_fit)
import random

'''Este script realiza un ajuste lineal para los datos de
Hubble(v=Ho*d) y un Histograma de Ho'''

np.random.seed(330)

datosHubble = np.loadtxt("data/hubble_original.dat")
dHubble = np.zeros(len(datosHubble))
vHubble = np.zeros(len(datosHubble))


for i in range(len(datosHubble)):
    dHubble[i] = datosHubble[i][0]
    vHubble[i] = datosHubble[i][1]


def f(H, x):
    return H * x


a_1, a_covarianza_1 = curve_fit(f, dHubble, vHubble)
a_2, a_covarianza_2 = curve_fit(f, vHubble, dHubble)



def bootstrap(datos):
    N_B=100000
    H = np.zeros(N_B)
    N=len(datos)
    for i in range(N_B):
        s = np.random.randint(low=0, high=N, size=N)
        dist = datos[s,0]
        vel = datos[s,1]
        a1 = curve_fit(f,dist, vel)
        a2 = curve_fit(f,vel, dist)
        a = (a1[0] + 1/a2[0]) / 2
        H[i] = a
    H = np.sort(H)
    H_025 = H[int(N_B*0.025)]
    H_975 = H[int(N_B*0.975)]
    return (H_025, H_975, H)



xH = np.linspace(0, 2.2, 50)
N0 = len(dHubble)
H0 = (a_1 + 1/a_2) / 2
H0_conf = bootstrap(datosHubble)
print "Datos usados por Edwin Hubble:"
print "El valor para la constante de hubble es de Ho =", H0, "[km/(s*Mpc)]"
print "El Intervalo de confianza es de : ", H0_conf[:2], "[km/(s*Mpc)]"


plt.figure(1)
plt.plot(dHubble, vHubble, 'ro', label='Mediciones originales')
plt.plot(xH, H0*xH, 'b', label='Mejor ajuste')
plt.plot(xH, H0_conf[0]*xH, 'y--', label='Zona de confianza')
plt.legend(loc=2)
plt.plot(xH, H0_conf[1]*xH, 'y--')
plt.title('Grafico Vel.recesion v/s Distancia + Datos originales de Edwin Hubble')
plt.ylabel('Velocidad de recesion de las nebulosas [km/s]')
plt.xlabel('Distancia a las nebulosas [Mpc]')
plt.xlim([0, 2.2])
plt.savefig('vvsdhubble.png')


plt.figure(2)
plt.hist(H0_conf[2], bins=80, normed=True, label='Histograma de H0')
plt.axvline(x=H0, label='Mejor ajuste', color='r')
plt.axvline(x=H0_conf[0], label='Zona de confianza', color='y')
plt.xlabel('$H_0$ [km/(s Mpc)]')
plt.title('Histograma de H0 usando bootstrap+ Datos originales de Hubble')
# lo mas frecuente sera el h0 que mejor se ajusta
plt.legend()
plt.axvline(x=H0_conf[1], color='y')
plt.savefig('histhubble.png')
plt.show()
