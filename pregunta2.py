from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import (leastsq, curve_fit)
import random

'''Este script realiza un ajuste lineal para los datos de
Freedman(v=Ho*d) y un Histograma de Ho'''


np.random.seed(330)


def f(H, x):
    return H * x


def bootstrap(datos):
    N_B = 100000
    H = np.zeros(N_B)
    N = len(datos)
    for i in range(N_B):
        s = np.random.randint(low=0, high=N, size=N)
        dist = datos[s, 1]
        vel = datos[s, 0]
        a1 = curve_fit(f, dist, vel)
        a2 = curve_fit(f, vel, dist)
        a = (a1[0] + 1/a2[0]) / 2
        H[i] = a
    H = np.sort(H)
    H_025 = H[int(0.025*N_B)]
    H_975 = H[int(0.975*N_B)]
    return (H_025, H_975, H)


datosSNI = np.loadtxt("data/SNIa.dat", usecols=(1, 2))
dSNI = np.zeros(len(datosSNI))
vSNI = np.zeros(len(datosSNI))

for i in range(len(datosSNI)):
    dSNI[i] = datosSNI[i][1]
    vSNI[i] = datosSNI[i][0]

# dist = data[:, 1]
# vel = data[:, 0]

a_1, a_covarianza_1 = curve_fit(f, dSNI, vSNI)
a_2, a_covarianza_2 = curve_fit(f, vSNI, dSNI)
a_prom = (a_1 + 1/a_2) / 2

N0 = len(dSNI)
H0 = a_prom
H0_conf = bootstrap(datosSNI)
print "El valor para la constante de Hubble es de Ho =", H0, "[km/(s*Mpc)]"
print "Intervalo de confianza: ", H0_conf[:2], "[km/(s*Mpc)]"

xS = np.linspace(0, 500, 50)
plt.figure(1)
plt.plot(dSNI, vSNI, 'ro', label='Mediciones')
plt.plot(xS, H0*xS, 'b', label='Mejor ajuste')
plt.plot(xS, H0_conf[0]*xS, 'y--', label='Zona de confianza')
plt.legend(loc=2)
plt.plot(xS, H0_conf[1]*xS, 'y--')
plt.title('Grafico Vel.recesion vs Distancia+Datos supernovas tipo I Freedman')
plt.xlabel('Distancia a las supernovas [Mpc]')
plt.ylabel('Velocidad de recesion de las nebulosas [km/s]')
plt.xlim([0, 500])
plt.savefig('vvsdFreedman.png')


plt.figure(2)
plt.hist(H0_conf[2], bins=100, normed=True, label='Histograma de H0')
plt.axvline(x=H0, label='Mejor ajuste', color='r')
plt.axvline(x=H0_conf[0], label='Zona de confianza', color='y')
plt.xlabel('$H_0$ [km/(s Mpc)]')
plt.title('Histograma de H0 usando bootstrap + Datos de Freedman 2000')
plt.legend()
plt.axvline(x=H0_conf[1], color='y')
plt.savefig('histFreedman.png')
plt.show()
