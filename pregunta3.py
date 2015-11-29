from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import (leastsq, curve_fit)
import random

''' Este script realiza el ajuste lineal a los flujos banda z y banda i
 ocupando monte carlo y polyfit,
 ademas de 2 histogramas para m y p de la recta que ajusta datos'''

np.random.seed(330)


def fit_datos(I, Z):
    coef1 = np.polyfit(I, Z, 1)
    coef2 = np.polyfit(Z, I, 1)
    x_c = (coef1[1]*coef2[0]+coef2[1])/(1.0-coef1[0]*coef2[0])
    y_c = coef1[0]*x_c + coef1[1]
    m = np.tan((np.arctan(coef1[0]) + np.arctan(1.0/coef2[0])) / 2.0)
    return (m, y_c - m*x_c)


# Lee los datos y realiza un ajuste lineal
bandI = 3.631*np.loadtxt("data/DR9Q.dat", usecols=(80,))
errI = 3.631*np.loadtxt("data/DR9Q.dat", usecols=(81,))
bandZ = 3.631*np.loadtxt("data/DR9Q.dat", usecols=(82,))
errZ = 3.631*np.loadtxt("data/DR9Q.dat", usecols=(83,))
coef = fit_datos(bandI, bandZ)


# Simulacion Monte-Carlo
N_mc = 10000
params = np.zeros((2, N_mc))
for i in range(N_mc):
    fake_bandI = np.random.normal(loc=bandI, scale=errI)
    # mean and standard desviation
    fake_bandZ = np.random.normal(loc=bandZ, scale=errZ)
    f_param = fit_datos(fake_bandI, fake_bandZ)
    # retorna m y posicion
    params[0][i] = f_param[0]
    params[1][i] = f_param[1]
m = np.sort(params[0])
p = np.sort(params[1])
m_025 = m[int(N_mc*0.025)]
p_025 = p[int(N_mc*0.025)]
m_975 = m[int(N_mc*0.975)]
p_975 = p[int(N_mc*0.975)]
print "Ajuste lineal entre flujo de banda i y banda z"
print "flujoz = m * flujoi + p"
print "m =", coef[0]
print "El Intervalo de confianza para m es: ", (m_025, m_975)
print "b =", coef[1], "[10^-6 Jy]"
print "El Intervalo de confianza para p es : ", (p_025, p_975), "[10^-6 Jy]"
xI = np.linspace(-20, 600, 500)
plt.figure(1)
plt.errorbar(bandI, bandZ, xerr=errI, yerr=errZ, fmt='g.', label='Medidas+Err')
plt.plot(xI, coef[0]*xI + coef[1], 'b', label='Mejor ajuste')
plt.plot(xI, m_975*xI + p_975, 'y--', label='Zona de confianza')
plt.legend(loc=2)
plt.plot(xI, m_025*xI + p_025, 'y--')
plt.title('Grafico de Flujo cuasares banda z v/s banda i')
plt.xlabel('Flujo en banda i [$10^{-6}$ Jy]')
plt.ylabel('Flujo en banda z [$10^{-6}$ Jy]')
plt.xlim([-20, 600])
plt.savefig('flujografico.png')

plt.figure(2)
plt.hist(m, bins=200, normed=True, label='Histograma de m')
plt.axvline(x=coef[0], label='Mejor ajuste', color='r')
plt.axvline(x=m_025, label='Zona de confianza', color='y')
plt.xlabel('$m$')
plt.ylabel('Frecuencia')
plt.title('Histograma de $m$ utilizando Monte-Carlo')
plt.legend()
plt.axvline(x=m_975, color='y')
plt.savefig('Histogramadem.png')

plt.figure(3)
plt.hist(p, bins=200, normed=True, label='Histograma de p')
plt.axvline(x=coef[1], label='Mejor ajuste', color='r')
plt.axvline(x=p_025, label='Zona de confianza', color='y')
plt.xlabel('$p \,[10^{-6} Jy]$')
plt.ylabel('Frecuencia')
plt.title('Histograma de $p$ utilizando Monte-Carlo')
plt.legend()
plt.axvline(x=p_975, color='y')
plt.savefig('Histogramadep.png')
plt.show()
