#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Thu Nov 18 10:51:11 2021

#@author: stevengebel
"""F-Praktikum Auswertung: Lebensdauer von Myonen
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as fit
# from Praktikum import niceplot
import matplotlib.ticker as ticker

plt.style.use('ggplot')
plt.rcParams.update({
    'font.family': 'serif',
    'text.usetex': False,
    'pgf.rcfonts': False,
})

"""Vorversuch"""
def linear(x, A, B):
    return A*x + B

def exponential(x, A, B, C, D):
    return A*np.exp(D*(x-C)) + B


# UPM1 in 50 V steps
U1 = np.array([1700, 1750, 1800, 1850, 1900, 1950, 2000,  2051, 2099, 2150,
               2200, 2251, 2301, 2350])
# UPM2
U2 = np.array([1657, 1710, 1755, 1799, 1835, 1878, 1908, 1942, 1974, 1997, 
               2016, 2024, 2037, 2039])
dU1 = np.array([1, 1, 2.5, 2.5, 5, 5,5, 5,5,10, 10, 10, 10, 10])
# Linearer Fit 
popt, cov = fit(linear, U1[:-2], U2[:-2], sigma = dU1[:-2], absolute_sigma=True)


# zweiter Vorversuch (bis 2200, weil danach nicht mehr lienar)
# extra measurements in U1 = [1900, 2000] -> 1925, 1975
U12 = np.array([1750, 1800, 1850, 1900, 1925, 1950, 1975, 2000, 2050, 2100, 
                2150, 2200])
U22 = linear(U12, *popt)
print(np.round(U22))
# Anzahl Ereignisse,  N>10000 is given because √N/N < 1% for U1=2000V
N = np.array([1280, 2718, 4634, 6583, 7101, 8415, 9100, 10017, 13250, 
              17919, 26795, 42167])
dN = np.sqrt(N)
print(np.round(dN))
Tmess = 435 # in s (soll konstant sein, s.d. für U1 = 2000V N>10000

"""Plot 1"""
# niceplot(x=U1, y=U2, error=True, xaxis=r'U$_{PM1}$/V', x2=U1[0:-2], y2=linear(U1[0:-2], *popt),
#         c2='tab:green', lw2=2.5, plot2=True, plotlabel2='fit to A$\cdot$U$_{PM1}$ + B',
#         yaxis=r'U$_{PM2}$/V', 
#         titel=r'Calibration curve for U$_{PM1}$ and U$_{PM2}$',
#         ls='', ls2='-.',c='k', legend=True, marker='.', yerr=dU1, lw=0, 
#         safename='vorversuch1', plotlabel='measurement points',
#         safefig=True)
print("Die Parameter des U2(U1)-Fits sind:")
print("A+-∆A={}+-{}".format(popt[0], np.sqrt(np.diag(cov))[0]))
print("B+-∆B={}+-{}".format(popt[1], np.sqrt(np.diag(cov))[1]))

#print(linear(U12, *popt))

"""Plot 2"""
# popt1, cov1= fit(linear, U12[0:3], N[0:3], sigma=dN[0:3], absolute_sigma=True)
# popt2, cov2 = fit(exponential, U12[7:-1], N[7:-1], p0=[5000, 10000, 2000, 0],
#                   sigma=dN[7:-1], absolute_sigma=True)
# print("Die Parameter des lin-Fits sind:")
# print("A+-∆A={}+-{}".format(popt1[0], np.sqrt(np.diag(cov1))[0]))
# print("B+-∆B={}+-{}".format(popt1[1], np.sqrt(np.diag(cov1))[1]))
# print("Die Parameter des exp-Fits sind:")
# print("A+-∆A={}+-{}".format(popt2[0], np.sqrt(np.diag(cov2))[0]))
# print("B+-∆B={}+-{}".format(popt2[1], np.sqrt(np.diag(cov2))[1]))
# print("C+-∆C={}+-{}".format(popt2[2], np.sqrt(np.diag(cov2))[2]))
# print("D+-∆D={}+-{}".format(popt2[3], np.sqrt(np.diag(cov2))[3]))
# niceplot(x=U12, y=N, yerr=dN, xaxis='U$_{PM1}$', 
#          yaxis=r'Coincidence count Z$_{12}$', error=True,
#           titel='Equal amplitude coincidence meauserement', ls='', safefig=True,
#           safename='vorversuch2', legend=True, plotlabel='measured data',
#           marker='.', c='k', plot2=True, plot3=True, c2='k', c3='k', 
#           x2=np.linspace(1750, 1900, 100), y2=linear(np.linspace(1750, 1900, 100), *popt1),
#           x3=np.linspace(2000, 2200, 100), y3=exponential(np.linspace(2000, 2200, 100), *popt2),
#           lw2=2.5,ls2='--', ls3='-.', lw3=2.5, plotlabel2=r'fit to A$\cdot$U$_{PM1}$+B', 
#           plotlabel3=r'fit to A$\cdot$ $\exp$((U$_{PM1}$-C)$\cdot$D)+B')


"""Kurzzeitversuch - Max-Likelihood-Methode"""

# channel = np.loadtxt('testLM_19.11.txt', usecols=0)
# N2 =  np.loadtxt('testLM_19.11.txt', usecols=1)

# # Berechnung der Halbwertszeit von Tau via max-likelihood-Methode: exp-Gesetz
# t1 =  (10**(-6)/24)*channel

# def tau1(x):
#     return (1/np.sum(N2[x:-1]))*np.sum(N2[x:-1]*t1[x:-1])

# def tau2(x):
#     return (1/np.sum(N2[25:x]))*np.sum(N2[25:x]*t1[25:x])

# intervall = np.arange(1,251)
# Tau1 = []
# Tau2 = []

# for i in intervall:
#     Tau1.append(tau1(i))
#     Tau2.append(tau2(i))


# plt.ylabel(r'N', size=20)
# plt.xlabel(r'channel', size=20)
# plt.title(r'',size=20)
# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# plt.plot(channel, N2, ms=ms, lw=0, marker='.', c='tab:blue')
# plt.show()

#plt.ylabel(r'$\tau$', size=20)
#plt.xlabel(r'initial channel', size=20)
#plt.title(r'',size=20)
#plt.ylim(0, 3e-6)
#plt.yticks(fontsize=18)
#plt.xticks(fontsize=18)
#plt.plot(intervall, Tau1, ms=ms, lw=0, marker='.', c='tab:blue')
#plt.plot(intervall, Tau2, ms=ms, lw=0, marker='.', c='tab:blue')
#plt.show()





