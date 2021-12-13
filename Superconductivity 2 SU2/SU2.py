#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 15:29:21 2021

@author: stevengebel
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as fit

# Daten einlesen
p = np.load('p.npy')
T = np.load('T.npy')
V_Ge = np.load('V_Ge.npy')
Trans = np.load('Transitions.npy')
V_Ge2 = np.loadtxt('high_T.txt', skiprows = 4, usecols = 4)
I_Ge2  = np.loadtxt('high_T.txt', skiprows = 3, usecols = 5)/10
T1 = np.loadtxt('high_T.txt', skiprows = 3, usecols = 1)

R_Ge2 = V_Ge2/I_Ge2 
R = np.mean(R_Ge2)
print(R)
t = 293.15
alpha = 0.208
    
print((R*np.pi*(0.25/2)**2)/((0.208)*(1+4.2e-3*(np.mean(T1)-t))))
###########################
plt.figure(figsize=(18,9))
lw = 3
ms = 13

Cfield = 23400*Trans[:,0]
uebergang = 23400*Trans[:,1]

# plt.ylabel(r'$\Delta_H \, / \, \frac{A}{m}$', size=20)
# plt.xlabel(r'T / K', size=20)
# # plt.title(r'$',size=20)
# # plt.ylim(,)
# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# plt.plot(T[:,0], uebergang, ms=ms, marker='.', lw=0)
# # plt.legend(fontsize=17, markerscale=1.5)
# plt.title(r'Übergangsbreite $\Delta_H$ in Abhängigkeit der Termperatur $T$',size=20)
# plt.savefig('uebergang.svg', dpi=300)
# plt.show()

I = [
np.loadtxt('Vorversuch1.txt', skiprows=3, usecols=5),
np.loadtxt('Vorversuch2_I_0.txt', skiprows=3, usecols=5),
np.loadtxt('Vorversuch3_I_0.txt', skiprows=3, usecols=5),
np.loadtxt('Vorversuch4_I_umgepolt_0.txt', skiprows=3, usecols=5),
np.loadtxt('Vorversuch5_umgepolt_I_0.txt', skiprows=3, usecols=5)]

mean_I = [np.mean(I[0]),np.mean(I[1]),np.mean(I[2]),np.mean(I[3]),np.mean(I[4])]
std_I = [np.std(I[0]),np.std(I[1]),np.std(I[2]),np.std(I[3]),np.std(I[4])]

# für Ge-therm (ersten beiden Temperaturen abgeändert, s.d. fit besser an 
# gemessene daten passt, in R(T) und T(R) sichtbar, vgl. mit Dampfdrucktabelle)
Tmod = np.array([1.65, 2.15, 2.20218519, 2.39051667, 2.58537097, 2.90315517,
 3.16122581, 3.45738028, 3.75964583, 4.02471875, 4.23158462, 4.50267213,
 4.80481538, 5.09777333, 5.39795,    5.69574,    6.00012308, 6.29858,
 6.59711667])

R_Ge = V_Ge[:,0]/100e-6
# Messunsicherheit für Widerstand (Ge-fit)
dR_Ge = R_Ge*(V_Ge[:,1]/V_Ge[:,0])

###################################
# Aus der dampfdruckkurve bestimmte Temperatur (oft annäherung)
Tdampf = [1.65, 2.15,2.2,2.4,2.6,3.0,3.3,3.65,3.95,4.25]
# wirklich gemessene
p_tabelle = [7.0476, 35.878, 40.576, 63.476, 93.988,175.57,
             277.63,425.28,587.63,787.43]

# diff = p_tabelle-p[:,0][0:len(p_tabelle)]
# print(p_tabelle, p[:,0])
# # Plot
# plt.ylabel(r'$(p_{Tabelle} - p_{Messung})/p_{Tabelle}$', size=20)
# plt.xlabel(r'$p_{Tabelle}$', size=20)
# # plt.title(r'$',size=20)
# # plt.ylim(,)
# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# plt.plot(p_tabelle, diff/p_tabelle, ms=ms, marker='.', lw=0, 
#           label=r'$(p_{Tabelle} - p_{Messung})/p_{Tabelle}$')
# plt.legend(fontsize=17, markerscale=1.5)
# plt.savefig('pgueteohnebetrag.svg', dpi=300)
# plt.show()


##############
# pdampf = np.array([3.0688, 4.5660, 5.7054, 7.0476, 8.6136, 10.424, 12.500, 14.861, 
#                    17.526, 20.511, 23.832, 27.498, 31.514, 35.878, 40.576, 45.655, 51.151,
#                    57.085, 63.476, 70.343, 77.704, 85.579, 93.988, 102.95, 112.48, 122.59,
#                    133.31, 144.65, 156.63, 169.26, 175.51, 196.57, 211.28, 
#                    226.73, 242.92, 259.89, 277.63, 296.18, 315.55, 335.76, 356.81, 
#                    378.74, 401.56, 425.28, 449.92, 475.51, 502.05, 529.57, 558.09,
#                    587.63, 618.21, 649.86, 682.59, 716.42, 751.37, 787.43, 824.65])

# Schränke Intervall zum fitten bis zum Phasenübergang bei T~4.2K ein
# print(T[:,0][0:11])
# parameters4, cov4 = fit(dampf_fit, p[:,0][0:11], T[:,0][0:11], 
#                         absolute_sigma=True)
# # plt.plot(np.linspace(1.5,980,1000), dampf_fit(np.linspace(1.5,980,1000), *parameters4), 
# #           c='tab:blue', lw=lw, 
# #           label = r'Fitfunktion T(p) = $\sqrt{A \cdot p}$ + B')
# # plt.plot(p[:,0], T[:,0], ms=ms, marker='*', lw=0, c='tab:orange', label = 'Messdaten')
# # plt.legend(fontsize=20, markerscale=2)
# # plt.xlabel(r'p/mmHg', size=20)
# # plt.ylabel(r'T/K', size=20)
# # plt.savefig('pdampf_messung.png', dpi=300)
# # plt.show()
# print("Die Parameter des T(p)-Fits sind A+-∆A={}+-{}, B+-∆B = {}+-{}".format(parameters4[0],
#       np.sqrt(np.diag(cov4))[0], parameters4[1], np.sqrt(np.diag(cov4))[1]))
# # Abweichung von gemessenen Daten
# print("Die Abweichung von den gemessenen Daten und denen des Fitmodells für |T(p)_fit - T(p)_mess|:")
# print(np.abs(T[:,0][0:11] - dampf_fit(p[:,0][0:11], *parameters4)))




# parameters3, cov3 = fit(dampf_fit, pdampf, Tdampf, absolute_sigma=True)
# plt.xlabel(r'p/mmHg', size=20)
# plt.ylabel(r'T/K', size=20)
# plt.title(r'Dampfdruckkurve von $^{4}He$',size=20)
# # plt.ylim(,)
# plt.yticks(fontsize=14)
# plt.xticks(fontsize=14)
# plt.plot(np.linspace(1.5,800,1000), dampf_fit(np.linspace(1.5,800,1000), *parameters3), 
#           c='tab:blue', lw=lw, 
#           label = r'Fitfunktion T(p) = $A\sqrt{B \cdot p}$ + C')
# plt.plot(pdampf, Tdampf, ms=ms, c='tab:orange', marker='.', label='Dampfdrucktabelle')
# plt.legend(fontsize=20, markerscale=2)
# plt.savefig('pdampf4png', dpi=300)

# print("Die Parameter des pdampf-Fits sind \
#       A+-∆A={}+-{}, B+-∆B={}+-{}, C+-∆C = {}+-{}".format(parameters3[0],
#       np.sqrt(np.diag(cov3))[0], parameters3[1], np.sqrt(np.diag(cov3))[1], 
#     parameters3[2], np.sqrt(np.diag(cov3))[2]))
# plt.show()

# angepasste Temperatur für Ge-Plot (Übernahme von Werten aus Dampfdrucktabelle für ^4He)
# Modelle, wie A*ln(B*T) 

# Fitfunktion für kritisches Magnetfeld
def Hc_fit(T, H_c, H_c0, T_c):
    """
    Parameters
    ----------
    T : array
        Temperaturen von 4.2 K bis 7 K.
    H_c : array
        kritische Magnetfelder.
    H_c0 : integer
        Amplitude, Fitparameter.
    T_c : integer
        kritische Temperatur, Fitparameter.
    """
    return H_c0*(1-(T/T_c)**2)

# Fitfunktion für Ge-Thermometer: T(R)
def Ge_fit(R, A, B):
    
    """
    Parameters
    ----------
    T : array
        Temperatur.
    R : array
        Widerstand.
    """
    return -1/B*np.log(R/A)

def Ge_fit2(T, A, B):
    return A*np.exp(-B*T)

###################
# sigma = T[:,1], absolute_sigma=True
# Fit-Daten erhalten
# parameters1, cov1 = fit(Ge_fit2, Tmod, R_Ge, 
#                                 absolute_sigma=True)
# parameters2, cov2 = fit(Ge_fit, R_Ge, Tmod, sigma = T[:,0], absolute_sigma=False)
# parameters3, cov3 = fit(Ge_fit2, R_Ge, T[:,0], 
#                                 absolute_sigma=True)
# parameters4, cov4 = fit(Ge_fit, R_Ge, T[:,0], 
#                                 absolute_sigma=True)

# print("Die Parameter des R(T)-Fits sind:")
# print("A+-∆A={}+-{}".format(parameters1[0], np.sqrt(np.diag(cov1))[0]))
# print("B+-∆B={}+-{}".format(parameters1[1], np.sqrt(np.diag(cov1))[1]))
# print("Die Parameter des T(R)-Fits sind:")
# print("A+-∆A={}+-{}".format(parameters2[0], np.sqrt(np.diag(cov2))[0]))
# print("B+-∆B={}+-{}".format(parameters2[1], np.sqrt(np.diag(cov2))[1]))
# print("C+-∆C={}+-{}".format(parameters2[2], np.sqrt(np.diag(cov2))[2]))



# ##############
# # Plot-Bereich
# T vs. R
# plt.xlabel(r'R/$\Omega$', size=20)
# plt.ylabel(r'T/K', size=20)
# plt.title(r'Temperatur-Widerstandskurve des Ge-Thermometers',size=20)
# plt.ylim(min(Tmod) - 0.3, max(Tmod) + 0.3)
# plt.xlim(0, 8500)
# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# # # plt.plot(np.linspace(min(R_Ge), max(R_Ge), 1000), 
# # #           Ge_fit(np.linspace(min(R_Ge), max(R_Ge), 1000), *parameters4), 
# # #           c='tab:red', lw=1.5, linestyle='--',
# # #           label = r'Fitfunktion ohne Korrektur$')
# plt.plot(np.linspace(min(R_Ge), max(R_Ge), 1000), 
#           Ge_fit(np.linspace(min(R_Ge), max(R_Ge), 1000), *parameters2), 
#           c='tab:blue', lw=lw, 
#           label = r'Fitfunktion T(R) = $-\frac{1}{B} \cdot \ln \left ( \frac{R}{A} \right)$')
# # # plt.plot(np.linspace(min(T[:,0]), max(T[:,0]), 1000), 
# # #           Ge_fit2(np.linspace(min(T[:,0]), max(T[:,0]), 1000), *parameters1), 
# # #           c='tab:blue', lw=lw, 
# # #           label = r'Fitfunktion R(T) = $A\cdot exp(-B\cdot T)$')
# plt.plot(R_Ge, T[:,0], ms=ms, c='tab:orange', marker='.', label='Messdaten', lw=0)
# # # plt.errorbar(R_Ge, Tmod, yerr=T[:,1], 
# # #              capsize = 10, c='tab:orange') VIEL ZU KLEINE FEHLER
# plt.legend(fontsize=15, markerscale=2)
# plt.savefig('Ge_T(R)FINAL.svg', dpi=300)
# plt.show()

# plt.xlabel(r'T/K ', size=20)
# plt.ylabel(r'R/$\Omega$', size=20)
# plt.title(r'Temperatur-Widerstandskurve des Ge-Thermometers',size=20)
# plt.ylim(0, 9000)
# plt.xlim(min(T[:,0]) - 0.1, max(T[:,0]) + 0.1)
# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# plt.plot(np.linspace(min(T[:,0]), max(T[:,0]), 1000), 
#           Ge_fit2(np.linspace(min(T[:,0]), max(T[:,0]), 1000), *parameters1), 
#           c='tab:blue', lw=lw, 
#           label = r'Fitfunktion R(T) = $A\cdot \exp(-B\cdot T)$')
# plt.plot(np.linspace(min(R_Ge), max(R_Ge), 1000), 
#           Ge_fit2(np.linspace(min(R_Ge), max(R_Ge), 1000), *parameters3), 
#           c='tab:red', lw=1.5, linestyle='--',
#           label = r'Fitfunktion ohne Korrektur$')
# plt.plot(T[:,0], R_Ge, ms=ms, c='tab:orange', marker='.', label='Messdaten',  lw=0)
# # plt.errorbar(R_Ge, Tmod, yerr=T[:,1], 
# #              capsize = 10, c='tab:orange') VIEL ZU KLEINE FEHLER
# plt.legend(fontsize=15, markerscale=2)
# plt.savefig('Ge_R(T)2.svg', dpi=300)
# plt.show()



##########################
#????
# T_fit(R)- T vs. T(p)
# print(R_Ge, p[:,0], T[:,0][0:11])
# plt.ylabel(r'$|(T_{fit}(p) - T_i)|$/K', size=20)
# plt.xlabel(r'p/mmHg', size=20)
# plt.title(r'',size=20)
# plt.ylim()
# plt.xlim()
# plt.yticks(fontsize=14)
# plt.xticks(fontsize=14)
# plt.plot(p[:,0][0:11], np.abs(T[:,0][0:11] - dampf_fit(p[:,0][0:11], *parameters4)),c='tab:blue', lw=0, 
#           label = r'$|T(p)_{fit} - T(p)_{mess}|$', marker='*', ms=ms)
# plt.plot(pdampf, np.abs(dampf_fit(pdampf, *parameters4) - Tdampf), c='tab:orange', lw=0, 
#           label = r'$|T(p)_{fit} - T(p)_{Tabelle}|$', marker='*', ms=ms)
# plt.errorbar(T[:,0], R_Ge, yerr=dR_Ge, capsize = 10, c='tab:orange')
# plt.legend(fontsize=20, markerscale=2)
# plt.savefig('T_fit.png', dpi=300)
# plt.show()


##########################
# T(R) - T(p) vs. T(p) T(R)-T_eich vs. T(R)

# # T(R) - T(p) vs. T(p)
# plt.ylabel(r'$(T_{fit}(R) - T(p))$/K', size=20)
# plt.xlabel(r'$T(p)/K, \,\, T_{fit}(R)$/K', size=20)
# plt.title('Abweichungen der jeweiligen Temperaturmessmethoden',size=20)
# # plt.ylim(0.025, 0.27)
# # plt.xlim(min(dampf_fit(p[:,0][0:11], *parameters4)), max(dampf_fit(p[:,0][0:11], *parameters4))+0.3)
# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# # T(R) - T(p)
# plt.plot(Tdampf, 
#           (Ge_fit(R_Ge[0:len(p_tabelle)], *parameters2) - Tdampf)/Tdampf,
#           c='tab:blue', lw=0, 
#           label = r'$(T(R)_{fit} - T(p))/T(p)$', marker='.', ms=ms)
# # T(R) - T_eich(R)
# T_eich = [2,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.59,3.81,4.0,4.2,4.5,5.0,5.5]
# R0 = np.array([8836.21, 6489.62, 4983.43, 3953.46, 3205.12, 2769.28, 2335.09, 
#                 1976.83, 1731.62, 1524.65, 2348.98, 1186.75, 1006.24, 770.88, 
#                 606.35])
# plt.plot(T_eich, (T_eich - Ge_fit(R0, *parameters2))/T_eich, c='tab:orange', 
#           marker='.', ms=ms, label=r'$(T(R)_{fit} - T(R)_{Eich})/T(R)_{Eich}$',lw=0)
# # plt.errorbar(T[:,0], R_Ge, yerr=dR_Ge, capsize = 10, c='tab:orange')
# plt.legend(fontsize=17, markerscale=1.5)
# plt.savefig('TvsT12relohnebetrag.svg', dpi=300)
# plt.show()

####################
# print(mean_I)
# plt.ylabel(r'$H$/ T', size=20)
# plt.xlabel(r'$I_{Probe}$/ A', size=20)
# plt.title(r'Abhängigkeit de kritischen Magnetfelds von dem Probenstrom',size=20)
# # plt.ylim()
# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# plt.plot(mean_I, Trans[:,0], ms=ms, marker='.', lw=0, c = 'tab:orange')
# plt.errorbar(mean_I, Trans[:,0], xerr=std_I, capsize = 10, c='tab:orange', lw=0)
# plt.legend(fontsize=17, markerscale=1.5)
# plt.savefig('HvsI.svg', dpi=300)
# plt.show()





