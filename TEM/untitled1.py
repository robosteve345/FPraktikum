#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Sat Dec  4 17:21:16 2021

#@author: stevengebel

"""F-Praktikum Auswertung:Transmissionselektronenmikroskop
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as fit

x, y =np.loadtxt('02122021goldkrass.txt', skiprows=1, usecols=(0,1), unpack=True)

plt.plot(x, y, markersize=3, marker='.', ls='-.')

tau=np.array([2.191e-6, 2.191e-6, 2.201e-6])/(6.58212e-16) # in sekunden

m = 1.88353162742e-28/(1.78266e-36)  # in kg

y = np.sqrt((192*np.pi**3)/(tau*m**5))

v = np.sqrt(1/(np.sqrt(2)*y))

dtau = np.array([0.06e-6, 0.034e-6, 0.034e-6])/(6.58212e-16)

dy = 0.5*np.sqrt((192*np.pi**3)/(m**5*tau**3))*dtau

dv = dy/(2**(5/4)*y**(3/2))

# print(y, v, dy, dv, dv/v)
print(14192+41146)