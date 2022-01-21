import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit as fit

"""
Auswertung FPraktikum Gaslaser, Fabry-Perot-Interferometer
"""

# Faserspektrum
# nu, I = np.loadtxt('Faserspec.txt', usecols=(0,1), unpack=True, skiprows=1)
# plt.plot(nu, I)

# Kommerzieller HeNe-Laser Modenbetrachtung

t1, I1 = np.loadtxt('mitPol.txt', usecols = (0,1), unpack=True, skiprows=4)

plt.plot(t1, I1)
