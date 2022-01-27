#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Fri Jan 21 18:31:07 2022
# @author: stevengebel

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
from scipy.signal import find_peaks

"""
Auswertung FPraktikum Gaslaser, Fabry-Perot-Interferometer
"""


def niceplot(x, y, c, lw=None, lw2=None, lw3=None, lw4=None, lw5=None, lw6=None,
             lw7=None, ls=None, ls6=None, ls7=None,
             ls2=None, ls3=None, ls4=None, ls5=None, plot2=None, plot3=None,
             plot4=None, plot5=None, plot6=None, plot7=None, x2=None, c2=None,
             y2=None, x3=None,
             y3=None, c3=None, x4=None, y4=None, x5=None, y5=None, c4=None,
             x6=None, y6=None, c6=None, x7=None, y7=None, c7=None,
             c5=None, marker=None, marker2=None, marker3=None, marker4=None,
             marker5=None, marker6=None, marker7=None, ylabel=None, xlabel=None,
             ms=10, cs=5, fs=15, ticks=6,
             size=(8, 8), safefig=None, error=None, errorlabel=None,
             safename=None, yaxis=None, xaxis=None, yerr=None, plotlabel=None,
             legend=None, plotlabel2=None, plotlabel3=None, plotlabel4=None,
             plotlabel5=None, plotlabel6=None, plotlabel7=None, titel=None,
             # xlim=None, ylim=None
             ):
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel(r'{}'.format(yaxis), size=fs)
    ax.set_xlabel(r'{}'.format(xaxis), size=fs)
    ax.set_title(r'{}'.format(titel), size=fs+3)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    #####
    # pinmp_ticks(ax.xaxis, ticks) # For special plotstyle
    # pinmp_ticks(ax.yaxis, ticks)
    # ax.grid(which='minor', alpha=.5)
    # ax.grid(which='major', alpha=.5)
    # ax.tick_params(which='both')
    #####
    # ax.set_ylim(ylim)
    # ax.set_xlim(xlim)

    if error == True:
        ax.errorbar(x, y, yerr=yerr, capsize=cs, c=c, ms=ms, ls=ls,
                    marker=marker, label=r'{}'.format(plotlabel))
    else:
        ax.plot(x, y, ms=ms, lw=lw, ls=ls, marker=marker, c=c,
                label=r'{}'.format(plotlabel))
    # Bei Bedarf mehrere Fehlerbalken   
    # if error2 == True:
    #     ax.errorbar(x, y, yerr=yerr2, capsize=cs, c=c2,
    #                  label = r'${}$'.format(errorlabel2))
    # if error3 == True:
    #     ax.errorbar(x, y, yerr=yerr3, capsize=cs, c=c3,
    #                  label = r'${}$'.format(errorlabel3))
    # if error4 == True:
    #     ax.errorbar(x, y, yerr=yerr4, capsize=cs, c=c4,
    #                  label = r'${}$'.format(errorlabel4))
    if plot2 == True:
        ax.plot(x2, y2, ms=ms, lw=lw2, ls=ls2, marker=marker2, c=c2,
                label=r'{}'.format(plotlabel2))
    if plot3 == True:
        ax.plot(x3, y3, ms=ms, lw=lw3, ls=ls3, marker=marker3, c=c3,
                label=r'{}'.format(plotlabel3))
    if plot4 == True:
        ax.plot(x4, y4, ms=ms, lw=lw4, ls=ls4, marker=marker4, c=c4,
                label=r'{}'.format(plotlabel4))
    if plot5 == True:
        ax.plot(x5, y5, ms=ms, lw=lw5, ls=ls5, marker=marker5, c=c5,
                label=r'{}'.format(plotlabel5))
    if plot6 == True:
        ax.plot(x6, y6, ms=ms, lw=lw6, ls=ls6, marker=marker6, c=c6,
                label=r'{}'.format(plotlabel6))
    if plot7 == True:
        ax.plot(x7, y7, ms=ms, lw=lw7, ls=ls7, marker=marker7, c=c7,
                label=r'{}'.format(plotlabel7))
    if legend == True:
        ax.legend(fontsize=fs, markerscale=ms/5, facecolor='white')
    if safefig == True:
        plt.savefig('{}.svg'.format(safename), dpi=300)

    plt.show()


# Faserspektrum
#nu, I = np.loadtxt('Faserspec.txt', usecols=(0,1), unpack=True, skiprows=1)
#plt.plot(nu, I)

"""
Hemisphärischer Resonator: Resonatorlänge
"""
p1 = np.array([1.23, 1.2, 1.17, 0.877, 128e-3, 1.32])  # in W
l1 = np.array([63.3, 69.8, 89.3, 95.6, 100.5, 79.5])  # in cm
sigma1 = np.array([2.8, 3.2, 3, 2.6, 1.8, 4.8]) * 1e-3
# niceplot(x=l1, y=p1, yerr=sigma1, c='tab:blue', ls='', marker='x',
        # titel='Einfluss der Resonatorlänge $L$', xaxis='L / cm', yaxis='P / mW',
        # error=True, size=(20, 10))
"""
Hemisphärischer Resonator: Ausgangsleistung mit Polarisator
"""
phi = np.array([96, 106, 116, 86, 76, 56, 36, 26, 16, 6, -4, -14])
p2 = np.array(
    [2.464e-6, 30.9e-6, 110e-6, 25.3e-6, 102e-6, 0.399e-3, 0.744e-3, 0.858e-3, 0.972e-3, 1.02e-3, 0.994e-3, 0.907e-3])
sigma2 = np.array([24e-9, 0.117e-6, 0.6e-6, 0.19e-6, 6.54e-6, 1.1e-3, 3.7e-3, 3.7e-3, 3 - 2e-3, 2.4e-3, 2.3e-3, 2.0e-3])
# plt.plot(phi, p2, ls='', marker='x')
# niceplot(x=phi, y=p2, yerr=sigma2, c='tab:blue', error=True, size=(20,10))

"""
Kommerzieller HeNe-Laser Modenbetrachtung
"""
def gaussian(x, A, mu, sigma):
    """Gauss-function with arbitrary amplitude A (not normed).
    Parameters
    ----------
    x : array.
    A : integer, amplitude.
    mu : integer, expectation value.
    sigma : integer, standard deviation.
    B: noise
    """
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

# Daten einlesen
t1, I1 = np.loadtxt('mitPol.txt', usecols=(0, 1), unpack=True, skiprows=4)
t2, I2 = np.loadtxt('ohnePol.txt', usecols=(0, 1), unpack=True, skiprows=4)
# get peaks for gaussian fit
indices1, parameters1 = find_peaks(I1, height=50)
indices2, parameters2 = find_peaks(I2, height=50)
# niceplot(x=t1[indices1], y=I1[indices1],
#          x2=t1, y2=I1, c='tab:red', c2='tab:blue',
#          plot2=True, ls='', marker='x', plotlabel='peaks', plotlabel2='Polarisationswinkel $\Theta=\pi/2$',
#          titel='Modenbetrachtung des kommerziellen HeNe-Lasers', legend=True,
#          xaxis='time', yaxis='intensity', size=(10,5))
# niceplot(x=t2[indices2], y=I2[indices2],
#          x2=t2, y2=I2, c='tab:red', c2='tab:blue',
#          plot2=True, ls='', marker='x', plotlabel='peaks', plotlabel2='Polarisationswinkel $\Theta=2\pi$',
#          titel='Modenbetrachtung des kommerziellen HeNe-Lasers', legend=True,
#          xaxis='time', yaxis='intensity', size=(10,5))

"""
Kommerzieller HeNe-Laser Modenbetrachtung, versch. Resonatorlänge
"""

# Daten einlesen
t3, I3 = np.loadtxt('4peak.txt', usecols=(0, 1), unpack=True, skiprows=5) # schlechte Peaks
t4, I4 = np.loadtxt('4peak2.txt', usecols=(0, 1), unpack=True, skiprows=5) # Bessere Peaks
t5, I5 = np.loadtxt('4peak61.txt', usecols=(0, 1), unpack=True, skiprows=5)
# get peaks for gaussian fit
# indices3, parameters3 = find_peaks(I3, height=25) schlechte daten
indices4, parameters4 = find_peaks(I4, height=40, distance=16)
print("i = {}".format(indices4))
print("i2 = {}".format(indices4[0:4]))
print("i3 = {}".format(indices4[4:8]))
indices5, parameters5 = find_peaks(I5, height=40, distance=16)
# Gauss-fits
# popt3, cov3 = fit(gaussian, t3[indices3], I3[indices3] # schlechte daten
                 # )
popt41, cov41 = fit(gaussian, t4[indices4[0:4]], I4[indices4[0:4]], p0=[120, 130, 1e3],
                   # bounds=[[1e2, 1e2, 1e3, 1],[1.5e2, 3e2, 1e4, 100]]
                 )
popt42, cov42 = fit(gaussian, t4[indices4[4:8]], I4[indices4[4:8]], p0=[120, 400, 1e3],
                   # bounds=[[1e2, 1e2, 1e3, 1],[1.5e2, 8e2, 1e4, 100]]
                 )
popt51, cov51 = fit(gaussian, t5[indices5[0:4]], I5[indices5[0:4]], p0=[120, 130, 1e3]
                 )
popt52, cov52 = fit(gaussian, t5[indices5[4:8]], I5[indices5[4:8]], p0=[120, 130, 1e3]
                 )
print("GAUSSIAN PARAMETERS41:")
print("A+-∆A={}+-{}".format(popt41[0], np.sqrt(np.diag(cov41))[0]))
print("mu+-∆mu={}+-{}".format(popt41[1], np.sqrt(np.diag(cov41))[1]))
print("sigma+-∆sigma={}+-{}".format(popt41[2], np.sqrt(np.diag(cov41))[2]))
# print("B+-∆B={}+-{}".format(popt41[3], np.sqrt(np.diag(cov41))[3]))
print("The gaussian-line-with is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt41[2],
                                                               2 * np.sqrt(2 * np.log(2)) * popt41[2] *
                                                               np.sqrt(np.diag(cov41))[2]))
print("GAUSSIAN PARAMETERS42:")
print("A+-∆A={}+-{}".format(popt42[0], np.sqrt(np.diag(cov42))[0]))
print("mu+-∆mu={}+-{}".format(popt42[1], np.sqrt(np.diag(cov42))[1]))
print("sigma+-∆sigma={}+-{}".format(popt42[2], np.sqrt(np.diag(cov42))[2]))
# print("B+-∆B={}+-{}".format(popt42[3], np.sqrt(np.diag(cov42))[3]))
print("The gaussian-line-with is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt42[2],
                                                               2 * np.sqrt(2 * np.log(2)) * popt42[2] *
                                                               np.sqrt(np.diag(cov42))[2]))

print("GAUSSIAN PARAMETERS51:")
print("A+-∆A={}+-{}".format(popt51[0], np.sqrt(np.diag(cov51))[0]))
print("mu+-∆mu={}+-{}".format(popt51[1], np.sqrt(np.diag(cov51))[1]))
print("sigma+-∆sigma={}+-{}".format(popt51[2], np.sqrt(np.diag(cov51))[2]))
print("The gaussian-line-with is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt51[2],
                                                               2 * np.sqrt(2 * np.log(2)) * popt51[2] *
                                                               np.sqrt(np.diag(cov51))[2]))
print("GAUSSIAN PARAMETERS52:")
print("A+-∆A={}+-{}".format(popt52[0], np.sqrt(np.diag(cov52))[0]))
print("mu+-∆mu={}+-{}".format(popt52[1], np.sqrt(np.diag(cov52))[1]))
print("sigma+-∆sigma={}+-{}".format(popt52[2], np.sqrt(np.diag(cov52))[2]))
print("The gaussian-line-with is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt52[2],
                                                               2 * np.sqrt(2 * np.log(2)) * popt52[2] *
                                                               np.sqrt(np.diag(cov52))[2]))
# Plot
niceplot(x=t4, x2=t3, y=I4, y2=I3, plotlabel=r'L=80cm second image',
         plotlabel2='L=80cm ', plotlabel3='L=61cm',
         plot2=False, size=(10, 5), c='tab:blue', c2='tab:olive', legend=True, plot4=True,
         plot5=True, x5=t4[indices4], y5=I4[indices4], c5='tab:red', ls5='', marker5='x', plotlabel5='peaks',
         x4=np.linspace(0,300, 1000), y4=gaussian(np.linspace(0,300,1000), *popt41), c4='tab:green',
         x6=np.linspace(300,600,1000), plot6=True, c6='tab:pink', y6=gaussian(np.linspace(300,600,1000), *popt42),
         plotlabel6='gaussian fit'
         )

niceplot(x=t5, y=I5, c='tab:blue', legend=True, plotlabel='L=61cm',
         x2=np.linspace(0,300,1000), y2=gaussian(np.linspace(0,300,1000), *popt51), plot2=True, c2='tab:green',
         x3=t5[indices5], y3=I5[indices5], plot3=True, ls3='', marker3='x', c3='tab:red',
         plot4=True, x4=np.linspace(300,600,1000), y4=gaussian(np.linspace(300,600,1000), *popt52) , c4='tab:pink',
         plotlabel4='gaussian fit', size=(10, 5)
         )