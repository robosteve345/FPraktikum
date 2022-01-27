#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Fri Jan 21 18:31:07 2022
# @author: stevengebel

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit as fit

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
             ms=10, cs=5, fs=20, ticks=6,
             size=(8, 8), safefig=None, error=None, errorlabel=None,
             safename=None, yaxis=None, xaxis=None, yerr=None, plotlabel=None,
             legend=None, plotlabel2=None, plotlabel3=None, plotlabel4=None,
             plotlabel5=None, plotlabel6=None, plotlabel7=None, titel=None,
             # xlim=None, ylim=None
             ):
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel(r'{}'.format(yaxis), size=22)
    ax.set_xlabel(r'{}'.format(xaxis), size=22)
    ax.set_title(r'{}'.format(titel), size=25)
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
        ax.legend(fontsize=fs, markerscale=ms - 8, facecolor='white')
    if safefig == True:
        plt.savefig('{}.svg'.format(safename), dpi=300)

    plt.show()


# Faserspektrum
# nu, I = np.loadtxt('Faserspec.txt', usecols=(0,1), unpack=True, skiprows=1)
# plt.plot(nu, I)

"""
Hemisphörischer Resonator: Resonatorlänge
"""
p1 = np.array([1.23, 1.2, 1.17, 0.877, 128e-3, 1.32])  # in W
l1 = np.array([63.3, 69.8, 89.3, 95.6, 100.5, 79.5])  # in cm
sigma1 = np.array([2.8, 3.2, 3, 2.6, 1.8, 4.8]) * 1e-3
niceplot(x=l1, y=p1, yerr=sigma1, c='tab:blue', ls='', marker='x',
         titel='Einfluss der Resonatorlänge $L$', xaxis='L / cm', yaxis='P / mW',
         error=True)
"""
Hemisphörischer Resonator: Ausgangsleistung mit Polarisator
"""
phi = np.array([96, 106, 116, 86, 76, 56, 36, 26, 16, 6, -4, -14])
p2 = np.array(
    [2.464e-6, 30.9e-6, 110e-6, 25.3e-6, 102e-6, 0.399e-3, 0.744e-3, 0.858e-3, 0.972e-3, 1.02e-3, 0.994e-3, 0.907e-3])
sigma2 = np.array([24e-9, 0.117e-6, 0.6e-6, 0.19e-6, 6.54e-6, 1.1e-3, 3.7e-3, 3.7e-3, 3 - 2e-3, 2.4e-3, 2.3e-3, 2.0e-3])
plt.plot(phi, p2, ls='', marker='x')
niceplot(x=phi, y=p2, yerr=sigma2, c='tab:blue', error=True, size=(10,10))

"""
Kommerzieller HeNe-Laser Modenbetrachtung
"""
t1, I1 = np.loadtxt('mitPol.txt', usecols=(0, 1), unpack=True, skiprows=4)
t2, I2 = np.loadtxt('ohnePol.txt', usecols=(0, 1), unpack=True, skiprows=4)
niceplot(x=t1, y=I1, x2=t2, y2=I2, plotlabel='Mit Polarisierung', safefig=False,
         c='tab:blue', c2='tab:orange', plot2=True, legend=True,
         plotlabel2='Ohne Polarisierung', titel='')

"""
Kommerzieller HeNe-Laser Modenbetrachtung, versch. Resonatorlänge
"""
t3, I3 = np.loadtxt('4peak.txt', usecols=(0, 1), unpack=True, skiprows=5)
t4, I4 = np.loadtxt('4peak2.txt', usecols=(0, 1), unpack=True, skiprows=5)
t5, I5 = np.loadtxt('4peak61.txt', usecols=(0, 1), unpack=True, skiprows=5)

niceplot(x=t3, x2=t4, x3=t5, y=I3, y2=I4, y3=I5, plotlabel=r'L=80cm',
         plotlabel2='L=80 cm DIFFERENT PICTURE', plotlabel3='L=61cm',
         plot2=True, plot3=True, size=(20, 10), c='tab:blue', c2='tab:red',
         c3='tab:orange', legend=True)