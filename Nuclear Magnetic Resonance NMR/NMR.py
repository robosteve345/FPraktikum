#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Thu Dec 16 11:20:08 2021

#@author: stevengebel

"""F-Praktikum Auswertung NMR"""


import numpy as np
import matplotlib.pyplot as plt
import glob 
from scipy.optimize import curve_fit as fit


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
             size=(8,8), safefig=None, error=None, errorlabel=None, 
             safename=None, yaxis=None, xaxis=None, yerr=None, plotlabel=None, 
             legend=None, plotlabel2=None, plotlabel3=None, plotlabel4=None,
             plotlabel5=None, plotlabel6=None, plotlabel7=None, titel=None):
    
    fig = plt.figure(figsize=size) 
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel(r'{}'.format(yaxis), size=22)
    ax.set_xlabel(r'{}'.format(xaxis), size=22)
    ax.set_title(r'{}'.format(titel),size=25)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    # ax.set_yticklabels(x, fontsize=fs)
    # ax.set_xticklabels(y, fontsize=fs)
    #####
    # pinmp_ticks(ax.xaxis, ticks)
    # pinmp_ticks(ax.yaxis, ticks)
    # ax.grid(which='minor', alpha=.5)
    # ax.grid(which='major', alpha=.5)
    # ax.tick_params(which='both')
    #####
    ########## 
    # For NMR:
    ax.set_ylim(0, max(y)+0.05*max(y))
    ax.set_xlim(-100, 100)
    ##########
    
    if error == True:
        ax.errorbar(x, y, xerr=yerr, capsize=cs, c=c, ms=ms, ls=ls,
                     marker=marker, label = r'{}'.format(plotlabel))
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
    if  plot2 == True:
        ax.plot(x2, y2, ms=ms, lw=lw2, ls=ls2, marker=marker2, c=c2, 
             label=r'{}'.format(plotlabel2))
    if  plot3 == True:
        ax.plot(x3, y3, ms=ms, lw=lw3, ls=ls3, marker=marker3, c=c3, 
             label=r'{}'.format(plotlabel3))
    if  plot4 == True:
        ax.plot(x4, y4, ms=ms, lw=lw4, ls=ls4, marker=marker4, c=c4, 
             label=r'{}'.format(plotlabel4))
    if  plot5 == True:
        ax.plot(x5, y5, ms=ms, lw=lw5, ls=ls5, marker=marker5, c=c5, 
             label=r'{}'.format(plotlabel5)) 
    if  plot6 == True:
        ax.plot(x6, y6, ms=ms, lw=lw6, ls=ls6, marker=marker6, c=c6,
                label=r'{}'.format(plotlabel6))   
    if  plot7 == True:
        ax.plot(x7, y7, ms=ms, lw=lw7, ls=ls7, marker=marker7, c=c7, 
                label=r'{}'.format(plotlabel7))
    if legend == True:
        ax.legend(fontsize=fs-5, markerscale=ms-2, facecolor='white')
    if safefig == True:
        plt.savefig('{}.svg'.format(safename), dpi=300)

    plt.show()


"""Part 1: Centre frequency and parameter estimation""" 
# files = {}
# filenames = glob.glob("*.txt")
n =  3 # skipheader rows
# col1, col2 = 1, 2   # Used columns
# for i in range(len(filenames)):
#     re, im, hz = np.loadtxt(filenames[i], skiprows=n, unpack=True)

#     files[i] = hz, np.sqrt(im**2+re**2)

# for i in range(len(filenames)):
#     plt.plot(files[i][0], files+[i][1]*1e-6, label=filenames[i])
    
# plt.legend()


# Fit to lorentian and gaussian lineshape for the resonance frequency plot
# Formulae from Wikipedia
def gaussian(x, mu, sigma):
    """Gauss-function, with amplitude dependent on sigma and mu.
    """
    return (sigma*np.sqrt(2*np.pi))**(-1)*np.exp(-(x-mu)**2/(2*sigma**2))


def gaussian2(x, A, mu, sigma):
    """Gauss-function with arbitrary amplitude A.
    Parameters
    ----------
    x : array.
    A : integer, amplitude.
    mu : integer, expectation value.
    sigma : integer, standard deviation.
    """
    return A*np.exp(-(x-mu)**2/(2*sigma**2))
    

def lorentzian(x, x0, gamma, A):
    """
    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    x0 : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    """
    return A*(gamma/2)/(np.pi*((x-x0)**2 + (gamma/2)**2))

"""Part 2"""
    
def t1fit(x, A, T1):
        """
        Fit-function for the spin-lattice-relaxation (Mz dependent)
        to the formula:      M_{Echo}(dt) \propto M_{Sat}*(1 - exp(-dt/T1))
        
        Parameters
        ----------
        x : array, tau in seconds.
        A : integer, M_{sat}
        T1 : integer, spin-lattice-relaxation time in seconds
        """
        return 
    

def t2fit(x, A, T2):
        """
        Fit-function for the spin-spin-relaxation (Mx, My dependent)
        to the formula:      M_{Echo}(tau) \propto M_{Sat}*exp(-2*tau/T2)
        
        Parameters
        ----------
        x : array, tau in seconds.
        A : integer, M_{sat}
        T2 : integer, spin-spin-relaxation time in seconds
        """
        return 

"""Part 2"""

tau = np.logspace(0, 2.3, 15)*5 # in micro seconds used for T2 spin-spin-relaxation

def main():
    """Hauptprogramm"""
    print(__doc__)
    
    """Part 1"""
    # Load data
    re45, im45, f45 = np.loadtxt("45.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    re452, im452, f452 = np.loadtxt("452.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    re454, im454, f454 = np.loadtxt("454.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    re456, im456, f456 = np.loadtxt("456.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    re458, im458, f458 = np.loadtxt("458.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    re46, im46, f46 = np.loadtxt("46.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    reres, imres, fres = np.loadtxt("resonance.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    
    # Plot all the different frequency measurements
    # niceplot(x=f45, y=np.sqrt(re45**2+im45**2), c="tab:blue", plotlabel=r'$\nu$ = 45MHz', legend=True,
    #          plot2=True, plot3=True, plot4=True, plot5=True, plot6=True, 
    #          plot7=True, xaxis=r'$\nu$ / MHz', yaxis=r'Intensity / $\times 10^6$',
    #          titel=r'Finding the resonance frequency $\nu_{res}$',
    #          x2=f452, y2=np.sqrt(re452**2+im452**2), c2="tab:orange", 
    #          safefig=True, safename='frequencies', plotlabel2=r'$\nu$ = 45.2MHz',
    #          x4=f454, y4=np.sqrt(re454**2+im454**2), c4="tab:pink", 
    #          plotlabel4=r'$\nu$ = 45.4MHz',
    #          x5=f456, y5=np.sqrt(re456**2+im456**2), c5="tab:red", 
    #          plotlabel5=r'$\nu$ = 45.6MHz',
    #          x6=f458, y6=np.sqrt(re458**2+im458**2), c6="tab:purple", 
    #          plotlabel6=r'$\nu$ = 45.8MHz',
    #          x7=f46, y7=np.sqrt(re46**2+im46**2), c7="tab:brown", 
    #          plotlabel7=r'$\nu$ = 46MHz',
    #          x3=fres, y3=np.sqrt(reres**2+imres**2), c3="tab:green", 
    #          plotlabel3=r'$\nu_{res}$ = 45.385MHz')

    # Generate parameters and errors for both pulse signal 
    # popt1, cov1 = fit(gaussian, fres, np.sqrt(reres**2+imres**2), sigma=None,
                      # absolute_sigma=True)
    popt3, cov3 = fit(gaussian2, fres, np.sqrt(reres**2+imres**2), sigma=None,
                      absolute_sigma=True)
    popt2, cov2 = fit(lorentzian, fres, np.sqrt(reres**2+imres**2), sigma=None,
                      absolute_sigma=True)
    
    
    # Plot of the resonance frequency with the two fits
    niceplot(safefig=True, safename='fitgausslorentz', xaxis=r'$\nu$ / MHz', 
             yaxis=r'Intensity / $\times 10^6$', 
             titel='Comparing different Fit models to the peak', legend=True, 
             plot2=False, plot3=True, plot4=True,
             x=fres, y=np.sqrt(reres**2+imres**2)*1e-6, c='k', 
             plotlabel=r'$\nu_{res}$ = 45.385MHz',
             #x2=fres, y2=gaussian(fres, *popt1), c2='tab:orange', 
             plotlabel2='gaussian',
             x3=np.linspace(-100, 100, 1000), 
             y3=lorentzian(np.linspace(-100, 100, 1000), *popt2)*1e-6, 
             c3='tab:blue', 
             plotlabel3='lorentzian',
             x4=np.linspace(-100, 100, 1000), 
             y4=gaussian2(np.linspace(-100, 100, 1000), *popt3)*1e-6, 
             c4='tab:orange', 
             plotlabel4='gaussian2',
             ls='', marker='s', lw2=0.5, lw3=3, lw4=3, ms=5)
             
    
    # print("Die Parameter des gaussian-Fits sind:")
    # print("mu+-∆mu={}+-{}".format(popt1[0], np.sqrt(np.diag(cov1))[0]))
    # print("sigma+-∆sigma={}+-{}".format(popt1[1], np.sqrt(np.diag(cov1))[1]))
    print("Die Parameter des gaussian2-Fits sind:")
    print("A+-∆A={}+-{}".format(popt3[0], np.sqrt(np.diag(cov3))[0]))
    print("mu+-∆mu={}+-{}".format(popt3[1], np.sqrt(np.diag(cov3))[1]))
    print("sigma+-∆sigma={}+-{}".format(popt3[2], np.sqrt(np.diag(cov3))[2]))
    print("The gaussian-line-with is: tau = {}".format(popt3[2]))
    print("Die Parameter des lorentzian-Fits sind:")
    print("x0+-∆x0={}+-{}".format(popt2[0], np.sqrt(np.diag(cov2))[0]))
    print("gamma+-∆gamma={}+-{}".format(popt2[1], np.sqrt(np.diag(cov2))[1]))
    print("A+-∆A={}+-{}".format(popt2[2], np.sqrt(np.diag(cov2))[2]))
    print("The lorentzian-line-with is: tau = {}".format(popt2[1]))
    
    
    """Part 2"""
    # popt4, cov4 = t1fit()
    
    
    plt.show()
    
if __name__ == "__main__":
    main()