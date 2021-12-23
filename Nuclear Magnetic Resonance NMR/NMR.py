#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Thu Dec 16 11:20:08 2021

#@author: stevengebel

"""F-Praktikum Auswertung NMR"""


import numpy as np
import matplotlib.pyplot as plt
import glob 
from scipy import integrate
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
             plotlabel5=None, plotlabel6=None, plotlabel7=None, titel=None,
             xlim=None, ylim=None
             ):
    
    fig = plt.figure(figsize=size) 
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel(r'{}'.format(yaxis), size=22)
    ax.set_xlabel(r'{}'.format(xaxis), size=22)
    ax.set_title(r'{}'.format(titel),size=25)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    #####
    # pinmp_ticks(ax.xaxis, ticks) # For special plotstyle
    # pinmp_ticks(ax.yaxis, ticks)
    # ax.grid(which='minor', alpha=.5)
    # ax.grid(which='major', alpha=.5)
    # ax.tick_params(which='both')
    #####
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if error:
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
    if  plot2:
        ax.plot(x2, y2, ms=ms, lw=lw2, ls=ls2, marker=marker2, c=c2, 
             label=r'{}'.format(plotlabel2))
    if  plot3:
        ax.plot(x3, y3, ms=ms, lw=lw3, ls=ls3, marker=marker3, c=c3, 
             label=r'{}'.format(plotlabel3))
    if  plot4:
        ax.plot(x4, y4, ms=ms, lw=lw4, ls=ls4, marker=marker4, c=c4, 
             label=r'{}'.format(plotlabel4))
    if  plot5:
        ax.plot(x5, y5, ms=ms, lw=lw5, ls=ls5, marker=marker5, c=c5, 
             label=r'{}'.format(plotlabel5)) 
    if  plot6:
        ax.plot(x6, y6, ms=ms, lw=lw6, ls=ls6, marker=marker6, c=c6,
                label=r'{}'.format(plotlabel6))   
    if  plot7:
        ax.plot(x7, y7, ms=ms, lw=lw7, ls=ls7, marker=marker7, c=c7, 
                label=r'{}'.format(plotlabel7))
    if legend:
        ax.legend(fontsize=fs, markerscale=ms-8, facecolor='white')
    if safefig:
        plt.savefig('{}.svg'.format(safename), dpi=300)

    plt.show()


"""Part 1: Centre frequency and parameter estimation""" 
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
def t1fit(x, A, T1, C):
        """
        Fit-function for the spin-lattice-relaxation (Mz dependent)
        to the formula:      M_{Echo}(dt) \propto M_{Sat}*(1 - exp(-dt/T1))
        
        Parameters
        ----------
        x : array, tau in seconds.
        A : integer, M_{sat}
        T1 : integer, spin-lattice-relaxation time in seconds
        """
        return A*(1 - np.exp(-x/T1)) + C 
    

def t2fit(x, A, T2, B):
        """
        Fit-function for the spin-spin-relaxation (Mx, My dependent)
        to the formula:      M_{Echo}(tau) \propto M_{Sat}*exp(-2*tau/T2)
        
        Parameters
        ----------
        x : array, tau in seconds.
        A : integer, M_{sat}
        T2 : integer, spin-spin-relaxation time in seconds
        """
        return A*np.exp(-2*x/T2) + B
    # Bemerkung: wenn die 2 aus dieser fitunktion entfernt wird, spuckt der fit für T2
    # den Wert aus, den auch NTNMR ausgegeben hat ~750. Im Skript steht die formel mit 
    # -2*tau, s.d. wir hier ~1500 als fitparameter bekommen...
    # Vergleich mit Wikipedia: https://en.wikipedia.org/wiki/Relaxation_%28NMR%29#T1_and_T2
    # -> T2: da ist es auch ohne die -2...
    
    
def intensityintegral(re, im, t, n, m):
        """
        Compute integrals over signal magnitude for each measurement
        
        Parameters
        ----------
        re : array_like, real part.
        im : array_like, imaginary part.
        t : array_like, abscissa data.
        n : integer, interval length.
        m : integer, amount of measurements.

        Returns
        -------
        x :  absciss data, sliced size = n x m.
        magn : magnitude of the complex signal, size = n x m
        integmagn : array_like, integrated signal, size = m
        """
        x = []
        magn = []
        integmagn = [] # integrated signal via simpson method
        maximum = [] # maximum magnitude of each measurement 
        for i in range(m):
            integmagn.append(integrate.simpson(y=np.sqrt(re[i*n:i*n+n]**2 + 
                                            im[i*n:i*n+n]**2),
                                            x=t[i*n:i*n+n]))
            magn.append(np.sqrt(re[i*n:i*n+n]**2 + im[i*n:i*n+n]**2))
            x.append(t[i*n:i*n+n])
            maximum.append(max(np.sqrt(re[i*n:i*n+n]**2 + im[i*n:i*n+n]**2)))
        
        # print("maximum = {}".format(maximum))
        # print("integmagn = {}".format(integmagn))
        return x, magn, integmagn


"""Part 2"""
tau = np.logspace(0, 2.3, 15)*5 # in micro seconds used for T2 spin-spin-relaxation

def main():
    """Hauptprogramm"""
    print(__doc__)
    n=3
    """Part 1"""
    # # Load data
    # re45, im45, f45 = np.loadtxt("45.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    # re452, im452, f452 = np.loadtxt("452.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    # re454, im454, f454 = np.loadtxt("454.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    # re456, im456, f456 = np.loadtxt("456.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    # re458, im458, f458 = np.loadtxt("458.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    # re46, im46, f46 = np.loadtxt("46.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    # reres, imres, fres = np.loadtxt("resonance.txt", unpack=True, skiprows=n, usecols=(0,1,2))
    
    # # Plot all the different frequency measurements
    # niceplot(x=f45, y=np.sqrt(re45**2+im45**2)*1e-6, c="tab:blue", 
    #           plotlabel=r'$\nu$ = $\,$45MHz', legend=True,
    #           plot2=True, plot3=True, plot4=True, plot5=True, plot6=True, 
    #           plot7=True, xaxis=r'$\nu$ / MHz', yaxis=r'Intensity / $\times 10^6$',
    #           titel=r'Finding the resonance frequency $\nu_{res}$',
    #           x2=f452, y2=np.sqrt(re452**2+im452**2)*1e-6, c2="tab:orange", 
    #           safefig=True, safename='frequencies', plotlabel2=r'$\nu$ = $\,$45.2MHz',
    #           x4=f454, y4=np.sqrt(re454**2+im454**2)*1e-6, c4="tab:pink", 
    #           plotlabel4=r'$\nu$ = $\,$45.4MHz',
    #           x5=f456, y5=np.sqrt(re456**2+im456**2)*1e-6, c5="tab:red", 
    #           plotlabel5=r'$\nu$ = $\,$45.6MHz',
    #           x6=f458, y6=np.sqrt(re458**2+im458**2)*1e-6, c6="tab:purple", 
    #           plotlabel6=r'$\nu$ = $\,$45.8MHz',
    #           x7=f46, y7=np.sqrt(re46**2+im46**2)*1e-6, c7="tab:brown", 
    #           plotlabel7=r'$\nu$ = $\,$46MHz',
    #           x3=fres, y3=np.sqrt(reres**2+imres**2)*1e-6, c3="tab:green", 
    #           plotlabel3=r'$\nu_{res}$ = 45.385MHz', xlim=(-650, 650),
    #           lw=3, lw2=3, lw3=3, lw4=3, lw5=3, lw6=3, lw7=3, fs=16,
    #           ylim=(0,2.8)
    # )

    # # Generate parameters and errors for both pulse signal 
    # # popt1, cov1 = fit(gaussian, fres, np.sqrt(reres**2+imres**2), sigma=None,
    #                   # absolute_sigma=True)
    # popt3, cov3 = fit(gaussian2, fres, np.sqrt(reres**2+imres**2), sigma=None,
    #                   absolute_sigma=True)
    # popt2, cov2 = fit(lorentzian, fres, np.sqrt(reres**2+imres**2), sigma=None,
    #                   absolute_sigma=True)
    
    
    # # Plot of the resonance frequency with the two fits:
    # niceplot(safefig=True, safename='fitgausslorentz', xaxis=r'$\nu$ / MHz', 
    #           yaxis=r'Intensity / $1\times 10^6$', 
    #           titel='Comparing different Fit models to the peak', legend=True, 
    #           plot2=False, plot3=True, plot4=True,
    #           x=fres, y=np.sqrt(reres**2+imres**2)*1e-6, c='k', 
    #           plotlabel=r'$\nu_{res}$ = $\,$45.385MHz',
    #           #x2=fres, y2=gaussian(fres, *popt1), c2='tab:orange', 
    #           plotlabel2='gaussian',
    #           x3=np.linspace(-100, 100, 1000), 
    #           y3=lorentzian(np.linspace(-100, 100, 1000), *popt2)*1e-6, 
    #           c3='tab:blue', 
    #           plotlabel3='Lorentzian',
    #           x4=np.linspace(-100, 100, 1000), 
    #           y4=gaussian2(np.linspace(-100, 100, 1000), *popt3)*1e-6, 
    #           c4='tab:orange', 
    #           plotlabel4='Gaussian',
    #           ls='', marker='s', lw2=0.5, lw3=3, lw4=3, ms=5, ls3='-.', ls4='--',
    #           xlim=(-100,100), ylim=(0,2.8)
    # )
    # # print("Die Parameter des gaussian-Fits sind:")
    # # print("mu+-∆mu={}+-{}".format(popt1[0], np.sqrt(np.diag(cov1))[0]))
    # # print("sigma+-∆sigma={}+-{}".format(popt1[1], np.sqrt(np.diag(cov1))[1]))
    # print("Die Parameter des gaussian2-Fits sind:")
    # print("A+-∆A={}+-{}".format(popt3[0], np.sqrt(np.diag(cov3))[0]))
    # print("mu+-∆mu={}+-{}".format(popt3[1], np.sqrt(np.diag(cov3))[1]))
    # print("sigma+-∆sigma={}+-{}".format(popt3[2], np.sqrt(np.diag(cov3))[2]))
    # print("The gaussian-line-with is: tau = {}".format(popt3[2]))
    # print("Die Parameter des lorentzian-Fits sind:")
    # print("x0+-∆x0={}+-{}".format(popt2[0], np.sqrt(np.diag(cov2))[0]))
    # print("gamma+-∆gamma={}+-{}".format(popt2[1], np.sqrt(np.diag(cov2))[1]))
    # print("A+-∆A={}+-{}".format(popt2[2], np.sqrt(np.diag(cov2))[2]))
    # print("The lorentzian-line-with is: tau = {}".format(popt2[1]))
    
    
    """Part 2"""
    """T1-measurement"""
    # "T1_raw.txt" ->range(100) (brauchbar)
    # "T1_manual.txt" -> range(10) (unbrauchbar)
    # "T1_manualext.txt" -> range(14) (brauchbar)
    # "T1_auto30.txt" -> range(30) (UNBRAUCHBAR, wird gar nicht benutzt)
    # Bemerkung: habe bei dne dateien 1-2 zeilen txt datei gelöscht, weil python sonst fehler
    # etc. code anzeigte 
    
    # Import data:
    re12, im12, t12 = np.loadtxt("T1_manualext.txt", usecols=(0,1,2), unpack=True, skiprows=4)
    dt12, k12 = np.loadtxt("T1_manualext_fit.txt", usecols=(0,1), unpack=True, skiprows=3)
    
    re14, im14, t14 = np.loadtxt("T1_raw.txt", usecols=(0,1,2), unpack=True, skiprows=4)
    dt14 = np.linspace(0,20,100) 
    
    # Compute integrals: 
    x2, magn2, integmagn2 = intensityintegral(re12, im12, t12, n=1024, m=14)  
    x4, magn4, integmagn4 = intensityintegral(re14, im14, t14, n=1024, m=100)
    
    # Obtain fit-parameters for each 3 measurement sets:
    popt7, cov7 = fit(t1fit, dt12, k12, p0=[2e7, 0.9, 2e7], 
                      sigma=None, absolute_sigma=True)
    popt9, cov9 = fit(t1fit, dt14[0:50], integmagn4[0:50], p0=[1e7, 0.9, 1],
                        sigma=None, absolute_sigma=True)
    
    # Fit parameters for spin-lattice relaxation time T1:
    print("MANUALEXT-T1-FIT PARAMETERS:")
    print("A+-∆A={}+-{}".format(popt7[0], np.sqrt(np.diag(cov7))[0])) 
    print("T1+-∆T1={}+-{} / ms??".format(popt7[1], np.sqrt(np.diag(cov7))[1]))
    print("C+-∆C={}+-{} ".format(popt7[2], np.sqrt(np.diag(cov7))[2]))
    print("exponent={} ".format(-dt12[3]/popt7[1]))
    print("RAW-T1-FIT PARAMETERS:")
    print("A+-∆A={}+-{}".format(popt9[0], np.sqrt(np.diag(cov9))[0]))
    print("T1+-∆T1={}+-{} / ms??".format(popt9[1], np.sqrt(np.diag(cov9))[1]))
    print("C+-∆C={}+-{} ".format(popt9[2], np.sqrt(np.diag(cov9))[2]))
    print("exponent={} ".format(-dt14[3]/popt9[1]))

    # Plot fits from 3 measurement sets:
    # niceplot(x=dt14[0:50], y=np.asarray(integmagn4)[0:50], c='k',
    #           x2=np.linspace(0,40,1000), y2=t1fit(np.linspace(0,40,1000), *popt9), 
    #           x3=dt12, y3=np.asarray(integmagn2), c3='g', c4='g',
    #           x4=np.linspace(0,40,1000), y4=t1fit(np.linspace(0,40,1000), *popt7), 
    #           c2='k', plotlabel=r'raw', plotlabel2=r'fit raw', 
    #           plotlabel3='manualext', plotlabel4='fit manualext',
    #           ls='',  marker='s', plot2=True, plot3=True, plot4=True, 
    #           lw2=3, lw4=3, ls3='', marker3='s', 
    #           xaxis=r'$\Delta$t / ms', yaxis=r'Intensity / $1\times 10^7$', 
    #           titel=r'Determining the spin-lattice relaxation time  T$_2$', legend=True,
    #           safefig=True, safename='T1_plotfits',
    #           xlim=(-0.5,11), ylim=(0.9*min(np.asarray(integmagn4)),
    #                                 max(np.asarray(integmagn2))+0.5e7), lw7=5,
    #           plot7=True, x7=np.linspace(0,40,1000), y7=814933 + 403261*(1 - np.exp(-np.linspace(0,40,1000)/0.28384812))
    #           )
    niceplot(x=dt14[0:50], y=np.asarray(integmagn4)[0:50], c='k', c2='k',
             x2=np.linspace(0,40,1000), y2=t1fit(np.linspace(0,40,1000), *popt9), 
             plotlabel='raw', plotlabel2='fit raw', legend=True,
             xlim=(-0.5,11), ylim=(0.9*min(np.asarray(integmagn4)),
                                    max(np.asarray(integmagn2))+0.5e7),
             ls='', marker='s', plot2=True, lw2=3)
    
    niceplot(x=dt12, y=np.asarray(integmagn2), c='k', c2='k',
             x2=np.linspace(0,10,1000), y2=t1fit(np.linspace(0,10,1000), *popt7), 
             plotlabel=r'manualext', plotlabel2=r'fit manualext', 
             xlim=(-0.5,11), ylim=(0.9*min(np.asarray(integmagn4)),
                                    max(np.asarray(integmagn2))+0.5e7),
             ls='', marker='s', plot2=True, lw2=3, 
             plot3=True, x3=np.linspace(0, 20, 100), y3=4e7+64403261*(1-np.exp(-np.linspace(0, 20, 100)/0.2))
             )
    
    
    """T2-measurement"""
    # # # NTNMR fit: T2= 746.389734 ± 28.430624 mikroseconds
    # tau, fit21, fit22 = np.loadtxt("T2fit.txt", unpack=True, 
    #                                     usecols=(0,1,2), skiprows=0)
    # re2, im2, t2 = np.loadtxt("T2_rawdata.txt", usecols=(0,1,2), 
    #                           skiprows=3, unpack=True)
    
    # # Compute integrals over signal magnitude for each measurement (12):
    # x5, magn5, integmagn5 = intensityintegral(re2, im2, t2, 512, 12)
    
    # # Plot for various tau (5 chosen spread values) (not mandatory):
    # niceplot(x=x5[0], y=magn5[0]*1e-6, c5='tab:blue',
    #           x2=x5[2], y2=magn5[2]*1e-6, c4='tab:orange',
    #           x3=x5[5], y3=magn5[5]*1e-6, c3='tab:green',
    #           x4=x5[8], y4=magn5[8]*1e-6, c2='tab:red',
    #           x5=x5[11], y5=magn5[11]*1e-6, c='tab:purple',
    #           plotlabel=r'$\tau=11\,\mu$s', plotlabel2=r'$\tau=48\,\mu$s',
    #           plotlabel3=r'$\tau=151\,\mu$s', plotlabel4=r'$\tau=468\,\mu$s',
    #           plotlabel5=r'$\tau=1000\,\mu$s',
    #           plot2=True, plot3=True, plot4=True, plot5=True,
    #           lw=3, lw2=3, lw3=3, lw4=3, lw5=3,
    #           xaxis=r'', yaxis=r'Intensity / $1\times 10^6$', 
    #           titel=r'Signals in dependence of $\tau$', legend=True,
    #           safefig=True, safename='T1plot',
    #           xlim=(0,100), ylim=(0,3.85)
    # )
    
    # # Fit parameters for spin-spin relaxation time T2 from the 
    # # integrated signals:
    # popt5, cov5 = fit(t2fit, tau, integmagn5, p0=[0.7e8, 746.389734], 
    #                   bounds=[[0.5e8, 500],[1e8, 1800]], sigma=None, 
    #                   absolute_sigma=True)
    # # Bemerkung: Fit scheint für mittlere tau und große tau rel. schlecht zu passen
    # print("The T2-fit parameters are:")
    # print("A+-∆A={}+-{}".format(popt5[0], np.sqrt(np.diag(cov5))[0]))
    # print("T2+-∆T2={}+-{} / microseconds".format(popt5[1], np.sqrt(np.diag(cov5))[1]))
    # print("B+-∆B={}+-{} ".format(popt5[2], np.sqrt(np.diag(cov5))[2]))
    # niceplot(x=np.linspace(0,1000, 1000), 
    #           y=t2fit(np.linspace(0,1000, 1000), *popt5)*1e-7, lw=3,
    #           c='k', xaxis=r'$\tau$ / $\mu$s', 
    #           yaxis=r'Intensity / $1\times 10^7$',
    #           x2=tau, y2=np.asarray(integmagn5)*1e-7, c2='k', ls2='', 
    #           marker2='s', plotlabel='fit', plotlabel2='experimental data', 
    #           legend=True, safefig=True, safename='T2fit', 
    #           titel=r'Determining the Spin-Spin Relaxation time T$_2$', 
    #           ylim=(2, (max(t2fit(tau, *popt5)) + 0.5)*1e-7), 
    #           xlim=(0,1000), plot2=True, ls='-.'
    # )

    plt.show()
    
    
if __name__ == "__main__":
    main()
    

    
    
    
    