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
             xlim=None, ylim=None):
    
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
        ax.legend(fontsize=fs, markerscale=ms-8, facecolor='white')
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
        return A*(1 - np.exp(-x/T1))
    

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
        return A*np.exp(-2*x/T2)
    # Bemerkung: wenn die 2 aus dieser fitunktion entfernt wird, spuckt der fit für T2
    # den Wert aus, den auch NTNMR ausgegeben hat ~750. Im Skript steht die formel mit 
    # -2*tau, s.d. wir hier ~1500 als fitparameter bekommen...
    # Vergleich mit Wikipedia: https://en.wikipedia.org/wiki/Relaxation_%28NMR%29#T1_and_T2
    # -> T2: da ist es auch ohne die -2...


"""Part 2"""
tau = np.logspace(0, 2.3, 15)*5 # in micro seconds used for T2 spin-spin-relaxation

def main():
    """Hauptprogramm"""
    print(__doc__)
    
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
    #           ls='', marker='s', lw2=0.5, lw3=3, lw4=3, ms=5, 
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
    
    
    # """Part 2"""
    # """T1-measurement"""
    # Text-datei zu ersetzen mit: "T1_raw.txt" ->range(100), 
    
    # NTNMR-Fit: T1 = 6792560298301635 ± 1841978447959600700000000000000
    re1, im1, t1 = np.loadtxt("T1_manualext.txt", usecols=(0,1,2), unpack=True, skiprows=3)
    dt = np.loadtxt("T1_manualext_fit.txt", usecols=0, unpack=True, skiprows=3)
    
    # Compute integrals over signal magnitude for each measurement (12):
    integmagntrapz1 = [] # integral via simpson and trapezoid method 
    integmagnsimps1 = []
    t1plotx = []
    t1ploty = []
    for i in range(14):
        integmagntrapz1.append(np.trapz(y=np.sqrt(re1[i*1024:i*1024+1024]**2 + 
                                            im1[i*1024:i*1024+1024]**2),
                                            x=t1[i*1024:i*1024+1024]))
        integmagnsimps1.append(integrate.simpson(y=np.sqrt(re1[i*1024:i*1024+1024]**2 + 
                                            im1[i*1024:i*1024+1024]**2),
                                            x=t1[i*1024:i*1024+1024]))
        t1ploty.append(np.sqrt(re1[i*1024:i*1024+1024]**2 + 
                                im1[i*1024:i*1024+1024]**2))
        t1plotx.append(t1[i*1024:i*1024+1024])
    
    print("intensity1 = {}".format(integmagntrapz1))  
    print("dt = {}".format(dt))
    
        
        
    # # Plot for various tau (5 chosen spread values) (not mandatory):
    # niceplot(x=t1plotx[0], y=t1ploty[0]*1e-6, c='tab:blue',
    #          x2=t1plotx[-1], y2=t1ploty[-1]*1e-6, c2='tab:orange',
    # #           x3=t2plotx[5], y3=t2ploty[5]*1e-6, c3='tab:green',
    # #           x4=t2plotx[8], y4=t2ploty[8]*1e-6, c2='tab:red',
    # #           x5=t2plotx[11], y5=t2ploty[11]*1e-6, c='tab:purple',
    #           plotlabel=r'$\Delta$t=11\,\mu$s', plotlabel2=r'$\Delta$t=48\,\mu$s',
    # #           plotlabel3=r'$\tau=151\,\mu$s', plotlabel4=r'$\tau=468\,\mu$s',
    # #           plotlabel5=r'$\tau=1000\,\mu$s',
    #           plot2=True,
    #           # plot3=True, plot4=True, plot5=True,
    # #           lw=3, lw2=3, lw3=3, lw4=3, lw5=3,
    # #           xaxis=r'$\Delta$t / $\mu$s', yaxis=r'Intensity / $1\times 10^6$', 
    #           titel=r'Signals in dependence of $\tau$', legend=True,
    # #           safefig=True, safename='T1plot',
    #           xlim=(0,200), ylim=(0, max(t1ploty[0])*1.05*1e-6)
    # )
    
    # Fit parameters for spin-lattice relaxation time T1:
    popt4, cov4 = fit(t1fit, dt, integmagntrapz1, 
                      sigma=None, absolute_sigma=True)
    print("The T1-fit parameters are:")
    print("A+-∆A={}+-{}".format(popt4[0], np.sqrt(np.diag(cov4))[0]))
    print("T1+-∆T1={}+-{} / seconds".format(popt4[1], np.sqrt(np.diag(cov4))[1]))
    
    niceplot(x=np.linspace(0,200, 500), ls='', marker='.',
              y=t1fit(np.linspace(0,200, 500), *popt4)*1e-9, 
              c='tab:blue', xaxis=r'$\Delta$t / $\mu$s', yaxis=r'Intensity / $1\times 10^6$',
    #          x2=fit2x, y2=np.asarray(integmagnsimps)*1e-7, c2='k', ls2='', 
    #          marker2='s', plotlabel='fit', plotlabel2='experimental data', 
    #          legend=True, safefig=True, safename='T2fit',  plot2=True
    #          titel=r'Determining the Spin-Spin Relaxation time T$_2$', 
              ylim=(0, (max(t1fit(dt, *popt4)) + 0.5)*1e-9), xlim=(0,10)
    )
    
    
    # """T2-measurement"""
    # # NTNMR fit: T2= 746.389734 ± 28.430624 mikroseconds
    # # fit21 entspricht maximum des realteils jeder messung, fit22 = ??
    # # fit2x sind tau-Werte
    # fit2x, fit21, fit22 = np.loadtxt("T2fit.txt", unpack=True, 
    #                                     usecols=(0,1,2), skiprows=0)
    # re2, im2, t2 = np.loadtxt("T2_rawdata.txt", usecols=(0,1,2), 
    #                           skiprows=3, unpack=True)
    
    # # Compute integrals over signal magnitude for each measurement (12):
    # integmagntrapz2 = [] # integral via simpson and trapezoid method 
    # integmagnsimps2 = []
    # t2plotx = []
    # t2ploty = []
    # for i in range(12):
    #     integmagntrapz2.append(np.trapz(y=np.sqrt(re2[i*512:i*512+512]**2 + 
    #                                         im2[i*512:i*512+512]**2),
    #                                         x=t2[i*512:i*512+512]))
    #     integmagnsimps2.append(integrate.simpson(y=np.sqrt(re2[i*512:i*512+512]**2 + 
    #                                         im2[i*512:i*512+512]**2),
    #                                         x=t2[i*512:i*512+512], dx=2))
    #     t2ploty.append(np.sqrt(re2[i*512:i*512+512]**2 + 
    #                            im2[i*512:i*512+512]**2))
    #     t2plotx.append(t2[i*512:i*512+512])
    # print("intensity2 = {}".format(integmagntrapz2))   
    # print("tau = {}".format(fit2x)) # used values for tau
    
    # # Plot for various tau (5 chosen spread values) (not mandatory):
    # niceplot(x=t2plotx[0], y=t2ploty[0]*1e-6, c5='tab:blue',
    #           x2=t2plotx[2], y2=t2ploty[2]*1e-6, c4='tab:orange',
    #           x3=t2plotx[5], y3=t2ploty[5]*1e-6, c3='tab:green',
    #           x4=t2plotx[8], y4=t2ploty[8]*1e-6, c2='tab:red',
    #           x5=t2plotx[11], y5=t2ploty[11]*1e-6, c='tab:purple',
    #           plotlabel=r'$\tau=11\,\mu$s', plotlabel2=r'$\tau=48\,\mu$s',
    #           plotlabel3=r'$\tau=151\,\mu$s', plotlabel4=r'$\tau=468\,\mu$s',
    #           plotlabel5=r'$\tau=1000\,\mu$s',
    #           plot2=True, plot3=True, plot4=True, plot5=True,
    #           lw=3, lw2=3, lw3=3, lw4=3, lw5=3,
    #           xaxis=r'$\tau$ / $\mu$s', yaxis=r'Intensity / $1\times 10^6$', 
    #           titel=r'Signals in dependence of $\tau$', legend=True,
    #           safefig=True, safename='T1plot',
    #           xlim=(0,100), ylim=(0,3.85)
    # )
    
    # # Fit parameters for spin-spin relaxation time T2 from the 
    # # integrated signals:
    # popt5, cov5 = fit(t2fit, fit2x, integmagnsimps, p0=[0.7e8, 746.389734], 
    #                   bounds=[[0.5e8, 500],[1e8, 1800]], sigma=None, 
    #                   absolute_sigma=True)
    # # Bemerkung: Fit scheint für mittlere tau und große tau rel. schlecht zu passen
    # print("The T2-fit parameters are:")
    # print("A+-∆A={}+-{}".format(popt5[0], np.sqrt(np.diag(cov5))[0]))
    # print("T2+-∆T2={}+-{} / microseconds".format(popt5[1], np.sqrt(np.diag(cov5))[1]))
    
    # niceplot(x=np.linspace(0,1000, 1000), 
    #          y=t2fit(np.linspace(0,1000, 1000), *popt5)*1e-7, lw=3,
    #          c='tab:blue', xaxis=r'$\tau$ / $\mu$s', yaxis=r'Intensity / $1\times 10^7$',
    #          x2=fit2x, y2=np.asarray(integmagnsimps)*1e-7, c2='k', ls2='', 
    #          marker2='s', plotlabel='fit', plotlabel2='experimental data', 
    #          legend=True, safefig=True, safename='T2fit', 
    #          titel=r'Determining the Spin-Spin Relaxation time T$_2$', 
    #          ylim=(2, (max(t2fit(fit2x, *popt5)) + 0.5)*1e-7), 
    #          xlim=(0,1000), plot2=True
    # )

    plt.show()
    
    
if __name__ == "__main__":
    main()
    

    
    
    
    