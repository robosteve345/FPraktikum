#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Wed Nov 24 13:23:06 2021

#@author: stevengebel

"""F-Praktikum Auswertung: Solar cell
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as fit
import matplotlib.ticker as ticker
plt.style.use('ggplot') # Plotstyle 
plt.rcParams.update({
    'font.family': 'serif',
    'text.usetex': False,
    'pgf.rcfonts': False,
})


def cosfit(x, A, B):
    return A*np.cos(B*x)


def linearfit(x, A, B):
    return A*x + B


# Funktion die versch. Parameter findet
def MPP(U, I, a):
    """
    Parameters
    ----------
    U : array
        Spannung.
    I : array
        Strom.
    a : array
        Anzahl der solarpanels für aufgenommene Messungen, die in den I,U-arrays sind.

    Returns
    -------
    I0 : Kurzschlussstrom
    V_OC : Leerlaufspannung
    PMLP = U_MPP*I_MPP
    FF = U_MPP*I_MPP / I0*V_OC
    eta 
    """
    # Konstanten
    A = 26 # area of an inorganic cell
    X = 1/3 # Anteil an i_0, i = X*i_0
    i_0 = 100e-3 # in W/cm^2, "Sonne!"
    # Parameter
    I0 = [] # in A 
    V_OC = [] # in V
    index_PMLP = []
    PMLP = [] # in mW
    FF = []
    eta = []
    for i in range(len(U)):
        P_in = a[i]*A*X*i_0 # in W, Leistung = Intensität*Fläche
        I0.append(I[i][np.argmin(abs(U[i]))])  
        V_OC.append(U[i][np.argmin(abs(I[i]))]) 
        index_PMLP.append(np.argmin(U[i]*I[i]))
        PMLP.append((U[i][np.argmin(U[i]*I[i])]*I[i][np.argmin(U[i]*I[i])])*1e3) 
        FF.append((((I[i][np.argmin(abs(U[i]))]*U[i][np.argmin(abs(I[i]))]) / 
                 (U[i][np.argmin(U[i]*I[i])]*I[i][np.argmin(U[i]*I[i])]))**(-1)))
        eta.append((((I[i][np.argmin(abs(U[i]))]*U[i][np.argmin(abs(I[i]))]) / 
                 (U[i][np.argmin(U[i]*I[i])]*I[i][np.argmin(U[i]*I[i])]))**(-1))*np.abs(I[i][np.argmin(abs(U[i]))])*U[i][np.argmin(abs(I[i]))]/P_in)
    return I0, V_OC, index_PMLP, PMLP, FF, eta


def pinmp_ticks(axis, ticks):
    """For Nice Plots"""
    axis.set_major_locator(ticker.MaxNLocator(ticks))
    axis.set_minor_locator(ticker.MaxNLocator(ticks * 10))
    return axis


def niceplot(x, y, c, lw=None, lw2=None, lw3=None, lw4=None, lw5=None, ls=None, 
             ls2=None, ls3=None, ls4=None, ls5=None, plot2=None, plot3=None, 
             plot4=None, plot5=None, x2=None, c2=None, y2=None, x3=None, 
             y3=None, c3=None, x4=None, y4=None, x5=None, y5=None, c4=None, 
             c5=None, marker=None, marker2=None, marker3=None , marker4=None,
             marker5=None, ylabel=None, xlabel=None, ms=10, cs=5, fs=20, ticks=6,
             size=(8,8), safefig=None, error=None, errorlabel=None, 
             safename=None, yaxis=None, xaxis=None, yerr=None, plotlabel=None, 
             legend=None, plotlabel2=None, plotlabel3=None, plotlabel4=None,
             plotlabel5=None, titel=None):
    """
    -Very ugly, but does its job-
    
    Parameters
    ----------
    Self-explanatory, plot(2,3,4)=True, safefig=True,
    if needed as well as corresponding 
    other parameters lw, ls etc.

    Returns
    -------
    Plot of the datasets (x_i, y_i), i in {1,2,3,4,5}
    """
    
    fig = plt.figure(figsize=size) 
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel(r'{}'.format(yaxis), size=22)
    ax.set_xlabel(r'{}'.format(xaxis), size=22)
    ax.set_title(r'{}'.format(titel),size=25)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    # ax.set_yticklabels(x, fontsize=fs)
    # ax.set_xticklabels(y, fontsize=fs)
    pinmp_ticks(ax.xaxis, ticks)
    pinmp_ticks(ax.yaxis, ticks)

    ax.grid(which='minor', alpha=.5)
    ax.grid(which='major', alpha=.5)

    ax.tick_params(which='both')
    if error == True:
        ax.errorbar(x, y, xerr=yerr, capsize=cs, c=c, ms=ms, ls=ls,
                     marker=marker, label = r'{}'.format(plotlabel))
    else:
        ax.plot(x, y, ms=ms, lw=lw, ls=ls, marker=marker, c=c, 
                 label=r'{}'.format(plotlabel))
    # ax.fill_between(x[0:4], y[0:4], color='red', alpha=.2)
    # ax.fill_between(x[3:8], y[3:8], color='green', alpha=.2)
    # ax.fill_between(x[7:12], y[7:12], color='blue', alpha=.2, )
    # ax.set_ylim(0)
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
        ax.plot(x5, y5, ms=ms, lw=lw5, ls=ls5, marker=marker5, c=c5)
               # , label=r'{}'.format(plotlabel5))    
    if legend == True:
        ax.legend(fontsize=fs, markerscale=ms-8, facecolor='white')
    if safefig == True:
        plt.savefig('{}.svg'.format(safename), dpi=300)

    plt.show()
    
    
"""C: Solar cell as power deliverer"""
############
# Datensätze
############
# A) 6 panels bright-curve
UCA1, ICA1 = np.loadtxt('Cfullbright.txt', usecols=(0,1), skiprows=0, unpack=True) 

BU, BI = np.loadtxt('B100.txt', usecols=(0,1), skiprows=0, unpack=True) 

# B) 
UCB2, ICB2 = np.loadtxt('C1par68ser.txt', usecols=(0,1), skiprows=0, unpack=True) 
UCB3, ICB3 = np.loadtxt('C68par1ser.txt', usecols=(0,1), skiprows=0, unpack=True ) 
UCB4, ICB4 = np.loadtxt('C33ser33par.txt', usecols=(0,1), skiprows=0, unpack=True ) 


# C)
UCC3, ICC3 = np.loadtxt('Cshademisc.txt', usecols=(0,1), skiprows=0, unpack=True ) 
UCC2, ICC2 = np.loadtxt('Cshadepar.txt', usecols=(0,1), skiprows=0, unpack=True ) 
UCC1, ICC1 = np.loadtxt('Cshadeser.txt', usecols=(0,1), skiprows=0, unpack=True ) 

# D)
UCD1, ICD1 = np.loadtxt('C12fbright.txt', usecols=(0,1), unpack=True)
UCD2, ICD2 = np.loadtxt('C12ebright.txt', usecols=(0,1), unpack=True)
UCD3, ICD3 = np.loadtxt('C12edark.txt', usecols=(0,1), unpack=True)
UCD4, ICD4 = np.loadtxt('C12fdark.txt', usecols=(0,1), unpack=True)

"""D: Influence of temperature"""
UD1, ID1 = np.loadtxt('D30.txt', usecols=(0,1), unpack=True)
UD2, ID2 = np.loadtxt('D35.txt', usecols=(0,1), skiprows=1, unpack=True)
UD3, ID3 = np.loadtxt('D60.txt', usecols=(0,1), skiprows=1, unpack=True)
UD4, ID4 = np.loadtxt('D65.txt', usecols=(0,1), skiprows=1, unpack=True)
dT = np.array([3,1,1,1,1,1,1,2]) # Stat. MU für Temperatur
T = np.array([30, 35.1, 40.1, 45, 50, 55, 60, 65]) # Temperatur in Celsius
V_OC = np.array([581, 561, 555, 547, 543, 535, 532, 510]) # in mV
    
"""E: Diffuse lighting"""
I_OC = np.array([197.5, 196, 194.8, 191.7, 191.9, 186.8, 184.3, 176.3, 
                 150.6, 123.6]) # in mA
angle = np.arange(0, 100, 10)


def main():
    """Hauptprogramm"""
    print(__doc__)
    # Erzeuge Arrays mit Spannungs und Strom-Werten, aller Messungen in C)
    U = np.array([UCA1, UCB2, UCB3, UCB4, UCC1, UCC2, UCC3, UCD1, UCD2, UCD3, UCD4])
    I = np.array([ICA1, ICB2, ICB3, ICB4, ICC1, ICC2, ICC3, ICD1, ICD2, ICD3, ICD4])
    a = np.array([6, 6, 6, 6, 3, 3, 4, 11, 11, 11, 11]) # hier kann auch verschattung mit berücksichtigt werden (weniger zellen beleuchtet)
    # Berechne alle Parameter für Aufgabenteile in C)
    parameter = MPP(U, I, a)
    """C: Solar cell as power deliverer"""
    # print(parameter) #  Hier kann man die Ausgabe bestimmt schöner formulieren, als ich
    # print("Parameter für Aufgabenteil C a): {}A,{}V,{}mW,{},{}".format(parameter[0][0],parameter[1][0],
    #                                                     parameter[3][0],parameter[4][0],
    #                                                     parameter[5][0]))
    # print("Parameter für Aufgabenteil C b1): {}A,{}V,{}mW,{},{}".format(parameter[0][1],parameter[1][1],
    #                                                     parameter[3][1],parameter[4][1],
    #                                                     parameter[5][1]))
    # print("Parameter für Aufgabenteil C b2): {}A,{}V,{}mW,{},{}".format(parameter[0][2],parameter[1][2],
    #                                                     parameter[3][2],parameter[4][2],
    #                                                     parameter[5][2]))
    # print("Parameter für Aufgabenteil C b3): {}A,{}V,{}mW,{},{}".format(parameter[0][3],parameter[1][3],
    #                                                     parameter[3][3],parameter[4][3],
    #                                                     parameter[5][3]))
    
    # print("Parameter für Aufgabenteil C c1): {}A,{}V,{}mW,{},{}".format(parameter[0][4],parameter[1][4],
    #                                                     parameter[3][4],parameter[4][4],
    #                                                     parameter[5][4]))
    # print("Parameter für Aufgabenteil C c2): {}A,{}V,{}mW,{},{}".format(parameter[0][5],parameter[1][5],
    #                                                     parameter[3][5],parameter[4][5],
    #                                                     parameter[5][5]))
    # print("Parameter für Aufgabenteil C c3): {}A,{}V,{}mW,{},{}".format(parameter[0][6],parameter[1][6],
    #                                                     parameter[3][6], parameter[4][6],
    #                                                     parameter[5][6]))
    # print("Parameter für Aufgabenteil C d1): {}A,{}V,{}mW,{},{}".format(parameter[0][7],parameter[1][7],
    #                                                     parameter[3][7],parameter[4][7],
    #                                                     parameter[5][7]))
    # print("Parameter für Aufgabenteil C d2): {}A,{}V,{}mW,{},{}".format(parameter[0][8],parameter[1][8],
    #                                                     parameter[3][8],parameter[4][8],
    #                                                     parameter[5][8]))
    # print("Parameter für Aufgabenteil C d3): {}A,{}V,{}mW,{},{}".format(parameter[0][9],parameter[1][9],
    #                                                     parameter[3][9],parameter[4][9],
    #                                                     parameter[5][9]))
    # print("Parameter für Aufgabenteil C d4): {}A,{}V,{}mW,{},{}".format(parameter[0][10],parameter[1][10],
    #                                                     parameter[3][10],parameter[4][10],
    #                                                     parameter[5][10]))

    # # # A) 6 modules
    # niceplot(x=UCA1, y=ICA1, c='tab:blue', ls='', marker='.', xaxis='$V$/V',
    #           plot5=True, x5=UCA1[parameter[2][0]], y5=ICA1[parameter[2][0]],
    #           marker5='v', c5='k', ls5='',
    #           yaxis=r'I/A', legend=True, safefig=True, safename='C1',
    #           plotlabel='measured data', titel='Light curve for setup of 6 cells')
    
    # # B) 6 modules with resistances
    # niceplot(x=UCB2, x2=UCB3, x3=UCB4, y=ICB2, y2=ICB3, y3=ICB4, plotlabel='R$_S$=68 $\Omega$, R$_P$=1 $\Omega$',   
    #           plotlabel2='R$_S$=1$\Omega$, R$_P$=68 $\Omega$',
    #           plotlabel3='R$_S$=R$_P$=3.3 $\Omega$',
    #           plot2=True, plot3=True, legend=True, xaxis='$V$/V', 
    #           yaxis='I/A', safefig=True, safename='C2', c='tab:green', 
    #           c2='tab:red', c3='tab:blue', ls='' ,plot5=True, 
    #           x5=np.array([UCB2[parameter[2][1]], UCB3[parameter[2][2]], UCB4[parameter[2][3]]]), 
    #           y5=np.array([ICB2[parameter[2][1]], ICB3[parameter[2][2]], ICB4[parameter[2][3]]]),
    #           marker5='v', c5='k', ls5='',
    #           ls2='', ls3='', marker='.', marker2='.', marker3='.',
    #           titel='Light curves for different wirings')
    
    # # # C) different shadings    
    niceplot(x=UCC1, x2=UCC2, x3=UCC3, y=ICC1, y2=ICC2, y3=ICC3, 
              plotlabel='1/3 shading, setup b)',   
              plotlabel2='1/2 shading, setup c)',
              plotlabel3='1/2 shading, setup a)',
              plot2=True, plot3=True, legend=True, xaxis='$V$/V', 
              yaxis='I/A', safefig=True, safename='C3', c='tab:green', 
              c2='tab:red', c3='tab:blue', ls='',
              ls2='', ls3='', marker='.', marker2='.', marker3='.', plot5=True, 
              x5=np.array([UCC3[parameter[2][4]], UCC2[parameter[2][5]], UCC1[parameter[2][6]]]), 
              y5=np.array([ICC3[parameter[2][4]], ICC2[parameter[2][5]], ICC1[parameter[2][6]]]),
              marker5='v', c5='k', ls5='',
              titel='Light curves for different shadings')
    # # # D)
    # # Light curve for 11 solar cells with/without ventilator
    # niceplot(x=UCD1, x2=UCD2, y=ICD1, y2=ICD2,
    #           plotlabel='ventilator',   plotlabel2='without ventilator',
    #           plot2=True, legend=True, xaxis='$V$/V', 
    #           yaxis='I/A', safefig=True, safename='C4', c='tab:green', 
    #           c2='tab:red', ls='', ls2='', marker='.', marker2='.', plot5=True, 
    #           x5=np.array([UCD1[parameter[2][7]], UCD2[parameter[2][8]]]), 
    #           y5=np.array([ICD1[parameter[2][7]], ICD2[parameter[2][8]]]),
    #           marker5='v', c5='k', ls5='',
    #           titel='Light curves for 11 solar cells measurement')
    # # Dark curve for 11 solar cells with ventilator
    # niceplot(x=UCD3, y=ICD3*1e2, plotlabel='ventilator',
    #           xaxis='$V$/V',  yaxis='I/mA',
    #           safefig=True, safename='C42_full', ls='', marker='.', 
    #           c='tab:green', legend=True, plot5=True, 
    #           x5=np.array([UCD3[parameter[2][9]]]), 
    #           y5=np.array([ICD3[parameter[2][9]]])*1e2,
    #           marker5='v', c5='k', ls5='',
    #           titel='Dark curve for 11 solar cells measurement')
    # # Dark curve for 11 solar cells without ventilator
    # niceplot(x=UCD4, y=ICD4, c='tab:green', plotlabel='without ventilator',
    #           safefig=True, safename='C42_empty', ls='', marker='.',
    #           x5=UCD4[parameter[2][10]], y5=ICD4[parameter[2][10]],
    #           plot5=True, marker5='v', c5='k', ls5='', xaxis='$V$/V', 
    #           yaxis='I/A', legend=True, titel='Dark curve for 11 solar cells measurement'
    #           )
    
    # """D: inluence of temperature"""
    # popt, cov = fit(linearfit, T[1:-2], V_OC[1:-2], sigma=None, absolute_sigma=True)
    # print("Die Parameter des lin-Fits sind:")
    # print("A+-∆A={}+-{}".format(popt[0], np.sqrt(np.diag(cov))[0]))
    # print("B+-∆B={}+-{}".format(popt[1], np.sqrt(np.diag(cov))[1]))
    # niceplot(T, V_OC, c='tab:green', ls='', marker='x', xaxis='T/$^{\circ}$C', 
    #           yaxis='V$_{OC}$/mV', safefig=True, plotlabel='measured data', 
    #           plot2=True, c2='k', plotlabel2=' fit to A$\cdot$T + B', 
    #           x2 = np.linspace(30,65, 300), y2=linearfit(np.linspace(30,65, 300), *popt),
    #           legend=True, titel='Temperature dependence of $V_{OC}$', 
    #           safename='D', ls2='-.', lw2=2.0, error=True, yerr=dT)
    
    # extra plot if neccessary
    # niceplot(x=UD1, x4=UD4[0:-3], y=ID1, y4=ID4[0:-3],
    #           plotlabel='T=$30^{\circ}$C',   plotlabel4='T=$65^{\circ}$C',
    #           plot4=True, legend=True, xaxis='$V$/V', 
    #           yaxis='I/A', safefig=True, safename='D2', c='tab:green', 
    #           c4='tab:orange', ls=''
    #           ,ls4='', marker='.', marker4='.',
    #           titel='Light curves for temperature extrema')
    
    # """E: Diffuse lighting"""
    # popt, cov = fit(cosfit, angle, I_OC/np.max(I_OC), p0 = [197.5, 0.001], sigma=None, 
    #             absolute_sigma=True)
    # A = 1 # in cm^2
    # print("Die Parameter des cos Fits sind:")
    # print("A+-∆A={}+-{}".format(popt[0], np.sqrt(np.diag(cov))[0]))
    # print("B+-∆B={}+-{}".format(popt[1], np.sqrt(np.diag(cov))[1]))
    # niceplot(angle, I_OC/np.max(I_OC), marker='v', ls='', ls2='-.', lw2=3, 
    #       c='tab:green', plot2=True, x2=np.linspace(0,100,901),
    #       y2=cosfit(np.linspace(0,100,901), *popt)/np.max(cosfit(np.linspace(0,100,901), *popt)), c2='k',  
    #       plotlabel2='fit to A $\cdot \cos($B$ \cdot \Phi)$', legend=True, 
    #       safefig=True, yaxis=r'I/max(I)',
    #       titel='Angular dependence of the current I($\Phi$)', 
    #       xaxis='$\Phi$/$^{\circ}$', plotlabel='measured data', safename='E1')
    # niceplot(x=np.cos(angle*np.pi/180), y=I_OC/np.max(I_OC), marker='v', ls='', 
    #       ls2='-.', lw2=3, 
    #       c='tab:green', plot2=True, x2=np.cos(angle*np.pi/180),
    #       c2='k',  y2=np.cos(angle*np.pi/180),
    #       plotlabel2='effective area', legend=True, safefig=True, 
    #       titel='Angular dependence of the current I($\Phi$)', 
    #       xaxis='$\cos(\Phi)$', 
    #       yaxis=r'I/max(I)',
    #       plotlabel='measured data', safename='E2')
    
    plt.show()
    
if __name__ == "__main__":
    main()
    