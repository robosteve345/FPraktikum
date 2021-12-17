import numpy as np
from scipy.optimize import curve_fit as fit
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib.ticker as ticker
# import pandas as pd
    
########################################################
"""Nützliche Programme für das phyikalische Praktikum"""
########################################################

# Konstanten
h = 6.6265*10**(-34) # in Js
eV = 1.6*10**(-19) # in J
c0 = 2.9979*10**(8) # in m/s
rad = np.pi/180

"""Auswertungsprogramme zur linearen Regression"""
# Steven neue Version
def lin_fit(x, m, t):
    """ Lineare Fitfunktion
    -----------------------
    x: Abszissenarray
    m: Steigung
    t: Ordinatenabschnitt
    """
    return m*x + t


def lin_fit2(x, m):
    """ Lineare Fitfunktion
    -----------------------
    x: Abszissenarray
    m: Steigung
    t: ORdinatenabschnitt
    """
    return m*x 

def lin_reg2(x, y, dy, sigma_y, Titel=None, dx=None, sigma_x=None, p=None, 
              xaxis=None, yaxis=None):
    """Berechnet linearen fit zu gegebenen 2D-Datensatz in x und y. Zeichnet
    Fehlerschläuche, Hilfsgeraden und berechnet syst. und stat. MU der
    Parameters m der fit-Funktion:        y = m*x 
    --------------
    x: Abszissenarray
    y: Ordinatenarray
    dx: syst. MU der Abszissenwerte. Default=None
    dy: syst. MU der Ordinatenwerte
    Title: Titel, als String einzugeben
    sigma_y: stat. MU der Ordinatenwerte
    sigma_x: stat. MU der Abszissenwerte. Default=None
    p: Initial guess für Parameter m: [m0]. Array-like
    Bei Bedarf könnte noch bounds mit angegeben werden
    -------------------------------------------------
    """
    
    # Bounds 
    # b = [[],[]]
    # Berechne Fit-Daten für zentralen Fit und Fehlerschläuche
    # mit Hilfe der syst. MUs
    par1, cov1 = fit(lin_fit2, x, y + dy, absolute_sigma=True, sigma=sigma_y) # Oben
    par2, cov2 = fit(lin_fit2, x, y, absolute_sigma=True, sigma=sigma_y, 
                      p0=p) # Zentral
    par3, cov3 = fit(lin_fit2, x, y - dy, absolute_sigma=True, sigma=sigma_y) # Unten
    
    par4, cov4 = fit(lin_fit2, x, y + sigma_y, absolute_sigma=True, sigma=sigma_y)
    par5, cov5 = fit(lin_fit2, x, y - sigma_y, absolute_sigma=True, sigma=sigma_y)
    
    # Ausgleichsgeraden
    # Steigungen
    m_max = (lin_fit2(x, *par3)[-1] - lin_fit2(x, *par1)[0]) / (x[-1] - x[0])
    m_min = (lin_fit2(x, *par1)[-1] - lin_fit2(x, *par3)[0]) / (x[-1] - x[0])
    dm = 0.5*(m_max-m_min) # syst. MU für Steigung m
    
    # Für stat. MU der Steigung:
    mmax = (lin_fit2(x, *par5)[-1] - lin_fit2(x, *par4)[0]) / (x[-1] - x[0])
    mmin = (lin_fit2(x, *par4)[-1] - lin_fit2(x, *par5)[0]) / (x[-1] - x[0])
    sigma_m2 = 0.5*(mmax-mmin)
    
    # m und t aus dem zentralen Fit und ihre stat. MUs via Covarianzmatrix
    sigma_m = np.sqrt(np.diag(cov2))[0]
    m = par2[0]
    # Ordinatenabschnitte der Ausgleichsgeraden 
    t_max = lin_fit2(x, *par1)[-1] - m_min*x[-1]  
    t_min = lin_fit2(x, *par3)[-1] - m_max*x[-1] 
    
    print("Die Daten wurden an die Funktion y = m*x gefittet")
    print("X = X_mean +- dstatX +- dsystX")
    print("m = {} +- {} +- {}".format(par2[0], sigma_m, dm))
    print("t_max = {}, t_min = {}".format(t_max, t_min))
    
    """Plotbereich"""
    plt.figure(figsize=(18,9))
    
    # Experimentelle Daten mit Fehlerbalken
    plt.plot(x, y, marker='x', c='k', ls='', ms=8, label = 'experimental data')
    plt.errorbar(x, y, yerr=sigma_y + dy,c='k', 
                  ls='', capsize=3.5, label='combined MU')
    
    """Bei Bedarf: (kombinierte?) MU in x-Richtung"""
    plt.errorbar(x, y, xerr=dx + sigma_x ,c='k', 
              ls='', capsize=3.5, label='combined MU')
              
    # Plot der linearen Fits
    plt.plot(x, lin_fit2(x, *par2), c='tab:red', lw=2, label='linear fit') # main fit
    plt.plot(x, lin_fit2(x, *par1), c='tab:green', ls = '--', 
                          label='error tubes', lw=0.5) # upper fit
    plt.plot(x, lin_fit2(x, *par3), c='tab:green', ls = '--', lw=0.5) # lower fit
    
    # Ausgleichsgeraden
    plt.plot(x, m_max*x + t_min, c='tab:blue', label = 'balance lines', lw=0.5) # Max.
    plt.plot(x, m_min*x + t_max, c='tab:blue', lw=0.5) # Min.
    
    # Weitere Plotbefehler
    plt.title("{}".format(Titel), size=18)
    plt.xlabel('{}'.format(xaxis), size=18) 
    plt.ylabel('{}'.format(yaxis), size=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=15)
    plt.show()
    plt.savefig("{}.png".format(Titel), dpi=300)
    
    return m, sigma_m, dm, sigma_m2


def lin_reg(x, y, dy, sigma_y, Titel=None, dx=None, sigma_x=None, p=None, 
            xaxis=None, yaxis=None):
    """Berechnet linearen fit zu gegebenen 2D-Datensatz in x und y. Zeichnet
    Fehlerschläuche, Hilfsgeraden und berechnet syst. und stat. MUs der
    Parameter m und ctder fit-Funktion:          y = m*x + t
    --------------
    x: Abszissenarray
    y: Ordinatenarray
    dx: syst. MU der Abszissenwerte. Default=None
    dy: syst. MU der Ordinatenwerte
    Title: Titel, als String einzugeben
    sigma_y: stat. MU der Ordinatenwerte
    sigma_x: stat. MU der Abszissenwerte. Default=None
    p: Initial guess für Parameter m, t: [m0, t0]. Array-like
    Bei Bedarf könnte noch bounds mit angegeben werden
    -------------------------------------------------
    """
    
    # Bounds 
    # b = [[],[]]
    # Berechne Fit-Daten für zentralen Fit und Fehlerschläuche
    # mit Hilfe der syst. MUs
    par1, cov1 = fit(lin_fit, x, y + dy, absolute_sigma=True, sigma=sigma_y) # Oben
    par2, cov2 = fit(lin_fit, x, y, absolute_sigma=True, sigma=sigma_y, 
                     p0=p) # Zentral
    par3, cov3 = fit(lin_fit, x, y - dy, absolute_sigma=True, sigma=sigma_y) # Unten

    
    # Ausgleichsgeraden
    # Steigungen
    m_max = (lin_fit(x, *par3)[-1] - lin_fit(x, *par1)[0]) / (x[-1] - x[0])
    m_min = (lin_fit(x, *par1)[-1] - lin_fit(x, *par3)[0]) / (x[-1] - x[0])
    dm = 0.5*(m_max-m_min) # syst. MU für Steigung m
    
    # m und t aus dem zentralen Fit und ihre stat. MUs via Covarianzmatrix
    sigma_m = np.sqrt(np.diag(cov2))[0]
    sigma_t = np.sqrt(np.diag(cov2))[1]
    m = par2[0]
    t = par2[1]
    
    # Ordinatenabschnitte der Ausgleichsgeraden
    t_max = lin_fit(x, *par1)[-1] - m_min*x[-1]  
    t_min = lin_fit(x, *par3)[-1] - m_max*x[-1] 
    dt = (t_max - t_min)*0.5
    
    print("Die Daten wurden an die Funktion y = m*x + t gefittet")
    print("X = X_mean +- dstatX +- dsystX")
    print("m = {} +- {} +- {}".format(par2[0], sigma_m, dm))
    print("t = {} +- {} +- {}".format(par2[1], sigma_t, dt))
    
    
    """Plotbereich"""
    plt.figure(figsize=(18,9))
    
    # Experimentelle Daten mit Fehlerbalken
    plt.plot(x, y, marker='.', c='k', ls='', ms=12, label = 'experimental data')
    plt.errorbar(x, y, yerr=sigma_y ,c='k', 
                 ls='', capsize=7, label='stat. MU', lw=2)
    
    """Bei Bedarf: (kombinierte?) MU in x-Richtung"""
    # plt.errorbar(x, y, xerr=dx+sigma_x ,c='k', 
    #           ls='', capsize=7, label='combined MU', lw=2)
              
    # Plot der linearen Fits
    plt.plot(x, lin_fit(x, *par2), c='tab:red', lw=2.5, label='linear fit') # main fit
    plt.plot(x, lin_fit(x, *par1), c='tab:green', ls = '--', 
                          label='error tubes', lw=2.5) # upper fit
    plt.plot(x, lin_fit(x, *par3), c='tab:green', ls = '--', lw=2.5) # lower fit
    
    # Ausgleichsgeraden
    plt.plot(x, t_min + m_max*x, c='tab:blue', label = 'balance lines', lw=2.0) # Max.
    plt.plot(x, t_max + m_min*x, c='tab:blue', lw=2.0) # Min.
    
    # Weitere Plotbefehler
    plt.title("{}".format(Titel), size=20)
    plt.xlabel('{}'.format(xaxis), size=20) 
    plt.ylabel('{}'.format(yaxis), size=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=17, markerscale=1.5)
    plt.show()
    plt.savefig("{}.png".format(Titel), dpi=300)
    
    return m, sigma_m, dm, t, sigma_t, dt


# Stevens alte Version mit linreg-Funktion aus Scipy
def lin_reg_oldversion(x_data, y_data, dy, dsyst_y, dx=None, dsyst_x=None):
    """Compute linear regression of 'x_data' and 'y_data', 
    with combined errorbars of 'dy, dsyst_y'. 
    Gives back slope, intercept and their stat. and syst. MUs. 
    Fits data to the function y = m*x+c.
    If needed, dx, dsyst_x can be manually implemented in the following.
    Plots stat. MUs as errorbars in y-direction.
    --------------------------------------------
    x_data: array-like
    y_data: array_like
    """
    (m, c, rvalue1, pvalue1, dstat_m) = linregress(x_data, y_data) # fit of the data
    (m1, c1, rvalue1, pvalue1, dstat_m1) = linregress(x_data, y_data + dsyst_y) # upper fit
    (m2, c2, rvalue2, pvalue2, dstat_m2) = linregress(x_data, y_data - dsyst_y) # lower fit
    
    fit = c + m*(x_data) # central fit
    fit1 = c1 + m1*(x_data) # upper fit
    fit2 = c2 + m2*(x_data) # lower fit
    # Ausgleichgeraden
    # Stiegungen m
    m_max = (fit1[-1] - fit2[0])/(x_data[-1]-x_data[0]) # calculate slopes
    m_min = (fit2[-1] - fit1[0])/(x_data[-1]-x_data[0])
    dsyst_m = (m_max - m_min)*0.5
    
    # Ordinatenabschnitte
    c_max = fit2[-1] - m_min*x_data[-1]  
    c_min = fit1[-1] - m_max*x_data[-1] 
    dsyst_c = (c_max - c_min)*0.5
    
    # Gebe St.abweichung des y-Achsenabschnitts an
    result = linregress(x_data,y_data)
    dstat_c = result.intercept_stderr  
    
    # Gebe Ergebnis aus
    print("Die Daten wurden an die Funktion y = m*x + c gefittet")
    print("m = {} +- {} (stat.) +- {} (syst)".format(m, dstat_m, dsyst_m))
    print("c = {} +- {} (stat.) +- {} (syst)".format(c, 
                                                      result.intercept_stderr, 
                                                      dsyst_c))
    
    # Plotbereich
    plt.figure(figsize=(18,9))
    
    # Experimentelle Daten mit Fehlerbalken
    plt.plot(x_data, y_data, marker='.', c='k', ls='', ms=4)
    plt.errorbar(x_data, y_data, yerr = dy ,c='tab:orange', 
              ls='', capsize=1.5, label='stat. MU')
    
    """Bei Bedarf: (kombinierte?) MU in x-Richtung"""
    # plt.errorbar(x_data, y_data, xerr=dx+dsyst_x ,c='tab:green', 
    #          ls='', capsize=2, label='combined MU')
              
    # Lineare Fits
    plt.plot(x_data, fit, c='tab:red', lw=1) # main fit
    plt.plot(x_data, fit1, c='tab:green', ls = '--', 
                          label='error tube', lw=0.9) 
    plt.plot(x_data, fit2, c='tab:green', ls = '--', lw=0.9) # main fit
    
    # Ausgleichsgeraden
    plt.plot(x_data, c_min + m_max*x_data, c='tab:blue', label = 'balance line', 
              lw=0.9) # Max.
    plt.plot(x_data, c_max + m_min*x_data, c='tab:blue', lw=0.9) # Min.
    
    # Weitere Plotbefehler
    plt.title("{}".format(Titel), size=15)
    plt.xlabel('x', size=15) 
    plt.ylabel('y', size=15)
    plt.legend()
    plt.savefig('Test.png', dpi=300)
    
    return m, dstat_m, dsyst_m, c, dstat_c, dsyst_c


# Lukas Koenigs Fitfunktion
def linear(x, a, b):
    return a*x+b


def linreg(x, y, ax=None, dx=None, dy=None, dysys=None, dxsys=None, plot=True):
    """
    führt lineare regression aus, mit minimal und maximalgerade für syst. 
    Fehler

    Parameters
    ----------
    x, y: Array-like, Daten
    ax : Plotbereich.
    dx, dy: statistische Fehler
    dxssys, dysys: systematischer Fehler in y-Richtung
    plot: bool, default True. Soll geplottet werden? Dann ax angeben!

    Returns
    -------
    popt, perr, da, db

    """
    # fit
    popt, pcov = curve_fit(linear, x, y, sigma = dy,
                        p0 = [0,0], absolute_sigma=True) 

    # Systematische Unsicherheiten
    if np.all(dxsys != None):
        if np.all(dysys != None):
            covh, v = curve_fit(linear, x + np.sign(popt[0])*dxsys, 
                                y - dysys, sigma=dy,
                                  p0 = [0,0], absolute_sigma=True)
            covl, v = curve_fit(linear, x - np.sign(popt[0])*dxsys, 
                                y + dysys, sigma=dy,
                                  p0 = [0,0], absolute_sigma=True)
        else:
            covh, v = curve_fit(linear, x - np.sign(popt[0])*dxsys, y, sigma=dy,
                                  p0 = [0,0], absolute_sigma=True)
            covl, v = curve_fit(linear, x+np.sign(popt[0])*dxsys, y, sigma=dy,
                                  p0 = [0,0], absolute_sigma=True)   
           
        lineh = linear(x, *covh)
        linel = linear(x, *covl)
        al = (-lineh[0]+linel[-1])/(x[-1]-x[0])
        ah = (-linel[0]+lineh[-1])/(x[-1]-x[0])
        da = (ah-al)/2
        linebl = linear(x[0], *covh) - al*x[0]
        linebh = linear(x[0], *covl) - ah*x[0]
        db=(linebl-linebh)/2
    elif np.all(dysys != None):
        covh, v = curve_fit(linear, x, y - dysys, sigma=dy,
                                  p0 = [0,0], absolute_sigma=True)
        covl, v = curve_fit(linear, x, y + dysys, sigma=dy,
                                  p0 = [0,0], absolute_sigma=True)
        lineh = linear(x, *covh)
        linel = linear(x, *covl)
        al = (-lineh[0]+linel[-1])/(x[-1]-x[0])
        ah = (-linel[0]+lineh[-1])/(x[-1]-x[0])
        da = (ah-al)/2
        linebl = linear(x[0], *covh) - al*x[0]
        linebh = linear(x[0], *covl) - ah*x[0]
        db=(linebl-linebh)/2
    else:
        da=db=0
        
    # Statistische Unsicherheiten
    perr = np.sqrt(np.diag(pcov))
    
    # Plot
    # if plot:
    #     xdraw = np.linspace(min(x)*1, max(x)*1.0, 20) 
    #     ax.errorbar(x, y, xerr=dx, yerr=dy, capsize=4, ls='', marker='x', 
    #                 label='Datapoints')
    #     ax.plot(xdraw, linear(xdraw,*popt), label='Regression')
    #     # Plot Fehlerschlauch 
    #     if np.all(dxsys != None) or np.all(dysys != None):
    #         ax.plot(xdraw, linear(xdraw, al, linebl), c='r', ls='--' ,
    #                  label='min/max slope')
    #         ax.plot(xdraw, linear(xdraw, ah, linebh), c='r', ls='--')
    #         ax.plot(xdraw, linear(xdraw,*covh), c='k', label='Borderlines')
    #         ax.plot(xdraw, linear(xdraw,*covl), c='k')
    return popt, perr, da, db


"""Plotfunktionen"""
def pinmp_ticks(axis, ticks):
    axis.set_major_locator(ticker.MaxNLocator(ticks))
    axis.set_minor_locator(ticker.MaxNLocator(ticks * 10))
    return axis


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
    
    if error == True:
        ax.errorbar(x, y, xerr=yerr, capsize=cs, c=c, ms=ms, ls=ls,
                     marker=marker, label = r'{}'.format(plotlabel))
    else:
        ax.plot(x, y, ms=ms, lw=lw, ls=ls, marker=marker, c=c, 
                 label=r'{}'.format(plotlabel))
    # ax.fill_between(x[0:4], y[0:4], color='red', alpha=.2)
    # ax.fill_between(x[3:8], y[3:8], color='green', alpha=.2)
    # ax.fill_between(x[7:12], y[7:12], color='blue', alpha=.2, )
    ax.set_ylim(0, max(y7+0.2e6))
    ax.set_xlim(-1250, 1250)
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
        ax.legend(fontsize=fs-4, markerscale=ms-2, facecolor='white')
    if safefig == True:
        plt.savefig('{}.svg'.format(safename), dpi=300)

    plt.show()


def nice_plot(x1, y1, c1='k', c2='k', fs=20, ts=25, yerr1=None, xerr1=None, yerr2=None, xerr2=None, 
              x2=None, y2=None, ax1=True, size=(18,9), ax2=None, 
              scatter=None, title=None, xaxis="x", yaxis="y", legend=None):
    """
    Plottet y gegen x mit adäuquaten Plotparametern für wiss. Arbeiten
    als Scatter-Plot oder zusammenhängen Graphen für bis zu zwei Plotbereiche.
    ------------------------------------------------
    Parameter:
        x1,2: Abszissenarray für Plotebereich 1,2
        y1,2: Ordinatenarray für Plotbereich 1,2
        size: Plotbereichgröße. Default: (18,9)
        ax1, ax2: Plotbereiche (Default: Nur ax1)
        yerr1,2: Fehlerbalken in Ordinatenrichtung, default: None
        xerr1,2: Fehlerbalken in Abszissenrichtung, default: None
        scatter: Daten als zusammenhängender Graph oder nicht (True oder False)
        title: Titel des Plots
        xaxis: x-Achsen Titel
        yaxis: y-Achsen Titel
        c: Farbe der Datenpunkte, Default: Schwarz (k)
        legend: Plotte legende. Default: None, plottet bei legend = True
    ------------------------------------------------
    Rückgabe:
        1-2 Plotbereiche mit Daten (samt Fehlerbalken)
    """
    if title is None:
        title = ""
    fig = plt.figure(figsize=(18,9)) 
    ax1 = fig.add_subplot(1, 2, 1)
    if ax2==True:
        ax2 = fig.add_subplot(1, 2, 2)
        if scatter==True:
            if legend==True:
                plt.legend()
            # Plotbereich 1
            plt.figure(figsize=size)
            #ax1.axis([np.min(x), np.max(x), np.min(v_werte), Emax])
            ax1.set_xlabel("{}".format(xaxis),fontsize=fs)
            ax1.set_ylabel("{}".format(yaxis),fontsize=fs)
            ax1.set_title("{}".format(title),fontsize=fs)
            plt.yticks(fontsize=ts)
            plt.xticks(fontsize=ts)
            ax1.errorbar(x1, y1, yerr=yerr1 ,color=c1, 
                 ls='', capsize=3.5, label='MU')
            ax1.errorbar(x1, y1, xerr=xerr1 ,color=c1, 
                 ls='', capsize=3.5, label='MU')
            ax1.scatter(x1, y1, linewidth=2, color=c1, label='')     
            # Plotbereich 2
            plt.figure(figsize=size)
            #ax1.axis([np.min(x), np.max(x), np.min(v_werte), Emax])
            ax2.set_xlabel("{}".format(xaxis),fontsize=fs)
            ax2.set_ylabel("{}".format(yaxis),fontsize=fs)
            ax2.set_title("{}".format(title))
            plt.yticks(fontsize=ts)
            plt.xticks(fontsize=ts)
            ax2.scatter(x2, y2, linewidth=2, color=c2, label='') 
            ax2.errorbar(x2, y2, yerr=yerr2 ,color=c2, 
                 ls='', capsize=3.5, label='MU')
            ax2.errorbar(x2, y2, xerr=xerr2 ,color=c2, 
                 ls='', capsize=3.5, label='MU')
    
        else:
            if legend==True:
                plt.legend()
            # Plotbereich 1
            plt.figure(figsize=size)
           # ax1.axis([np.min(x), np.max(x), np.min(v_werte), Emax])
            ax1.set_xlabel("{}".format(xaxis),fontsize=fs)
            ax1.set_title("{}".format(title),fontsize=fs)
            ax1.set_ylabel("{}".format(yaxis),fontsize=fs)
            plt.yticks(fontsize=ts)
            plt.xticks(fontsize=ts)
            ax1.plot(x1, y1, linewidth=2, color=c1, label='')  
            ax1.errorbar(x1, y1, yerr=yerr1 ,color=c1, 
                 ls='', capsize=3.5, label='MU')
            ax1.errorbar(x1, y1, xerr=xerr1 ,color=c1, 
                 ls='', capsize=3.5, label='MU')
            # Plotbereich 2
            plt.figure(figsize=size)
            #ax1.axis([np.min(x), np.max(x), np.min(v_werte), Emax])
            ax2.set_xlabel("{}".format(xaxis),fontsize=fs)
            ax2.set_ylabel("{}".format(yaxis),fontsize=fs)
            ax2.set_title("{}".format(title),fontsize=fs)
            plt.yticks(fontsize=ts)
            plt.xticks(fontsize=ts)
            ax2.plot(x2, y2, linewidth=2, color=c2, label='')   
            ax2.errorbar(x2, y2, yerr=yerr2 ,color=c2, 
                 ls='', capsize=3.5, label='MU')
            ax2.errorbar(x2, y2, xerr=xerr2 ,color=c2, 
                 ls='', capsize=3.5, label='MU')
            
    else:
        if scatter==True:
            if legend==True:
                plt.legend()
            # Plotbereich 1
            plt.figure(figsize=size)
            #ax1.axis([np.min(x1), np.max(x1), np.min(), Emax])
            ax1.set_xlabel("{}".format(xaxis),fontsize=fs)
            ax1.set_ylabel("{}".format(yaxis),fontsize=fs)
            ax1.set_title("{}".format(title),fontsize=fs)
            plt.yticks(fontsize=ts)
            plt.xticks(fontsize=ts)
            ax1.scatter(x1, y1, linewidth=2, color=c1, label='')  
            ax1.errorbar(x1, y1, yerr=yerr1 ,color=c1, 
                 ls='', capsize=3.5, label='MU')
            ax1.errorbar(x1, y1, xerr=xerr1 ,color=c1, 
                 ls='', capsize=3.5, label='MU')
    
        else:
            if legend==True:
                plt.legend()
            # Plotbereich 1
            plt.figure(figsize=size)
            #ax1.axis([np.min(x), np.max(x), np.min(v_werte), Emax])
            ax1.set_xlabel("{}".format(xaxis),fontsize=fs)
            ax1.set_ylabel("{}".format(yaxis), fontsize=fs)
            ax1.set_title("{}".format(title),fontsize=fs)
            ax1.plot(x1, y1, linewidth=2, color=c2, label='')     
            ax2.errorbar(x2, y2, yerr=yerr2 ,color=c2, 
                 ls='', capsize=3.5, label='MU')
            plt.yticks(fontsize=ts)
            plt.xticks(fontsize=ts)
            ax2.errorbar(x2, y2, xerr=xerr2 ,color=c2, 
                 ls='', capsize=3.5, label='MU')


# x = np.linspace(0, 2*np.pi, 1000)
# y1 = np.exp(-x)*np.cos(x)
# y2 = np.exp(-(x-4)**2)
# yerr1 = 0.01*y1
# yerr2 = 0.01*y2
# nice_plot(x, y1, c1='tab:red', x2=x, y2=y2, scatter=False, ax2=True, yerr1=yerr1, yerr2=yerr2)
        
# plt.figure(figsize=(18,9))
# plt.xlabel(r'', size=20)
# plt.ylabel(r'', size=20)
# plt.title(r'',size=20)
# plt.ylim(min() - ([]-[]),max() + ([]-[]]))
# plt.yticks(fontsize=14)
# plt.xticks(fontsize=14)
# # Fit
# plt.plot(xDSL, DS_fit(xDSL, *parameters4), label = 'Fitted function', 
#           c='tab:green', lw=2.5)
# plt.errorbar(xDSL, y, yerr=dy, ls="", capsize=4, lw=2, c='tab:blue')
# # experimental
# plt.plot(, , c='tab:blue', ls='', marker='x', ms=12, 
#           label = 'measured data') 
# plt.legend(fontsize=13)
# plt.savefig('XXX.png', dpi=300)
# # plt.show()
