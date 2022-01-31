import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
from scipy.signal import find_peaks
import scipy.constants as sc
from praktikum import niceplot

"""
Auswertung FPraktikum Gaslaser, Fabry-Perot-Interferometer
"""
# Determine background noise of the measured signals
def noisecalc(x, limit):
    """
    :parameters:
    x: data
    limit: up to which value data is sconsidered noise
    :returns: mean noise value
    """
    noise = []
    for i in x:
        if i < limit:
            noise.append(i)
        else:
            i =+ 1
    return np.mean(np.asarray(noise))


def fitparameters_calc(x, y, p0=None, name=None):
    """
    Calculate gaussian distribution fit parameters
    """
    popt, cov = fit(gaussian, x, y, p0, absolute_sigma=True)
    print("GAUSSIAN PARAMETERS {}:".format(name))
    print("A+-∆A={}+-{}".format(popt[0], np.sqrt(np.diag(cov))[0]))
    print("mu+-∆mu={}+-{}".format(popt[1], np.sqrt(np.diag(cov))[1]))
    print("sigma+-∆sigma={}+-{}".format(popt[2], np.sqrt(np.diag(cov))[2]))
    print("FWHM: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt[2],
                                              2 * np.sqrt(2 * np.log(2)) *
                                              np.sqrt(np.diag(cov))[2]))

    return popt[0], popt[1], popt[2], np.sqrt(np.diag(cov))[0], np.sqrt(np.diag(cov))[1], np.sqrt(np.diag(cov))[2]


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

def lorentzian(x, x0, gamma, A):
    """
    Parameters
    ----------
    x : array, abscissa data.
    x0 : integer, resonance frequency.
    gamma : integer, linewdith parameter.
    A: integer, arbitrary amplitude.
    """
    return A/(np.pi*gamma*(1 + ((x - x0)/(gamma))**2))

"""
Faserspektrum
"""
lambd, I = np.loadtxt('Faserspec.txt', usecols=(0,1), unpack=True, skiprows=1)
popt, cov = fit(gaussian, lambd[960:980], I[960:980] / np.max(I), p0 = [1, 630, 10], absolute_sigma=True)
fitparameters_calc(lambd[960:980], I[960:980] / np.max(I), p0 = [1, 630, 10], name='faser')
# niceplot(c='tab:blue', x=lambd, y=I / np.max(I), plotlabel='Spektrum', titel='Spektrum des offenen HeNe-Lasers',
#          xaxis='$\lambda$ / nm', yaxis='$I / I_0$', xlim=(620,645), safefig=True,
#          size=(10,5), safename='faserspec_gauss', ls='', marker='.', lw2=2,
#          plot2=True, c2='k', ls2='--', x2=np.linspace(625, 640, 1000), y2=gaussian(np.linspace(625, 640, 1000), *popt),
#          plotlabel2='Gauß-fit', legend=True)

# Unsicherheit ∆lambda=0.5nm (Auflösung 1nm)
print("Freier Laser Wellenlänge: lambda +- ∆lambda = ({} +- {}) nm".format(lambd[np.argmax(I)], 0.5*1))
print("Freier Laser Frequenz: nu +- ∆nu = ({} +- {}) 1/m".format(sc.c/(lambd[np.argmax(I)]*1e-9),
                                                                 sc.c*0.5*1e-9/(lambd[np.argmax(I)]*1e-9)**2))

"""
FPI: Kommerzieller HeNe-Laser, FSR, longitudinale Moden mit versch. Polarisationen
"""
# Daten einlesen
t1, I1 = np.loadtxt('mitPol.txt', usecols=(0, 1), unpack=True, skiprows=4)
t2, I2 = np.loadtxt('ohnePol.txt', usecols=(0, 1), unpack=True, skiprows=4)
# get peaks for gaussian fit
indices1, parameters1 = find_peaks(I1, height=50)
indices2, parameters2 = find_peaks(I2, height=50)
# get noise
noise1 = noisecalc(I1, 50)
noise2 = noisecalc(I2, 8)
print("Hintergrund = {}mV".format(np.array([noise1, noise2])))
print("Hintergrund Schwelle = {}".format(np.array([50,8])))
niceplot(x=t1[indices1], y=I1[indices1] - noise1,
         x2=t1, y2=I1 - noise1, c='tab:red', c2='tab:blue',
         plot2=True, ls='', marker='x', plotlabel='Peaks', plotlabel2='Polarisationswinkel $\Theta=\pi/2$',
         titel='Polarisation der Moden im kommerziellen Laser', legend=True,
         xaxis='Kanal / u', yaxis='Spannung U / mV', size=(10,5), safefig=True, safename='mitpol', xlim=(0,600)
         )
popt7, cov7 = fit(gaussian, t2[60:120], I2[60:120] - noise2, p0=[122, 95, 10])
popt8, cov8 = fit(gaussian, t2[160:220], I2[160:220] - noise2, p0=[60, 190, 10])
popt9, cov9 = fit(gaussian, t2[380:460], I2[380:460] - noise2, p0=[125, 415, 10])
popt10, cov10 = fit(gaussian, t2[480:540], I2[480:540] - noise2, p0=[60, 515, 10])
fitparameters_calc(t2[60:120], I2[60:120] - noise2, p0=[122, 95, 10], name='1st peak')
fitparameters_calc(t2[160:220], I2[160:220] - noise2, p0=[60, 190, 10], name='2nd peak')
fitparameters_calc(t2[380:460], I2[380:460] - noise2, p0=[125, 415, 10], name='3rd peak')
fitparameters_calc(t2[480:540], I2[480:540] - noise2, p0=[60, 515, 10], name='4th peak')
niceplot(x=t2[indices2], y=I2[indices2] - noise2,
         x2=t2, y2=I2 - noise2, c='tab:red', c2='tab:blue',
         plot2=True, ls='', marker='x', plotlabel='Peaks', plotlabel2='Messwerte',
         titel='Longitudinales Modenspektrum des kommerziellen Lasers', legend=True,
         xaxis='Kanal / u', yaxis='Spannung U / mV', size=(10,5), safefig=True, safename='ohnepol',
         x3=np.linspace(-np.sqrt(2 * np.log(2)) * popt7[2] + popt7[1], np.sqrt(2 * np.log(2)) * popt7[2] + popt7[1], len(t2[60:120])),
         y3=max(gaussian(t2[60:120], *popt7))*np.ones(len(t2[60:120])) / 2, c3='tab:red', plot3=True, plotlabel3='FWHM',
         plot4=True, y4=max(gaussian(t2[160:220], *popt8))*np.ones(len(t2[160:220])) / 2,
         x4= np.linspace(-np.sqrt(2 * np.log(2)) * popt8[2] + popt8[1], np.sqrt(2 * np.log(2)) * popt8[2] + popt8[1], len(t2[160:220])),
         c4='tab:red',
         plot5=True, y5=max(gaussian(t2[380:460], *popt9))*np.ones(len(t2[380:460])) / 2,
         x5= np.linspace(-np.sqrt(2 * np.log(2)) * popt9[2] + popt9[1], np.sqrt(2 * np.log(2)) * popt9[2] + popt9[1], len(t2[380:460])),
         c5='tab:red',
         plot6=True, y6=max(gaussian(t2[480:540], *popt10))*np.ones(len(t2[480:540])) / 2,
         x6= np.linspace(-np.sqrt(2 * np.log(2)) * popt10[2] + popt10[1], np.sqrt(2 * np.log(2)) * popt10[2] + popt10[1], len(t2[480:540])),
         c6='tab:red', xlim=(0,600)
         )


u = sc.c/(2*0.075) / (np.mean(np.array([t2[indices2][3] - t2[indices2][1], t2[indices2][2] - t2[indices2][0]])))
du = u*np.sqrt((sc.c*0.01/(2*(0.075)**2)/(sc.c/(2*0.075)))**2 + ((1)/(324.5))**2)
print("Theoretischer FSR = ({}+-{})1/s".format(sc.c/(2*0.075), sc.c*0.01/(2*(0.075)**2)))  # FPI: FSR=c/2d
print("dt = {}".format(np.array([t2[indices2][3] - t2[indices2][1], t2[indices2][2] - t2[indices2][0]])))
print("mean FSR = ({} +- {})u".format(np.mean(np.array([t2[indices2][3] - t2[indices2][1], t2[indices2][2] - t2[indices2][0]])),
                         np.absolute((t2[indices2][3] - t2[indices2][1]) - (t2[indices2][2] - t2[indices2][0])) / 2)) # max-min/2
print("1u = ({} +- {})Hz".format(sc.c/(2*0.075) / (np.mean(np.array([t2[indices2][3] - t2[indices2][1], t2[indices2][2] - t2[indices2][0]]))),
                                 u*np.sqrt((sc.c*0.01/(2*(0.075)**2)/(sc.c/(2*0.075)))**2 +
                                ((1)/(324.5))**2)))
print("UNBEKANNTE RESONATORLÄNGE:")
print("dnu={}u".format(np.array([indices2[1]-indices2[0], indices2[3]-indices2[2]])))
print("mean dnu = ({} +- {})u".format(np.mean(np.array([indices2[1]-indices2[0], indices2[3]-indices2[2]])),
      (np.abs((indices2[1]-indices2[0]) - (indices2[3]-indices2[2])))/2))
print("RESONATORLÄNGE L = ({}+-{})".format(sc.c/(2*np.mean(np.array([indices2[1]-indices2[0], indices2[3]-indices2[2]])*u)),
                            sc.c/(2*np.mean(np.array([indices2[1]-indices2[0],
                            indices2[3]-indices2[2]])*u)) * np.sqrt(((du)/(u))**2 + (0.5/(97.5))**2)))

"""
Offener HeNe-Laser longitudinale Modenbetrachtung, versch. Resonatorlänge
"""
# Daten einlesen
t3, I3 = np.loadtxt('4peak.txt', usecols=(0, 1), unpack=True, skiprows=5) # schlechte Peaks
t4, I4 = np.loadtxt('4peak2.txt', usecols=(0, 1), unpack=True, skiprows=5) # Bessere Peaks
t5, I5 = np.loadtxt('4peak61.txt', usecols=(0, 1), unpack=True, skiprows=5)
# get peaks for gaussian fit
# indices3, parameters3 = find_peaks(I3, height=25) schlechte daten
indices4, parameters4 = find_peaks(I4, height=40, distance=16)
indices5, parameters5 = find_peaks(I5, height=40, distance=16)
# estimate noise
noise4, noise5 = noisecalc(I4, 32), noisecalc(I5, 30)
print("Hintergrund = {}mV".format(np.array([noise4, noise5])))
print("Hintergrund Schwelle = {}".format(np.array([30,30])))
# Gauss-fits
# popt3, cov3 = fit(gaussian, t3[indices3], I3[indices3] # schlechte daten
                 # )
fitparameters_calc(t4[indices4[0:4]], I4[indices4[0:4]] - noise4, p0=[150, 139, 100], name='41')
fitparameters_calc(t4[indices4[4:8]], I4[indices4[4:8]] - noise4, p0=[150, 434, 100], name='42')
fitparameters_calc(t5[indices5[0:4]], I5[indices5[0:4]] - noise5, p0=[150, 133, 100], name='51')
fitparameters_calc(t5[indices5[4:8]], I5[indices5[4:8]] - noise5, p0=[150, 454, 100], name='52')
popt41, cov41 = fit(gaussian, t4[indices4[0:4]], I4[indices4[0:4]] - noise4, p0=[150, 139, 100],
                 absolute_sigma=True)
popt42, cov42 = fit(gaussian, t4[indices4[4:8]], I4[indices4[4:8]] - noise4, p0=[150, 434, 100],
                 absolute_sigma=True)
popt51, cov51 = fit(gaussian, t5[indices5[0:4]], I5[indices5[0:4]] - noise5, p0=[150, 133, 100],
                 absolute_sigma=True)
popt52, cov52 = fit(gaussian, t5[indices5[4:8]], I5[indices5[4:8]] - noise5, p0=[150, 454, 100],
                 absolute_sigma=True)
sigma = np.array([popt41[2], popt42[2], popt51[2], popt52[2]])
gamma = 2 * np.sqrt(2 * np.log(2)) * sigma # FWHM for gaussian distributions
dsigma = np.array([np.sqrt(np.diag(cov41))[2], np.sqrt(np.diag(cov42))[2],
                   np.sqrt(np.diag(cov51))[2], np.sqrt(np.diag(cov52))[2] ])
dgamma = 2 * np.sqrt(2 * np.log(2)) * dsigma

# Plot
niceplot(x=t4, x2=t3, y=I4 - noise4, y2=I3, plotlabel=r'L=80cm second image',
         plotlabel2='L=80cm ', plotlabel3='L=61cm',
         plot2=False, size=(10, 5), c='tab:blue', c2='tab:olive', legend=True, plot4=True,
         plot5=True, x5=t4[indices4], y5=I4[indices4] - noise4, c5='tab:red', ls5='', marker5='x', plotlabel5='peaks',
         x4=np.linspace(0,300, 1000), y4=gaussian(np.linspace(0,300,1000), *popt41), c4='tab:green',
         x6=np.linspace(300,600,1000), plot6=True, c6='tab:green', y6=gaussian(np.linspace(300,600,1000), *popt42),
         plotlabel6='Gauß-fit'
         )
# Bei Breiteren Peaks: Besserer Fit möglich
niceplot(x=t5, y=I5 - noise5, c='tab:blue', legend=True, plotlabel='L=61cm',
         x2=np.linspace(0,300,1000), y2=gaussian(np.linspace(0,300,1000), *popt51), plot2=True, c2='tab:green',
         x3=t5[indices5], y3=I5[indices5] - noise5, plot3=True, ls3='', marker3='x', c3='tab:red',
         plot4=True, x4=np.linspace(300,600,1000), y4=gaussian(np.linspace(300,600,1000), *popt52), c4='tab:green',
         plotlabel4='Gauß-fit', size=(10, 5)
         )