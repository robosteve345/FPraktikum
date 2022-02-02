import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
from scipy.signal import find_peaks
import scipy.constants as sc
from praktikum import niceplot
from scipy.ndimage import gaussian_filter1d

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
    return np.mean(np.asarray(noise))


def fitparameters_calc(x, y, ordinate=False, p0=None, name=None):
    """
    Calculate gaussian distribution fit parameters
    """
    if ordinate == True:
        popt, cov = fit(gaussian1, x, y, p0, absolute_sigma=True)
    else:
        popt, cov = fit(gaussian2, x, y, p0, absolute_sigma=True)

    print("GAUSSIAN PARAMETERS {}:".format(name))
    print("A+-∆A={}+-{}".format(popt[0], np.sqrt(np.diag(cov))[0]))
    print("mu+-∆mu={}+-{}".format(popt[1], np.sqrt(np.diag(cov))[1]))
    print("sigma+-∆sigma={}+-{}".format(popt[2], np.sqrt(np.diag(cov))[2]))
    print("FWHM: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt[2],
                                              2 * np.sqrt(2 * np.log(2)) *
                                              np.sqrt(np.diag(cov))[2]))

    return popt[0], popt[1], popt[2], np.sqrt(np.diag(cov))[0], np.sqrt(np.diag(cov))[1], np.sqrt(np.diag(cov))[2]


def gaussian1(x, A, mu, sigma, B):
    """Gauss-function with arbitrary amplitude A (not normed).
    Parameters
    ----------
    x : array.
    A : integer, amplitude.
    mu : integer, expectation value.
    sigma : integer, standard deviation.
    B: noise
    """
    return A*np.exp(-(x-mu)**2/(2*sigma**2)) + B


def gaussian2(x, A, mu, sigma):
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
# popt, cov = fit(gaussian, lambd[960:980], I[960:980] / np.max(I), p0 = [1, 630, 10], absolute_sigma=True)
# fitparameters_calc(lambd[960:980], I[960:980] / np.max(I), p0 = [1, 630, 10], name='faser')
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
# niceplot(x=t1[indices1], y=I1[indices1] - noise1,
#          x2=t1, y2=I1 - noise1, c='tab:red', c2='tab:blue',
#          plot2=True, ls='', marker='x', plotlabel='Peaks', plotlabel2='Polarisationswinkel $\Theta=\pi/2$',
#          titel='Polarisation der Moden im kommerziellen Laser', legend=True,
#          xaxis='Kanal / u', yaxis='Spannung U / mV', size=(10,5), safefig=True, safename='mitpol', xlim=(0,600)
#          )
popt7, cov7 = fit(gaussian1, t2[60:120], I2[60:120], p0=[122, 95, 10, 10])
popt8, cov8 = fit(gaussian1, t2[160:220], I2[160:220], p0=[60, 190, 10, 10])
popt9, cov9 = fit(gaussian1, t2[380:460], I2[380:460], p0=[125, 415, 10, 10])
popt10, cov10 = fit(gaussian1, t2[480:540], I2[480:540], p0=[60, 515, 10, 10])
fitparameters_calc(t2[60:120], I2[60:120] , ordinate=True, p0=[122, 95, 10, 10], name='1st peak')
fitparameters_calc(t2[160:220], I2[160:220], ordinate=True,p0=[60, 190, 10, 10], name='2nd peak')
fitparameters_calc(t2[380:460], I2[380:460], ordinate=True,p0=[125, 415, 10, 10], name='3rd peak')
fitparameters_calc(t2[480:540], I2[480:540], ordinate=True,p0=[60, 515, 10, 10], name='4th peak')
u = sc.c/(2*0.075) / (np.mean(np.array([t2[indices2][3] - t2[indices2][1], t2[indices2][2] - t2[indices2][0]]))) # in Hz
du = u*np.sqrt((sc.c*0.01/(2*(0.075)**2)/(sc.c/(2*0.075)))**2 + ((1)/(324.5))**2)
A2 = np.array([popt7[0], popt8[0], popt9[0], popt10[0]])
dA2 = np.array([ np.sqrt(np.diag(cov7))[0], np.sqrt(np.diag(cov8))[0],
                   np.sqrt(np.diag(cov9))[0], np.sqrt(np.diag(cov10))[0] ])
mu2 = np.array([popt7[1], popt8[1], popt9[1], popt10[1]])*u
dmu2 = np.array([ np.sqrt(np.diag(cov7))[2], np.sqrt(np.diag(cov8))[2],
                   np.sqrt(np.diag(cov9))[2], np.sqrt(np.diag(cov10))[2]])*u
sigma2 = np.array([popt7[2], popt8[2], popt9[2], popt10[2]])*u
FWHM2 = 2 * np.sqrt(2 * np.log(2)) * sigma2
dsigma2 = np.array([np.sqrt(np.diag(cov7))[2], np.sqrt(np.diag(cov8))[2],
                   np.sqrt(np.diag(cov9))[2], np.sqrt(np.diag(cov10))[2] ])*u
dFWHM2 = 2 * np.sqrt(2 * np.log(2)) * dsigma2
print("GEWICHTETER WERT: FWHM=({}+-{})Hz".format(np.dot(FWHM2, dFWHM2), np.sqrt(1/(np.sum(dFWHM2**2)))))
print("PARAMETER KOMMERZIELLER LASER")
print("A={}mV".format(A2))
print("dA={}mV".format(dA2))
print("mu={}Hz".format(mu2))
print("dmu={}Hz".format(dmu2))
print("sigma2={}Hz".format(sigma2))
print("dsigma2={}Hz".format(dsigma2))
print("FWHM={}Hz".format(FWHM2))
print("dFWHM={}Hz".format(dFWHM2))
print("sigma2+-sigma2=({}+-{})Hz".format(np.mean(sigma2), np.mean(dsigma2)))
print("FWHM+-FWHM=({}+-{})Hz".format(np.mean(FWHM2), np.mean(dFWHM2)))
# niceplot(x=t2[indices2], y=I2[indices2],
#          x2=t2, y2=I2, c='tab:red', c2='tab:blue',
#          plot2=True, ls='', marker='x', plotlabel='Peaks', plotlabel2='Messwerte',
#          titel='Longitudinales Modenspektrum des kommerziellen Lasers', legend=True,
#          xaxis='Kanal / u', yaxis='Spannung U / mV', size=(10,5), safefig=True, safename='ohnepol',
#          x3=np.linspace(-np.sqrt(2 * np.log(2)) * popt7[2] + popt7[1], np.sqrt(2 * np.log(2)) * popt7[2] + popt7[1], len(t2[60:120])),
#          y3=max(gaussian1(t2[60:120], *popt7))*np.ones(len(t2[60:120])) / 2, c3='tab:red', plot3=True, plotlabel3='FWHM',
#          plot4=True, y4=max(gaussian1(t2[160:220], *popt8))*np.ones(len(t2[160:220])) / 2,
#          x4= np.linspace(-np.sqrt(2 * np.log(2)) * popt8[2] + popt8[1], np.sqrt(2 * np.log(2)) * popt8[2] + popt8[1], len(t2[160:220])),
#          c4='tab:red',
#          plot5=True, y5=max(gaussian1(t2[380:460], *popt9))*np.ones(len(t2[380:460])) / 2,
#          x5= np.linspace(-np.sqrt(2 * np.log(2)) * popt9[2] + popt9[1], np.sqrt(2 * np.log(2)) * popt9[2] + popt9[1], len(t2[380:460])),
#          c5='tab:red',
#          plot6=True, y6=max(gaussian1(t2[480:540], *popt10))*np.ones(len(t2[480:540])) / 2,
#          x6= np.linspace(-np.sqrt(2 * np.log(2)) * popt10[2] + popt10[1], np.sqrt(2 * np.log(2)) * popt10[2] + popt10[1], len(t2[480:540])),
#          c6='tab:red', xlim=(0,600)
#          )
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
print("RESONATORLÄNGE L = ({}+-{})m".format(sc.c/(2*np.mean(np.array([indices2[1]-indices2[0], indices2[3]-indices2[2]])*u)),
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
I41 = gaussian_filter1d(I4, sigma=0.5)  # make graph smooth by gaussian filter
I51 = gaussian_filter1d(I5, sigma=0.5)  # make graph smooth by gaussian filter
indices4, parameters4 = find_peaks(I41, height=30, distance=16)
indices5, parameters5 = find_peaks(I51, height=36, distance=16)
indices = np.concatenate([np.asarray(indices4)+1, np.asarray(indices5)+1])
print(indices)
print("x1={}".format(t4[indices4]))
print("x2={}".format(t5[indices5]))
print("y1={}".format(I41[indices4]))
print("y2={}".format(I51[indices5]))
##############
# estimate noise
# noise4, noise5 = noisecalc(I4, 32), noisecalc(I5, 30)
# print("Hintergrund = 4:{}mV, 5:{}mV".format(noise4, noise5))
# print("Hintergrund Schwelle = {}".format(np.array([20, 30, 30])))
##############
# Gauss-fits
# Vorgehensweise: für 80cm: fit für jeweils letzten 3 peaks der 4er-peak-serien, addiere händisch ordinatenabschnitt,
# weil curve_fit sonst nicht mitmacht.
# Für 61cm: normaler fit mit ordinatenabschnitt
fitparameters_calc(t4[indices4[0:4]], I41[indices4[0:4]], ordinate=True, p0=[120, 136, 20, 20], name='41')
fitparameters_calc(t4[indices4[4:8]], I41[indices4[4:8]], ordinate=False, p0=[150, 434, 20], name='42')
# GAUSS FIT WIRD NUR FÜR EINE RESONATORLÄNGE BENÖTIGT
fitparameters_calc(t5[indices5[0:4]], I51[indices5[0:4]], ordinate=True, p0=[110, 133, 20, 20], name='51')
fitparameters_calc(t5[indices5[4:8]], I5[indices5[4:8]],  ordinate=False, p0=[150, 454, 20], name='52')
popt41, cov41 = fit(gaussian1, t4[indices4[0:4]], I41[indices4[0:4]], p0=[140, 136, 28, 20],
                 absolute_sigma=True)
popt42, cov42 = fit(gaussian1, t4[indices4[4:8]], I41[indices4[4:8]], p0=[95, 425, 20, 20],
                 absolute_sigma=True)

popt51, cov51 = fit(gaussian1, t5[indices5[0:4]], I51[indices5[0:4]], p0=[110, 133, 20, 20],
                 absolute_sigma=True)
popt52, cov52 = fit(gaussian2, t5[indices5[4:8]], I51[indices5[4:8]], p0=[150, 455, 20],
                 absolute_sigma=True)
A = np.array([popt41[0], popt42[0], popt51[0], popt52[0]])
dA = np.array( [np.sqrt(np.diag(cov41))[0], np.sqrt(np.diag(cov42))[0],
                   np.sqrt(np.diag(cov51))[0], np.sqrt(np.diag(cov52))[0]])
mu = np.array([popt41[1], popt42[1], popt51[1], popt52[1]])*u
dmu = np.array( [np.sqrt(np.diag(cov41))[1], np.sqrt(np.diag(cov42))[1],
                   np.sqrt(np.diag(cov51))[1], np.sqrt(np.diag(cov52))[1]])
sigma = np.array([popt41[2], popt42[2], popt51[2], popt52[2]])*u
FWHM = 2 * np.sqrt(2 * np.log(2)) * sigma # FWHM for gaussian distributions
dsigma = np.array([np.sqrt(np.diag(cov41))[2], np.sqrt(np.diag(cov42))[2],
                   np.sqrt(np.diag(cov51))[2], np.sqrt(np.diag(cov52))[2] ])*u
dFWHM = 2 * np.sqrt(2 * np.log(2)) * dsigma
print("A={}".format(A))
print("dA={}".format(dA))
print("mu={}Hz".format(mu))
print("sigma={}Hz".format(sigma))
print("dsigma={}Hz".format(dsigma))
print("FWHM={}Hz".format(FWHM))
print("dFWHM={}Hz".format(dFWHM))
print("sigma+-sigma=({}+-{})Hz".format(np.mean(sigma), np.mean(dsigma)))
print("FWHM+-FWHM=({}+-{})Hz".format(np.mean(FWHM), np.mean(dFWHM)))
# Plot
# 80 cm plot 2(bessere daten)
# niceplot(x=t4, y=I4 - noise4,
#          plotlabel='L=80cm (2)', plotlabel3='L=61cm',
#          plot2=False, size=(10, 5), c='tab:blue', c2='tab:olive', legend=True, plot4=True,
#          plot5=True, x5=t4[indices4], y5=I41[indices4], c5='tab:red', ls5='', marker5='x', plotlabel5='peaks',
#          x4=np.linspace(0,300, 1000), y4=gaussian1(np.linspace(0,300,1000), *popt41), c4='tab:green',
#          x6=np.linspace(300,600,1000), plot6=True, c6='tab:orange', y6=gaussian1(np.linspace(300,600,1000), *popt42),
#          plotlabel6='Gauß-fit', xlim=(0,600), plot7=True, x7=np.linspace(300,600,1000), y7=gaussian1(np.linspace(300,600,1000),
#                                                                                                   95, 430, 28, 16), c7='tab:olive',
#          plot3=True, x3=t4, y3=I41, c3='k'
#          )

# 41: 100, 135, 30, 18, 42: 95, 430, 28, 16, 51:
# # Bei Breiteren Peaks: Besserer Fit möglich
# niceplot(x=t5, y=I5, c='tab:blue', legend=True, plotlabel='L=61cm',
#          x2=np.linspace(0,300,1000), y2=gaussian1(np.linspace(0,300,1000), *popt51), plot2=True, c2='tab:green',
#          x3=t5[indices5], y3=I51[indices5], plot3=True, ls3='', marker3='x', c3='tab:red',
#          plot4=True, x4=np.linspace(300,600,1000), y4=gaussian2(np.linspace(300,600,1000), *popt52), c4='tab:green',
#          plot5=True, x5=t5, y5=I51, lw5=2.3, c5='k',
#          plotlabel4='Gauß-fit', size=(10, 5), xlim=(0,600), plot7=True, x7=np.linspace(0,300,1000), y7=gaussian2(np.linspace(0,300,1000),
#                                                                                                   140, 134, 27) + 20, c7='tab:olive',
#          plot6=True, x6=np.linspace(300,600,1000), y6=gaussian2(np.linspace(300,600,1000), 110, 452, 25) + 20, c6='tab:orange'
#          )
plt.show()