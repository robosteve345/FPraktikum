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


"""
Faserspektrum
"""
lambd, I = np.loadtxt('Faserspec.txt', usecols=(0,1), unpack=True, skiprows=1)
# niceplot(c='tab:blue', x=lambd, y=I/np.max(I), plotlabel='Spektrum', titel='Spektrum des offenen HeNe-Lasers',
#          xaxis='$\lambda$ / nm', yaxis='$I / I_0$', xlim=(550,700), safefig=True,
#          size=(10,5), y2=np.linspace(0, 1, 1000), safename='faserspec')
# Unsicherheit ∆lambda=0.5nm (Auflösung 1nm)
print("Freier Laser Wellenlänge: lambda +- ∆lambda = ({} +- {}) nm".format(lambd[np.argmax(I)], 0.5*1))
print("Freier Laser Frequenz: nu +- ∆nu = ({} +- {}) nm".format(sc.c/lambd[np.argmax(I)], sc.c*0.5/lambd[np.argmax(I)]**2))

"""
Hemisphärischer Resonator: Resonatorlänge
"""
p1 = np.array([1.23, 1.2, 1.17, 0.877, 128e-3, 1.32])  # in W
l1 = np.array([63.3, 69.8, 89.3, 95.6, 100.5, 79.5])  # in cm
sigma1 = np.array([2.8, 3.2, 3, 2.6, 1.8, 4.8]) * 1e-3
# niceplot(x=l1, y=p1, yerr=sigma1, c='tab:blue', ls='', marker='x',
#         titel='Einfluss der Resonatorlänge $L$', xaxis='L / cm', yaxis='P / mW',
#         error=True, size=(10, 5), safefig=True, safename='resonatorlength')

"""
Hemisphärischer Resonator: Ausgangsleistung mit Polarisator (Malus law)
"""
# Hat Lukas gemacht

"""
Kommerzieller HeNe-Laser Modenbetrachtung
"""
# Daten einlesen
t1, I1 = np.loadtxt('mitPol.txt', usecols=(0, 1), unpack=True, skiprows=4)
t2, I2 = np.loadtxt('ohnePol.txt', usecols=(0, 1), unpack=True, skiprows=4)
# get peaks for gaussian fit
indices1, parameters1 = find_peaks(I1, height=50)
indices2, parameters2 = find_peaks(I2, height=50)
# get noise
noise1 = noisecalc(I1, 3)
noise2 = noisecalc(I2, 8)
niceplot(x=t1[indices1], y=I1[indices1] - noise1,
         x2=t1, y2=I1 - noise1, c='tab:red', c2='tab:blue',
         plot2=True, ls='', marker='x', plotlabel='peaks', plotlabel2='Polarisationswinkel $\Theta=\pi/2$',
         titel='Modenbetrachtung des kommerziellen HeNe-Lasers', legend=True,
         xaxis='time', yaxis='intensity', size=(10,5),
         plot3=True, x3=t2, y3=I2 - noise2, c3='tab:green', plotlabel3='Polarisationswinkel $\Theta=2\pi$',
         plot4=True, x4=t2[indices2], y4=I2[indices2] - noise2, ls4='', marker4='x', c4='tab:red')
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
indices5, parameters5 = find_peaks(I5, height=40, distance=16)
# estimate noise
noise4, noise5 = noisecalc(I4, 30), noisecalc(I5, 30)
# Gauss-fits
# popt3, cov3 = fit(gaussian, t3[indices3], I3[indices3] # schlechte daten
                 # )
popt41, cov41 = fit(gaussian, t4[indices4[0:4]], I4[indices4[0:4]] - noise4, p0=[150, 139, 100],
                  #  bounds=[[50, 100, 10, 1],[1000, 200, 500, 50]]
                 )
popt42, cov42 = fit(gaussian, t4[indices4[4:8]], I4[indices4[4:8]] - noise4, p0=[150, 434, 100],
                 #   bounds=[[50, 400, 10, 1],[1000, 500, 500, 50]]
                 )
popt51, cov51 = fit(gaussian, t5[indices5[0:4]], I5[indices5[0:4]] - noise5, p0=[150, 133, 100],
                  #  bounds=[[50, 100, 10, 1],[1000, 200, 500, 50]]
                 )
popt52, cov52 = fit(gaussian, t5[indices5[4:8]], I5[indices5[4:8]] - noise5, p0=[150, 454, 100],
                   # bounds=[[100, 400, 10, 1],[1000, 500, 500, 50]]
                 )
print("GAUSSIAN PARAMETERS41:")
print("A+-∆A={}+-{}".format(popt41[0], np.sqrt(np.diag(cov41))[0]))
print("mu+-∆mu={}+-{}".format(popt41[1], np.sqrt(np.diag(cov41))[1]))
print("sigma+-∆sigma={}+-{}".format(popt41[2], np.sqrt(np.diag(cov41))[2]))
print("The gaussian-line-width is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt41[2],
                                                               2 * np.sqrt(2 * np.log(2)) *
                                                               np.sqrt(np.diag(cov41))[2]))
print("GAUSSIAN PARAMETERS42:")
print("A+-∆A={}+-{}".format(popt42[0], np.sqrt(np.diag(cov42))[0]))
print("mu+-∆mu={}+-{}".format(popt42[1], np.sqrt(np.diag(cov42))[1]))
print("sigma+-∆sigma={}+-{}".format(popt42[2], np.sqrt(np.diag(cov42))[2]))
print("The gaussian-line-width is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt42[2],
                                                               2 * np.sqrt(2 * np.log(2))  *
                                                               np.sqrt(np.diag(cov42))[2]))

print("GAUSSIAN PARAMETERS51:")
print("A+-∆A={}+-{}".format(popt51[0], np.sqrt(np.diag(cov51))[0]))
print("mu+-∆mu={}+-{}".format(popt51[1], np.sqrt(np.diag(cov51))[1]))
print("sigma+-∆sigma={}+-{}".format(popt51[2], np.sqrt(np.diag(cov51))[2]))
print("The gaussian-line-width is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt51[2],
                                                               2 * np.sqrt(2 * np.log(2)) *
                                                               np.sqrt(np.diag(cov51))[2]))
print("GAUSSIAN PARAMETERS52:")
print("A+-∆A={}+-{}".format(popt52[0], np.sqrt(np.diag(cov52))[0]))
print("mu+-∆mu={}+-{}".format(popt52[1], np.sqrt(np.diag(cov52))[1]))
print("sigma+-∆sigma={}+-{}".format(popt52[2], np.sqrt(np.diag(cov52))[2]))
print("The gaussian-line-width is: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt52[2],
                                                               2 * np.sqrt(2 * np.log(2)) *
                                                               np.sqrt(np.diag(cov52))[2]))
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