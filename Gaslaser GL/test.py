import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit

def gaussian(x, A, mu, sigma, B):
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

indices4=np.array([ 84, 114, 135, 164, 382, 409, 426, 457])
indices5=np.array([ 80, 116, 151, 185, 393, 435, 475, 500])
height4=np.array([ 46,  75, 119,  82,  41,  74, 113,  83])
height5=np.array([ 42, 133, 128,  43,  45, 128, 134,  43])

popt, cov=fit(gaussian, indices4[0:4], height4[0:4], absolute_sigma=True)
print("GAUSSIAN PARAMETERS:")
print("A+-∆A={}+-{}".format(popt[0], np.sqrt(np.diag(cov))[0]))
print("mu+-∆mu={}+-{}".format(popt[1], np.sqrt(np.diag(cov))[1]))
print("sigma+-∆sigma={}+-{}".format(popt[2], np.sqrt(np.diag(cov))[2]))
print("B+∆B=({}+-{})mV".format(popt[3], np.sqrt(np.diag(cov))[3]))

print("FWHM: ∆tau +-dtau = {}+-{}".format(2 * np.sqrt(2 * np.log(2)) * popt[2],
                                              2 * np.sqrt(2 * np.log(2)) *
                                              np.sqrt(np.diag(cov))[2]))

