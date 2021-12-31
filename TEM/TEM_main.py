import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy.ndimage import gaussian_filter1d
from scipy import interpolate

x, intensity = np.genfromtxt("02122021goldkrass.txt", skip_header=1, usecols=(0,1), unpack=True) # read stuff

x_new = np.linspace(min(x), max(x), 200)                    # make stuff more finely grained
intensity_new = interpolate.interp1d(x, intensity)(x_new)   # by linear interpolation
ysmoothed = gaussian_filter1d(intensity_new, sigma=1)       # make graph smooth by gaussian filter
peaks = signal.find_peaks(ysmoothed, x_new)                 # find peaks of smooth stuff

# plot stuff
plt.fill_between(x, intensity, 350, alpha=0.7)
plt.plot(x_new, ysmoothed)
for i in peaks[0]:
    plt.axvline(x_new[i])
plt.show()