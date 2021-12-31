import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy.ndimage import gaussian_filter1d

x, intensity = np.genfromtxt("02122021goldkrass.txt", skip_header=1, usecols=(0,1), unpack=True)
plt.fill_between(x, intensity, 350, alpha=0.7)
ysmoothed = gaussian_filter1d(intensity, sigma=1)
plt.plot(x, ysmoothed)
peaks = signal.find_peaks(ysmoothed, x)
for i in peaks[0]:
    plt.axvline(x[i])
plt.show()