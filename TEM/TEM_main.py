

#++üüüüüüüüüüüüüüüüüüüüüüüüüüüüüüüü+#
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy.ndimage import gaussian_filter1d
from scipy import interpolate

def weird_even(x):
    if x%2 == 0:
        return 1
    else:
        return -1

intensity = np.genfromtxt("BeugungProfile Of Au_difrf3Pol.txt", skip_header=0, unpack=False)  # read stuff
maximum = 16.3
a_gold = 4.078 * 1e-1  # Gitterkonstante in nm
x = np.linspace(0, maximum, len(intensity))
x_new = np.linspace(min(x), max(x), 1000)  # make stuff more finely grained
intensity_new = interpolate.interp1d(x, intensity)(x_new)  # by linear interpolation
ysmoothed = gaussian_filter1d(intensity_new, sigma=9)  # make graph smooth by gaussian filter
peaks = signal.find_peaks(ysmoothed)  # find peaks of smooth stuff
d = np.zeros((5, 5, 5))
for h in range(5):
    for k in range(5):
        for l in range(5):
            if l == k == h == 0:
                pass
            elif (weird_even(h+k)+weird_even(h+l)+weird_even(l+k))>0:
                d[h, k, l] = a_gold / (h ** 2 + k ** 2 + l ** 2) ** 0.5
                #plt.axvline(1 / (a_gold / (h ** 2 + k ** 2 + l ** 2) ** 0.5), color ='red')
                plt.text(1 / (a_gold / (h ** 2 + k ** 2 + l ** 2) ** 0.5), h*10+k, (h,k,l))
d_experimental = (1/x_new[peaks[0]])
millers = []
for d_i in d_experimental:
    millers.append(np.argmin(np.abs(d_i-d), keepdims=True))
print(d_experimental)
# plot stuff
plt.ylim(min(intensity), max(intensity) + 1)
plt.xlim(0, max(x))
plt.xlabel(r"$r$[1/nm]")
plt.ylabel("Intensity")
plt.fill_between(x, intensity, 350, alpha=0.7)
plt.plot(x, intensity, ls='', marker='x', label="Measured data", color='blue')
plt.plot(x_new, ysmoothed, ls='--', marker='', label='Smoothed data', color='red')
for h in peaks[0][:-1]:
    plt.axvline(x_new[h])
plt.axvline(x_new[peaks[0][-1]], label='Peaks')
plt.legend()
plt.savefig('diffrac_plot.pdf')
plt.show()
