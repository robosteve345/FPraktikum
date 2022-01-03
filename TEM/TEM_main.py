import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy.ndimage import gaussian_filter1d
from scipy import interpolate
import os
import pandas as pd


def weird_even(x):
    if x % 2 == 0:
        return 1
    else:
        return -1


def smoothing(x, y, steps=1000, sigma=6):
    x_new = np.linspace(min(x), max(x), steps)  # make stuff more finely grained
    y_new = interpolate.interp1d(x, y)(x_new)  # by linear interpolation
    ysmoothed = gaussian_filter1d(y_new, sigma=sigma)  # make graph smooth by gaussian filter
    print(sigma)
    return x_new, ysmoothed


def real_image(x, intensity, start, ax):
    x_new, intensity_smoothed = smoothing(x, intensity)
    peaks = x_new[signal.find_peaks(intensity_smoothed)[0]]
    d = peaks[start] - peaks[start + 1]
    for i in np.arange(start + 1, len(peaks) - 3):
        d_new = peaks[i + 1] - peaks[i]
        if d - d_new < 0.25 * d:
            d = (peaks[i] - peaks[start]) / (i - start)  # new distance is average of old distances
        else:
            i_final = i - 1
            break
        i_final = i - 1
    if ax:
        ax.axvspan(peaks[start], peaks[i_final], color='red', label='Avaraged area')

    d = (peaks[i_final] - peaks[start]) / (i_final - start)
    print([i_final])
    return intensity_smoothed, peaks, d,


def diffraction(x_diff, intensity_diff):
    x_new, ysmoothed = x_diff, intensity_diff # smoothing(x_diff, intensity_diff, sigma=5.5, steps=len(intensity_diff))  # make stuff more finely grained
    peaks = signal.find_peaks(ysmoothed, width=5, threshold=0.001*ysmoothed)  # find peaks of smooth stuff
    a_gold = 4.078 * 1e-1  # Gitterkonstante in nm
    d = {}
    for h in range(5):
        for k in range(5):
            for l in range(5):
                if l == k == h == 0:
                    pass
                elif (weird_even(h + k) + weird_even(h + l) + weird_even(l + k)) > 0:
                    d[(h, k, l)] = a_gold / (h ** 2 + k ** 2 + l ** 2) ** 0.5
    d_experimental = (1 / x_new[peaks[0]])
    return ysmoothed, x_new[peaks[0]], d_experimental, d


def plot(ax, x, intensity, intensity_smoothed, peaks, real_image=True):
    # plot stuff
    x_new = np.linspace(min(x), max(x), len(intensity_smoothed))
    ax.set_ylim(min(intensity) - 2, max(intensity) + 1)
    ax.set_xlim(0, max(x))
    if real_image:
        ax.set_xlabel(r"$d$[nm]")
    else:
        ax.set_xlabel(r"$r$[1/nm]")
    ax.set_ylabel("Intensity")
    ax.fill_between(x, intensity, ax.get_ylim()[1], alpha=0.7)
    ax.plot(x, intensity, ls='', marker='x', ms=1, label="Measured data", color='blue')
    #ax.plot(x_new, intensity_smoothed, ls='--', marker='', label='Smoothed data', color='red')
    ax.vlines(peaks, ax.get_ylim()[0], ax.get_ylim()[1], label='peaks')


def save_fig(fig, title, folder='unsorted', size=(9, 7)):
    fig.set_size_inches(*size)
    fig.tight_layout()
    try:
        os.makedirs(f'./figs/{folder}/')
    except OSError as exc:
        pass
    fig.savefig(f'./figs/{folder}/{title}.pdf')


def isfull(x):
    if pd.isna(x):
        return False
    else:
        return True


def main():
    intensity_diff = np.genfromtxt("BeugungProfile Of Au_difrf3Pol.txt", skip_header=0, unpack=False)  # read stuff
    """
    data = pd.read_excel("datahighresolution.xlsx", usecols=isfull)
    data.dropna(how='all', axis=1, inplace=True)  # remove non-filled columns between useful ones
    datas = {}
    d_list = []
    fig_real = plt.figure()
    j = 0
    starters = [3, 6, 3, 2, 3]
    for i in [2, 3, 6]:
        rough = data.iloc[:, [2 * i, 2 * i + 1]].copy()  # take pairs of columns, copy them
        rough.dropna(how='any', axis=0, inplace=True)  # remove all non filled rows (headers)
        datas[i] = ([rough.to_numpy(dtype=np.dtype(np.float32))[:, 0],
                     rough.to_numpy(dtype=np.dtype(np.float32))[:, 1]])  # turn it into numpy array
        if j < 3:
            ax_real = fig_real.add_subplot(1, 3, j + 1)
            intensity_real, peaks_real, d = real_image(*datas[i], starters[j], ax_real)
            plot(ax_real, *datas[i], intensity_real, peaks_real, j)
        else:
            intensity_real, peaks_real, d = real_image(*datas[i], starters[j], None)
        d_list.append(d)

        j = j + 1

    d_list = np.array(d_list)
    d_avg = np.average(d_list)
    d_sigma = np.std(d_list, ddof=1)
    print(d_list)
    df = pd.DataFrame(d_list)
    ax_real.legend(loc=1)
    save_fig(fig_real, 'analysis_real', size=(14, 5))
    """
    maximum_diff = 17.4875
    x_diff = np.linspace(0, maximum_diff, len(intensity_diff))

    intensity_smooth, peaks_diff, d_exp, d_theo = diffraction(x_diff, intensity_diff)
    print(peaks_diff)
    fig_diff = plt.figure()
    ax_diff = fig_diff.subplots()
    plot(ax_diff, x_diff, intensity_diff, intensity_smooth, peaks_diff, real_image=False)
    miller_doubles = []
    for i in d_theo:
        millersum = np.sum(np.array([*i]) ** 2)
        if 1 / d_theo[i] < ax_diff.get_xlim()[1] and millersum not in miller_doubles:
            ax_diff.text(1 / d_theo[i] + 0.1, 40 - np.sqrt(millersum), str([i[0], i[1], i[2]]), rotation='vertical', color='red')
            ax_diff.axvline(1 / d_theo[i], color='red')
            miller_doubles.append(millersum)
    # ax_diff.vlines(1/d_theo, ymin=0, ymax=max(intensity_diff)+1, color='red', label='theoretical values')
    ax_diff.legend()
    save_fig(fig_diff, "differential")
    plt.show()


if __name__ == "__main__":
    main()
