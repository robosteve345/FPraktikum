import numpy as np
import glob
import matplotlib.pyplot as plt
from pathlib import Path
import re
import pandas as pd
from scipy.optimize import curve_fit
import functools
from matplotlib import ticker

# Plot stuff
plt.style.use('ggplot')
plt.rcParams.update({
    'font.family': 'serif',
    'text.usetex': False,
    'pgf.rcfonts': False,
})

def kaustik(x, b, f, s, lam):
    A = 1 - x/f
    B = s + x - x*s/f
    C = -1/f
    D = 1 - s/f
    b_p = b*(A*D - C*B)/(A**2 + b**2*B**2)
    return np.sqrt(lam/(np.pi*b_p))

def malus(x, a, b):
    return a*np.cos(x + b)**2

def read_files(path):
    names = glob.glob(f"{Path().resolve()}\{path}\*.csv")
    data = np.array([])
    for i in names:
        x_data = int(re.search(r'\d+', Path(i).stem).group())
        mess = pd.read_csv(i, encoding_errors='ignore', delimiter=';')
        mess = mess.to_numpy()
        mess_ind = np.indices(np.shape(mess))

        y_index1_r = mess_ind[0][mess == 'FWHM vert']
        y_index1_c = mess_ind[1][mess == 'FWHM vert'] + 1
        y_data1 = mess[y_index1_r, y_index1_c]
        y_data1 = y_data1[~np.isnan(y_data1.astype(float))]
        y_index2_r = mess_ind[0][mess == 'FWHM horiz']
        y_index2_c = mess_ind[1][mess == 'FWHM horiz'] + 1
        y_data2 = mess[y_index2_r, y_index2_c]
        y_data = np.append(y_data1, y_data2[~np.isnan(y_data2.astype(float))])

        y_data = np.mean(y_data)/np.sqrt(2*np.log(2))*5.6e-3
        data = np.append(data, [x_data, y_data])

    data = np.reshape(data, (len(names), 2))
    x_data = data[:, 0]
    y_data = data[:, 1]
    return x_data, y_data

def pinmp_ticks(axis, ticks):
    axis.set_major_locator(ticker.MaxNLocator(ticks))
    axis.set_minor_locator(ticker.MaxNLocator(ticks * 10))
    return axis

def set_up_plot(ticks=10):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    pinmp_ticks(ax.xaxis, ticks)
    pinmp_ticks(ax.yaxis, ticks)

    ax.grid(which='minor', alpha=.3)
    ax.grid(which='major', alpha=.8)

    ax.tick_params(which='both')

    return fig, ax

def save_fig(fig, title, folder='unsorted', size=(5, 4)):
    fig.set_size_inches(*size)
    fig.tight_layout()
    fig.savefig(f'{title}.pdf')
    fig.savefig(f'{title}.svg')

def main():
    x_data, y_data = read_files("Kaustik")
    marker = 'x'
    dist = 750
    lam = 632.8e-6
    print(lam)
    kaustik_lam = functools.partial(kaustik, lam=lam)
    kaustik_free = functools.partial(kaustik, lam=lam, x=0, f=np.inf)
    params = curve_fit(kaustik_lam, x_data, y_data, p0=[0.002, 150, 750])[0]
    print(params)
    dist = params[-1]
    plot_x = np.linspace(0, max(x_data), 100)
    fig1, ax1 = set_up_plot()

    ax1.errorbar(x_data + dist, y_data, ls='', marker=marker, xerr=5, label='Experimentelle Daten', c='tab:blue')
    ax1.plot(plot_x+params[-1], kaustik_lam(plot_x, *params), color='tab:red', label='Gefitteter Theoretischer Verlauf')
    ax1.plot(np.linspace(0, params[-1]), kaustik_free(s=np.linspace(0, params[-1]), b=params[0]), color='tab:red')
    ax1.set_ylabel(r'$\omega^\prime$/mm')
    ax1.set_xlabel(r'Abstand zum Endspiegel im Resonator/mm')
    ax1.legend()
    ax1.set_title('Kaustik des Laserstrahls mit Linse')
    save_fig(fig1, 'kaustik', size=(7, 5))

    phi = np.array([96, 106, 116, 86, 76, 56, 36, 26, 16, 6, -4, -14])/360*2*np.pi
    p2 = np.array(
        [2.464e-6, 30.9e-6, 110e-6, 25.3e-6, 102e-6, 0.399e-3, 0.744e-3, 0.858e-3, 0.972e-3, 1.02e-3, 0.994e-3,
         0.907e-3])*1000
    sigma2 = np.array([24e-9, 0.117e-6, 0.6e-6, 0.19e-6, 6.54e-6, 1.1e-3, 3.7e-3, 3.7e-3, 3 - 2e-3, 2.4e-3, 2.3e-3, 2.0e-3])
    fig2, ax2 = set_up_plot()

    params2 = curve_fit(malus, phi, p2)[0]#, sigma=sigma2, absolute_sigma=True)[0]
    ax2.errorbar(phi, p2, yerr=0, ls='', marker=marker, label='Experimentelle Daten', c='tab:blue')
    plot_phi = np.linspace(min(phi), max(phi))
    print(params2)
    ax2.plot(plot_phi, malus(plot_phi, *params2), color='tab:red', label=r'Fit an Malus Gesetz: $I=I_0\cos(\phi)$')
    ax2.set_xlabel(r'$\phi$/rad')
    ax2.set_ylabel(r'$I$/mW')
    ax2.set_title('Malus Gesetz anhand des Laserlichtes')
    ax2.legend()
    save_fig(fig2, 'malus', size=(7,5))
    plt.show()

if __name__ == "__main__":
    main()