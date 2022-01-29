import numpy as np
import glob
import matplotlib.pyplot as plt
from pathlib import Path
import re
import pandas as pd
from scipy.optimize import curve_fit
import functools
from matplotlib import ticker
import matplotlib

# Plot stuff
plt.style.use('ggplot')
plt.rcParams.update({
    'font.family': 'serif',
    'text.usetex': False,
    'pgf.rcfonts': False,
})

def kaustik(x, b, s, f, lam):
    A = 1 - x/f
    B = s + x - x*s/f
    b_p = b*(1)/(A**2 + b**2*B**2)
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
        # Daten finden
        y_index1_r = mess_ind[0][mess == 'FWHM vert']
        y_index1_c = mess_ind[1][mess == 'FWHM vert'] + 1
        y_data1 = mess[y_index1_r, y_index1_c][~np.isnan(mess[y_index1_r, y_index1_c].astype(float))]
        y_index2_r = mess_ind[0][mess == 'FWHM horiz']
        y_index2_c = mess_ind[1][mess == 'FWHM horiz'] + 1
        y_data2 = mess[y_index2_r, y_index2_c][~np.isnan(mess[y_index2_r, y_index2_c].astype(float))]
        # Durchschnitt bilden
        y_data = np.mean([y_data1, y_data2])/np.sqrt(2*np.log(2))*5.6e-3
        s_y = np.abs(y_data1 - y_data2)[0]/2/np.sqrt(2*np.log(2))*5.6e-3
        data = np.append(data, [x_data, y_data, s_y])

    data = np.reshape(data, (len(names), 3))
    x_data = data[:, 0]
    y_data = data[:, 1]
    s_y = data[:, 2]
    return x_data, y_data, s_y

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

def center_axis(ax, fig):
    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    # Show ticks in the left and lower axes only
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    def offset(x, y):
        return matplotlib.transforms.ScaledTranslation(x/72., y/72., fig.dpi_scale_trans)
    for label in ax.xaxis.get_majorticklabels()[1:5]:
        label.set_transform(label.get_transform() + offset(0, 16))
    for label in ax.yaxis.get_majorticklabels()[1:]:
        label.set_transform(label.get_transform() + offset(22, 0.5))
    ax.yaxis.get_majorticklabels()[5].set_transform(label.get_transform() + offset(-17, 7))
    ax.xaxis.get_majorticklabels()[5].set_transform(label.get_transform() + offset(7, -10))



def hyperbola(x, ax):
    ax.plot(x, 1 / x, c='tab:blue')
    ax.fill_between(x, 0, 1/x, alpha=0.6, color='royalblue')

def save_fig(fig, title, folder='unsorted', size=(5, 4)):
    fig.set_size_inches(*size)
    fig.tight_layout()
    fig.savefig(f'{title}.pdf')
    fig.savefig(f'{title}.svg')

def main():
    x_data, y_data, s_y = read_files("Kaustik")
    marker = 'x'
    dist = 750
    lam = 632.8e-6
    print(lam)
    kaustik_lam = functools.partial(kaustik, lam=lam, f=150) # Problem: distance sehr stark korreliert mit f (und b?)2
    kaustik_free = functools.partial(kaustik, lam=lam, x=0, f=np.inf)
    params, pcov = curve_fit(kaustik_lam, x_data, y_data, p0=[0.002, 150], sigma=s_y, absolute_sigma=True)
    print(params)
    print(np.sqrt(np.diag(pcov)))
    #dist = params[-1]
    plot_x = np.linspace(0, max(x_data), 100)
    fig1, ax1 = set_up_plot()

    ax1.errorbar(x_data + dist, y_data, ls='', marker=marker, xerr=5, yerr=s_y, label='Experimentelle Daten', c='tab:blue')
    ax1.plot(plot_x+dist, kaustik_lam(plot_x, *params), color='tab:red', label='Gefitteter Theoretischer Verlauf')
    ax1.plot(np.linspace(0, dist), kaustik_free(s=np.linspace(0, dist), b=params[0]), color='tab:red')
    ax1.set_ylabel(r'$\omega^\prime$/mm')
    ax1.set_xlabel(r'Abstand zum Endspiegel im Resonator/mm')
    ax1.legend()
    ax1.set_title('Kaustik des Laserstrahls mit Linse')
    save_fig(fig1, 'kaustik', size=(7, 5))

    phi = np.array([96, 106, 116, 86, 76, 56, 36, 26, 16, 6, -4, -14])/360*2*np.pi
    p2 = np.array(
        [2.464e-6, 30.9e-6, 110e-6, 25.3e-6, 102e-6, 0.399e-3, 0.744e-3, 0.858e-3, 0.972e-3, 1.02e-3, 0.994e-3,
         0.907e-3])*1000
    sigma2 = np.array([24e-9, 0.117e-6, 0.6e-6, 0.19e-6, 6.54e-6, 1.1e-6, 3.7e-6, 3.7e-6, 3.2e-6, 2.4e-6, 2.3e-6, 2.0e-6])*1000
    fig2, ax2 = set_up_plot()

    params2 = curve_fit(malus, phi, p2, sigma=sigma2, absolute_sigma=True)[0]
    ax2.errorbar(phi, p2, yerr=sigma2,xerr=2/360*2*np.pi, ls='',marker=marker, label='Experimentelle Daten', c='tab:blue')
    plot_phi = np.linspace(min(phi), max(phi))
    print(params2)
    ax2.plot(plot_phi, malus(plot_phi, *params2), color='tab:red', label=r'Fit an Malus Gesetz: $I=I_0\cos(\phi)$')
    ax2.set_xlabel(r'$\phi$/rad')
    ax2.set_ylabel(r'$I$/mW')
    ax2.set_title('Malus Gesetz anhand des Laserlichtes')
    ax2.legend()
    save_fig(fig2, 'malus', size=(7,5))

    plt.rcParams.update({'figure.autolayout': True})
    fig3 = plt.figure()
    ax3 = fig3.subplots()
    stab_plot_x_r = np.linspace(0, 4.5, 200)
    center_axis(ax3, fig3)
    hyperbola(stab_plot_x_r, ax3)
    hyperbola(-stab_plot_x_r, ax3)
    ax3.vlines(1, 0, 1, color='tab:red', label='Arbeitsbereich')
    ax3.set_ylabel(r"$g_2$", labelpad=0, loc='top')
    ax3.set_xlabel(r"$g_1$", labelpad=0, loc='right')
    ax3.set_ylim(-4.5, 4.5)
    ax3.set_xlim(-4.5, 4.5)
    ax3.legend()
    save_fig(fig3, "stab", size=(5,5))
    plt.show()

if __name__ == "__main__":
    main()