"""
Funktionen für die Auswertung von Versuchen
"""

import numpy as np
import matplotlib.pyplot as plt

def niceplot(x, y, c, lw=None, lw2=None, lw3=None, lw4=None, lw5=None, lw6=None,
             lw7=None, ls=None, ls6=None, ls7=None,
             ls2=None, ls3=None, ls4=None, ls5=None, plot2=None, plot3=None,
             plot4=None, plot5=None, plot6=None, plot7=None, x2=None, c2=None,
             y2=None, x3=None,
             y3=None, c3=None, x4=None, y4=None, x5=None, y5=None, c4=None,
             x6=None, y6=None, c6=None, x7=None, y7=None, c7=None,
             c5=None, marker=None, marker2=None, marker3=None, marker4=None,
             marker5=None, marker6=None, marker7=None, ylabel=None, xlabel=None,
             ms=10, cs=5, fs=15, ticks=6,
             size=(8, 8), safefig=None, error=None, errorlabel=None,
             safename=None, yaxis=None, xaxis=None, yerr=None, plotlabel=None,
             legend=None, plotlabel2=None, plotlabel3=None, plotlabel4=None,
             plotlabel5=None, plotlabel6=None, plotlabel7=None, titel=None,
             xlim=None,
             # ylim=None
             ):
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel(r'{}'.format(yaxis), size=fs)
    ax.set_xlabel(r'{}'.format(xaxis), size=fs)
    ax.set_title(r'{}'.format(titel), size=fs+3)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    #####
    # pinmp_ticks(ax.xaxis, ticks) # For special plotstyle
    # pinmp_ticks(ax.yaxis, ticks)
    # ax.grid(which='minor', alpha=.5)
    # ax.grid(which='major', alpha=.5)
    # ax.tick_params(which='both')
    #####
    # ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if error == True:
        ax.errorbar(x, y, yerr=yerr, capsize=cs, c=c, ms=ms, ls=ls,
                    marker=marker, label=r'{}'.format(plotlabel))
    else:
        ax.plot(x, y, ms=ms, lw=lw, ls=ls, marker=marker, c=c,
                label=r'{}'.format(plotlabel))
    # Bei Bedarf mehrere Fehlerbalken
    # if error2 == True:
    #     ax.errorbar(x, y, yerr=yerr2, capsize=cs, c=c2,
    #                  label = r'${}$'.format(errorlabel2))
    # if error3 == True:
    #     ax.errorbar(x, y, yerr=yerr3, capsize=cs, c=c3,
    #                  label = r'${}$'.format(errorlabel3))
    # if error4 == True:
    #     ax.errorbar(x, y, yerr=yerr4, capsize=cs, c=c4,
    #                  label = r'${}$'.format(errorlabel4))
    if plot2 == True:
        ax.plot(x2, y2, ms=ms, lw=lw2, ls=ls2, marker=marker2, c=c2,
                label=r'{}'.format(plotlabel2))
    if plot3 == True:
        ax.plot(x3, y3, ms=ms, lw=lw3, ls=ls3, marker=marker3, c=c3,
                label=r'{}'.format(plotlabel3)
                 )
    if plot4 == True:
        ax.plot(x4, y4, ms=ms, lw=lw4, ls=ls4, marker=marker4, c=c4,
             #   label=r'{}'.format(plotlabel4)
                )
    if plot5 == True:
        ax.plot(x5, y5, ms=ms, lw=lw5, ls=ls5, marker=marker5, c=c5,
              #  label=r'{}'.format(plotlabel5)
        )
    if plot6 == True:
        ax.plot(x6, y6, ms=ms, lw=lw6, ls=ls6, marker=marker6, c=c6,
               # label=r'{}'.format(plotlabel6)
        )
    if plot7 == True:
        ax.plot(x7, y7, ms=ms, lw=lw7, ls=ls7, marker=marker7, c=c7,
            #    label=r'{}'.format(plotlabel7)
        )
    if legend == True:
        ax.legend(fontsize=fs-4, markerscale=ms/10, facecolor='white')
    if safefig == True:
        plt.savefig('{}.svg'.format(safename), dpi=300)

    plt.show()