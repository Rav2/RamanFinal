from __future__ import unicode_literals
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc  # used to customize all kinds of properties(size,dpi,color,style,etc.)
import os
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator  # do formatowania osi

# rc('font', family='Consolas')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator


def fit_function(x, loc_par, gamma, height, offset):
    return height * gamma ** 2 / ((x - loc_par) ** 2 + gamma ** 2) + offset


def perform_fitting(raman_shift, intensity, min_raman_shift, max_raman_shift):
    X_bool = (raman_shift < max_raman_shift) * (raman_shift > min_raman_shift)
    # print X_bool
    X_fit = raman_shift[X_bool]
    Y_fit = intensity[X_bool]
    #    print raman_shift
    #    print X_fit
    #    print intensity
    #    print Y_fit
    n = len(X_fit)

    mean = sum(X_fit) / n  #loc_par
    sigmal = math.sqrt(sum((X_fit - mean) ** 2) / n)

    popt, pcov = curve_fit(fit_function, X_fit, Y_fit, p0=[mean, sigmal, max(Y_fit), 0])
    #NAJWAZNIEJSZA LINIJKA
    loc_par_fit, gamma_fit, height_fit, offset1_fit = popt

    figure = plt.Figure()
    my_plot = figure.add_subplot(111)
    my_plot.plot(X_fit, Y_fit, "o", label="experimental data", color="red")
    my_plot.plot(X_fit, fit_function(X_fit, popt[0], popt[1], popt[2], popt[3]), linewidth="2.0", color="green",
                 label='fit')
    my_plot.legend()
    #plt.text(100, 20200, r'$\ x_0= %.3f,\u0263 = %.3f$' %(popt[0], math.fabs(popt[1])))
    #plt.savefig(name+".png", bbox_inches='tight')
    #plt.show()
    return (popt, figure)

