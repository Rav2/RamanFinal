#########################################################################
####################### IMPORTING REQUIRED MODULES ######################

import numpy
import pylab
from scipy.optimize import leastsq  # Levenberg-Marquadt Algorithm #


#########################################################################
########################### DEFINING FUNCTIONS ##########################

def lorentzian(x, coefficients):
    y = 0
    for p in coefficients:
        numerator = (p[0] ** 2)
        denominator = (x - (p[1])) ** 2 + p[0] ** 2
        y += p[2] * (numerator / denominator)
    return y


def residuals(p, y, x):
    elements_number = len(p)
    coefficients = []
    if elements_number % 3 == 0:
        for i in range(0, elements_number / 3):
            coefficients.append(p[3 * i: 3 * (i + 1)])
    else:
        raise Exception('Number of coefficients provided for residuals(p,y,x) is wrong!!!')
    err = y - lorentzian(x, coefficients)
    return err


def perform_fitting(X, Y):
    ind_bg_low = (X > min(X)) & (X < 450.0)
    ind_bg_high = (X > 600.0) & (X < max(X))
    x_bg = numpy.concatenate((X[ind_bg_low], X[ind_bg_high]))
    y_bg = numpy.concatenate((Y[ind_bg_low], Y[ind_bg_high]))
    # remove background
    m, c = numpy.polyfit(x_bg, y_bg, 1)
    background = m * X + c
    y_bg_corr = Y - background
    y_bg_peak_corr = Y - background
    # fit triple lorentz
    coefficients = [[3.5, 521.6, 12e3], [70, 480, 6000], [25, 500,
                                                          14000]]  # [[3.5,521.6,12e3],[60,480,30],[30,496,3000]]  # [hwhm, peak center, intensity] for three peaks
    pbest = leastsq(residuals, coefficients, args=(y_bg_corr, X), full_output=1)
    best_parameters = [pbest[0][:3], pbest[0][3:6], pbest[0][6:]]
    # subtract peaks at 480cm^-1 and 496cm^-1
    best_params = iter(best_parameters)
    next(best_params)
    for triple in best_params:
        y_bg_peak_corr -= (triple[0] ** 2 / ((X - triple[1]) ** 2 + triple[0] ** 2) * triple[2])
    return y_bg_peak_corr


#########################################################################
def main():
    #########################################################################
    ############################# LOADING DATA ##############################    
    a = numpy.loadtxt("../content/semiamorphSi1-514-20s-02.txt")
    x = a[:, 0]
    y = a[:, 1]
    #########################################################################
    ######################## BACKGROUND SUBTRACTION #########################    
    # defining the 'background' part of the spectrum #
    ind_bg_low = (x > min(x)) & (x < 450.0)
    ind_bg_high = (x > 600.0) & (x < max(x))

    x_bg = numpy.concatenate((x[ind_bg_low], x[ind_bg_high]))
    y_bg = numpy.concatenate((y[ind_bg_low], y[ind_bg_high]))

    # fitting the background to a line # 
    m, c = numpy.polyfit(x_bg, y_bg, 1)
    # m,c-polynomial coefficients
    #########################################################################
    # removing fitted background # 
    background = m * x + c
    y_bg_corr = y - background
    y_bg_peak_corr = y - background
    #########################################################################
    ############################# FITTING DATA ##############################

    # dla malych grubosci ziaren(<10nm) widac znaczace roznice, dla innych- widma sa takie same
    # initial values #
    coefficients = [[3.5, 521.6, 12e3], [60, 480, 5000],
                    [32, 500, 8000]]  # [hwhm, peak center, intensity] for three peaks
    # print (coefficients)
    # optimization # 
    pbest = leastsq(residuals, coefficients, args=(y_bg_corr, x), full_output=1)
    best_parameters = [pbest[0][:3], pbest[0][3:6], pbest[0][6:]]

    # fit to data #
    fit = lorentzian(x, best_parameters)
    #########################################################################
    ############################## REMOVING 480 AND 496 PEAKS #################################
    best_params = iter(best_parameters)
    next(best_params)
    for triple in best_params:
        y_bg_peak_corr -= (triple[0] ** 2 / ((x - triple[1]) ** 2 + triple[0] ** 2) * triple[2])
    #new lorentzian to cut function
    # initial values #
    coefficients1 = [3.5, 521.6, 12e3]  # [hwhm, peak center, intensity] #
    pbest1 = leastsq(residuals, coefficients1, args=(y_bg_peak_corr, x), full_output=1)
    best_parameters1 = pbest1[0]
    fit1 = lorentzian(x, [best_parameters1])  #changed from:     fit1 = lorentzian(x,best_parameters)
    #########################################################################
    ############################## PLOTTING #################################  
    pylab.plot(x, y_bg_corr, 'o', label="data without bcg")
    pylab.plot(x, y_bg_peak_corr, 'o', label="data with only main peak")
    pylab.plot(x, fit, '-', lw=2, label="lorentz1")
    pylab.plot(x, fit1, '-', label="lorentz2")
    #pylab.plot(x,y,"ro", label="original data")
    pylab.legend(loc=1)

    pylab.xlabel(r'$\omega$ (cm$^{-1}$)', fontsize=18)
    pylab.ylabel('Intensity (a.u.)', fontsize=18)

    pylab.show()
    #########################################################################
    ############################## READING FIT PARAMETERS FROM 2ND LORENTZIAN #################################    
    height = best_parameters1[2]
    print "height", height
    peak_center = best_parameters1[1]
    print "peak_center", peak_center
    hwhm = best_parameters1[0]
    print "hwhm", hwhm

#########################################################################
if __name__ == "__main__":
    main()