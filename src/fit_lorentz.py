#########################################################################
####################### IMPORTING REQUIRED MODULES ######################

import numpy
import pylab
from scipy.optimize import leastsq  # Levenberg-Marquadt Algorithm #




#########################################################################
########################### DEFINING FUNCTIONS ##########################

def lorentzian(x, coefficients):
    y = 0
    if len(coefficients)==2:
        # within_bounds = (coefficients[0][1]>515)*(coefficients[1][1]<490)*(515>coefficients[2][1]>480) * \
        within_bounds =(coefficients[1][2]>0)*(coefficients[1][1]>490)#*(475<coefficients[1][1]<485)*(505>coefficients[2][1]>490)*(525>coefficients[0][1]>518)
    elif len(coefficients)==3:
        within_bounds =(coefficients[2][2]>0)*(coefficients[1][2]>0)*(475<coefficients[2][1]<485)*(500>coefficients[1][1]>494)
    else:
        within_bounds=True
    for p in coefficients:
        numerator = (p[0] ** 2)
        denominator = (x - (p[1])) ** 2 + p[0] ** 2
        y += p[2] * (numerator / denominator)+coefficients[0][3]
    if within_bounds:
        return y
    else:
        return -10e10


def residuals(p, y, x):
    elements_number = len(p)
    coefficients = []
    if elements_number % 4 == 0:
        for i in range(0, elements_number / 4):
            coefficients.append(p[4 * i: 4 * (i + 1)])
    else:
        raise Exception('Number of coefficients provided for residuals(p,y,x) is wrong!!!')
    err = y - lorentzian(x, coefficients)

    return err


def perform_fitting(X, Y, omega_0, gamma_0, amplitude, sample_type):
    ind_bg_low = (X > min(X)) & (X < 400.0)
    ind_bg_high = (X > 600.0) & (X < max(X))
    x_bg = numpy.concatenate((X[ind_bg_low], X[ind_bg_high]))
    y_bg = numpy.concatenate((Y[ind_bg_low], Y[ind_bg_high]))
    # remove background
    m, c = numpy.polyfit(x_bg, y_bg, 1)
    background = m * X + c
    y_bg_corr = Y - background
    y_bg_peak_corr = Y - background
    # fit triple lorentz
    if sample_type==2:
        coefficients = [[3.5, 521.6, max(y_bg_corr), 0],
                        [19, 497, 1000, 0]]  # [hwhm, peak center, intensity] for three peaks
    elif sample_type==3:
        coefficients = [[3.5, 521.6, max(y_bg_corr), 0],[15, 497, 1000, 0], [15, 480, 1200, 0]]
    fitting_range=(X>400)*(X<800)
    pbest = leastsq(residuals, coefficients, args=(y_bg_corr[fitting_range], X[fitting_range]), full_output=1)
    #print pbest[0]
    best_parameters = []
    for i in range(0, int(len(pbest[0])/4)):
        best_parameters.append(pbest[0][4*i : 4*(i+1)])
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
    sample_type=3
    a = numpy.loadtxt("../content/SingleFiles/10.0.txt")
    x = a[:, 2]
    y = a[:, 3]
    #########################################################################
    ######################## BACKGROUND SUBTRACTION #########################
    # defining the 'background' part of the spectrum #
    ind_bg_low = (x > min(x)) & (x < 400.0)
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
    if sample_type==2:
        coefficients = [[3.5, 521.6, max(y_bg_corr), 0],# [70, 480, 100, 0],
                        [19, 497, 1000, 0]]  # [hwhm, peak center, intensity] for three peaks
    elif sample_type==3:
        coefficients = [[3.5, 521.6, max(y_bg_corr), 0],[15, 497, 1000, 0], [15, 480, 1200, 0]]
    # print (coefficients)
    # optimization #
    fitting_range=(x>400)*(x<800)
    pbest = leastsq(residuals, coefficients, args=(y_bg_corr[fitting_range], x[fitting_range]), full_output=1)
    print pbest[0]
    best_parameters = []
    for i in range(0, int(len(pbest[0])/4)):
        best_parameters.append(pbest[0][4*i : 4*(i+1)])
    best_parameters521=[pbest[0][:4]]
    best_parameters497=[pbest[0][4:8]]
    #best_parameters480=[pbest[0][8:]]
    # fit to data #
    fit = lorentzian(x, best_parameters)
    fit2 = lorentzian(x, best_parameters497)
    #fit3 = lorentzian(x, best_parameters480)
    #########################################################################
    ############################## REMOVING 480 AND 496 PEAKS #################################
    best_params = iter(best_parameters)
    next(best_params)
    for triple in best_params:
        y_bg_peak_corr -= (triple[0] ** 2 / ((x - triple[1]) ** 2 + triple[0] ** 2) * triple[2])
    # for i in range(0, len(y_bg_peak_corr)):
    #     if y_bg_peak_corr[i] < 0:
    #         y_bg_peak_corr[i]=0
    #new lorentzian to cut function
    # initial values #
    coefficients1 = [3.5, 521.6, 12e3, 0]  # [hwhm, peak center, intensity] #
    pbest1 = leastsq(residuals, coefficients1, args=(y_bg_peak_corr, x), full_output=1)
    best_parameters1 = pbest1[0]
    fit1 = lorentzian(x, [best_parameters1])  #changed from:     fit1 = lorentzian(x,best_parameters)
    #########################################################################
    ############################## PLOTTING #################################
    pylab.plot(x, y_bg_corr, 'o', label="data without bcg")
    pylab.plot(x, y_bg_peak_corr, 'o', label="data with only main peak")
    pylab.plot(x, fit, '-', lw=2, label=(""))
    pylab.plot(x, fit2, '-', lw=2, label=(""))
    #pylab.plot(x, fit3, '-', lw=2, label=(""))
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
    yTest=perform_fitting(x, y, 0, 0, 0, 2)
    pylab.plot(x, yTest, '.')
    pylab.show()
#########################################################################
if __name__ == "__main__":
    main()