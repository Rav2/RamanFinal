from pylab import *
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import fit_lorentz
import load_files
import fit_old

#########################################################################
############################# INITIAL VALUES FROM FILE "fit_double_lorentz.py" ##############################
omega0 = 521.604599812  # [cm^(-1)]
# gamma0=3.5 #[cm^(-1)]
gamma0 = 4.42010161931  # [cm^(-1)]
# a0-lattice const of bulk Si
a0 = 0.543102  # [nm]
q0 = 2 * math.pi / a0
offset_from_fit = 0  # z pliku tekstowego
height_from_fit = 71535.2586963  # z pliku tekstowego


#########################################################################
###########################DEFINING FUNCTIONS##########################
def alfa(d):
    return 7.87 * d


def w(q, omega_0):
    return omega_0 - 120 * (q / q0) ** 2


def c(q, d):
    return math.exp(((-1) * q ** 2 * d ** 2) / (2 * alfa(d)))


def intensity(omega, q, d, omega_0, gamma_0):
    return ((4 * pi ** 2 * q ** 2 * c(q, d)) / ((omega - w(q, omega_0)) ** 2 + (0.5 * gamma_0) ** 2))


# total intensity should be independent of q
def total_intensity(omega, d, N, omega_0, gamma_0):
    # calkujemy po dq,N-wybieramy poczatkowaintensity liczbe krokow
    # metoda trapezow
    a = -pi / a0  # 0.0
    b = pi / a0
    h = (b - a) / N
    inten = 0.5 * intensity(omega, a, d, omega_0, gamma_0) + 0.5 * intensity(omega, b, d, omega_0, gamma_0)
    for k in range(1, N):  # w nawiasach dolna i gorna granica
        inten += intensity(omega, a + k * h, d, omega_0, gamma_0)
    return h * inten


def find_max_intensity(d, omega_0, gamma_0):
    omega1 = arange(500, 530, 0.01)
    return max(total_intensity(omega1, d, 100, omega_0, gamma_0))


def normalized_intensity(omega, d, N, omega_0, gamma_0):
    return total_intensity(omega, d, N, omega_0, gamma_0) / find_max_intensity(d, omega_0, gamma_0)


def find_grain_diameter(X, Y, min_size, max_size, size_step, omega_0, gamma_0, peak_height, return_figure, sample_type):
    X_ex = X
    Y_ex = fit_lorentz.perform_fitting(X_ex, Y, omega_0, gamma_0, peak_height, sample_type)

    fitting_result = fit_old.perform_fitting(X_ex, Y_ex, 500, 580)
    single_lorentz_parameters = fitting_result[0]
    #print single_lorentz_parameters
    x_constraint = (X_ex < 535) * (X_ex > 510)
    X_fit = X_ex[x_constraint]
    Y_fit = Y_ex[x_constraint]

    list = arange(min_size, max_size, size_step)  # lista grubosci nanoziaren
    suma = []
    omega_exp = arange(450, 600, 1)
    for diameter in list:
        delta = 0
        delta = np.sum(np.abs(Y_fit - (single_lorentz_parameters[2] * normalized_intensity(X_fit, diameter, 100, omega_0, gamma_0) + single_lorentz_parameters[3])))
        suma.append(delta)

    print min(suma), suma.index(min(suma)) * size_step + min_size
    nanograin_size = suma.index(min(suma)) * size_step + min_size

    if return_figure:
        figure = plt.Figure()
        my_plot = figure.add_subplot(111)
        my_plot.plot(omega_exp,
                     [single_lorentz_parameters[2] * normalized_intensity(omega, nanograin_size, 100, omega_0, gamma_0) + single_lorentz_parameters[3]
                      for omega in omega_exp], label="fitted for %snm" % nanograin_size)
        my_plot.plot(X_fit, Y_fit, ".", label="exp data", color="c")
        figure.suptitle("d:{0:.1f}[nm]   omega:{1:.4f}[cm^-1]   hwhm:{2:.4f}[cm^-1]   inten:{3:.4}[arb]".format(nanograin_size, single_lorentz_parameters[0], single_lorentz_parameters[1], single_lorentz_parameters[2]))
        my_plot.legend(loc=1)
        return nanograin_size, single_lorentz_parameters[0], abs(single_lorentz_parameters[1]), single_lorentz_parameters[2], figure
    else:
        return nanograin_size, single_lorentz_parameters[0], abs(single_lorentz_parameters[1]), single_lorentz_parameters[2]

#########################################################################    
def main():
    #########################################################################
    ############################# LOADING DATA ##############################   
    result = load_files.load_mapping_file("../content/SemiAmorphS1/25.5.txt")
    X_ex = result[2]
    Y_ex = fit_lorentz.perform_fitting(X_ex, result[3], 522, 4, 0, 3)
    d = find_grain_diameter(X_ex, Y_ex, 1., 20., 0.1, 522, 4, 2.13334158e+04, 0, 3)
    print d[0]

    x_constraint = (X_ex < 600) * (X_ex > 450)
    X_fit = X_ex[x_constraint]
    Y_fit = Y_ex[x_constraint]


    # print "Y_fit",Y_fit
    # print "X_fit",X_fit
    # dla malych grubosci ziaren(<10nm) widac znaczace roznice, dla innych- widma sa takie same
    list = arange(1, 10, 0.1)  # lista grubosci nanoziaren
    suma = []
    omega_exp = arange(450, 600, 2)

    for diameter in list:
        delta = 0
        delta = np.sum(np.abs(Y_ex - (
            height_from_fit * total_intensity(X_ex, diameter, 100, omega0, gamma0) / find_max_intensity(diameter,
                                                                                                        omega0,
                                                                                                        gamma0) + offset_from_fit)))
        suma.append(delta)

    print min(suma), suma.index(min(suma)) * 0.1 + 1
    nanograin = suma.index(min(suma)) * 0.1 + 1
    plt.plot(omega_exp, [
        height_from_fit * total_intensity(omega, nanograin, 100, omega0, gamma0) / find_max_intensity(nanograin, omega0,
                                                                                                      gamma0) + offset_from_fit
        for omega in omega_exp], label="fit for %snm" % (nanograin))
    plt.plot(X_fit, Y_fit, "o", label="data", color="c")
    plt.legend(loc=1)
    plt.show()


#########################################################################
if __name__ == "__main__":
    main()
