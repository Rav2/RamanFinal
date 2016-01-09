from pylab import *
import numpy as np
from math import pi,sqrt
import matplotlib.pyplot as plt
import fit_lorentz
import load_files
import fit_old
from lmfit import Model, CompositeModel
from lmfit.models import LorentzianModel, LinearModel
from lmfit.lineshapes import lorentzian
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
A=171400 # [cm^(-2)]
B=100000 # [cm^(-2)]


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

def Fano_omega(k):
    return sqrt(A+B*cos(pi*k/2))

def epsilon(omega,k,gamma_0):
    return (omega-Fano_omega(k))/gamma_0

def Fano_epsilon(omega,omega_0,gamma_0):
    return (omega-omega_0)/(gamma_0/2)

def intensity_Fano(x,k,q,gamma_0):
    return ((q+epsilon(x,k,gamma_0) **2)/(1+epsilon(x,k,gamma_0)) **2 )

# |1/q| - Fano factor
def Fano_intensity(x,q,omega_0,gamma_0,amplitude):
    return amplitude*(((q+Fano_epsilon(x,omega_0,gamma_0)) **2)/(q+(Fano_epsilon(x,omega_0,gamma_0) )**2 ))

def total_intensity_fit_Fano(x, q, N,  gamma_0, amplitude):
    return amplitude*normalized_Fano_intensity(x, q, N,  gamma_0)

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

def total_Fano_intensity(omega, q, N, gamma_0):
    # calkujemy po dq,N-wybieramy poczatkowaintensity liczbe krokow
    # metoda trapezow
    a = -pi / a0  # 0.0
    b = pi / a0
    h = (b - a) / N
    inten = 0.5 * intensity_Fano(omega, a,  q, gamma_0) + 0.5 * intensity_Fano(omega, b, q,  gamma_0)
    for k in range(1, N):  # w nawiasach dolna i gorna granica
        inten += intensity_Fano(omega, a + k * h, q, gamma_0)
    return h * inten


def total_intensity_fit_func(x, d, N, omega_0, gamma_0, amplitude):
    return amplitude*normalized_intensity(x, d, N, omega_0, gamma_0)


def find_max_intensity(d, omega_0, gamma_0):
    omega1 = arange(500, 530, 0.01)
    return max(total_intensity(omega1, d, 100, omega_0, gamma_0))

def find_max_Fano_intensity(q,  gamma_0):
    omega1 = arange(500, 530, 0.01)
    return max(total_Fano_intensity(omega1, q, 100,  gamma_0))


def normalized_intensity(x, d, N, omega_0, gamma_0):
    return total_intensity(x, d, N, omega_0, gamma_0) / find_max_intensity(d, omega_0, gamma_0)

def normalized_Fano_intensity(x, q, N,  gamma_0):
    return total_Fano_intensity(x, q, N,  gamma_0) / find_max_Fano_intensity(q,  gamma_0)

def find_grain_diameter(X, Y, min_size, max_size, size_step, omega_0, gamma_0, peak_height, return_figure, sample_type,  peaks_params):
    """
    New parameters are stored in peaks_params. The old ones are still present to provide compatibility with old code.
    :param X:
    :param Y:
    :param min_size:
    :param max_size:
    :param size_step:
    :param omega_0:
    :param gamma_0:
    :param peak_height:
    :param return_figure:
    :param sample_type:
    :param peaks_params: [[omega_0, gamma_0, amplitude, from, to, is_visible, is_lorentz], ...]
    :return:
    """
    all_lorentz = True
    #for i in range(0, len(peaks_params)):

        #all_lorentz *= peaks_params[i][6]

    X_ex = X #dlaczego
    Y_ex = fit_lorentz.perform_fitting(X_ex, Y, peaks_params)

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

    print (min(suma), suma.index(min(suma)) * size_step + min_size)
    nanograin_size = suma.index(min(suma)) * size_step + min_size


        ###mozna przerobic program tak, zeby sprawdzal czy dany pik nie jest lorentzowski i wtedy dodawal instancje Model
        ### utworzona w oparciu o nasza funkcje do szukania rozmiaru ziaren
    max_y = max(Y)
    lin_mod = LinearModel(prefix="background_", min=0, max=0.5*max_y)
    pars = lin_mod.guess(Y, x=X)
    final_mod = lin_mod
    print "to peaks_params",peaks_params
    for ii in range(0,len(peaks_params)):
        if peaks_params[ii][6]==True:
            mod = LorentzianModel(prefix="l" + str(ii) + '_')
            pars.update(mod.make_params())
            pars['l'+str(ii)+'_center'].set(peaks_params[ii][0], min=peaks_params[ii][3], max=peaks_params[ii][4])
            ###TODO: sigma ponizej jest tozsama z gamma podanym przez uzytkownika, a powinno byc inaczej (patrz dokumentacja lmfit)!!! Uzytkownik powinien podawac par. sigma. (Dla poprawnosci dopasowania to wszystko jedno ale dla zrozumienia wynikow nie.)
            pars['l'+str(ii)+'_sigma'].set(peaks_params[ii][1])
            pars['l'+str(ii)+'_amplitude'].set(peaks_params[ii][2], min=0)#, max=max_y)
        elif peaks_params[ii][7]==True:
            mod = Model(Fano_intensity)
            print "fano"
            pars.update(mod.make_params())
            pars['q'].set(-5.25, min=-5.5,max=-5) #q may be negative as well as positive, the most common values for q: [-1,5]
            #pars['omega_0'].set(peaks_params[ii][0], vary=False)
            #pars['gamma_0'].set(peaks_params[ii][1], vary=False)
            #pars['amplitude'].set(peaks_params[ii][2], min=0)
            pars['omega_0'].set(515.6,min=510,max=525)
            pars['gamma_0'].set(2.3,min=0)
            #pars['amplitude'].set(max_y, min=0, max=max_y)
            pars['amplitude'].set(max_y, min=0,max=max_y)
        else:
            mod = Model(total_intensity_fit_func)#total_intensity_fit_func(omega, d, N, omega_0, gamma_0, amplitude)
            pars.update(mod.make_params())
            pars['d'].set(10, min=0, max=30)
            pars['N'].set(100, vary=False)
            pars['omega_0'].set(peaks_params[ii][0], vary=False)
            pars['gamma_0'].set(peaks_params[ii][1], vary=False)
            pars['amplitude'].set(0.8*max_y, min=0, max=max_y)
        final_mod = final_mod + mod
                # print("adding a peak!")
            ####PLOT###
    my_figure = plt.Figure()
    x_range = (X>400)*(X<600)
    my_plot = my_figure.add_subplot(111)
    legend_labels = []
    for ax in my_figure.axes:
        ax.set_xlabel(u'\u03C9'+' [cm^-1]')
        ax.set_ylabel('Intensity [arb.]')

    my_plot.plot(X[x_range], Y[x_range], "ro")
    legend_labels.append('Data')
            # init = final_mod.eval(pars, x=X)
            # my_plot.plot(X[x_range], init[x_range], 'k--')

    out = final_mod.fit(Y, pars, x=X)
    if not out.success:
        print("!!!!!!!!!!!!!FIT HAVE FAILED!!!!!!!!!!!")
        legend_labels = ['FIT FAILED!']
        my_plot.legend(legend_labels)
        return nanograin_size, single_lorentz_parameters[0], abs(single_lorentz_parameters[1]), single_lorentz_parameters[2], my_figure

    else:
        print(out.fit_report(min_correl=0.5))
        legend_labels.append('Best fit')

    my_plot.plot(X[x_range], out.best_fit[x_range], "b-")
        
    my_plot.plot(X[x_range], out.best_values['background_slope']*X[x_range]+out.best_values['background_intercept'], '--')
    legend_labels.append('background')
    for jj in range(0, len(peaks_params)):
        if(peaks_params[jj][6]==True):
            if peaks_params[jj][6]==True:
                sigma = out.best_values['l'+str(jj)+'_sigma']
                center = out.best_values['l'+str(jj)+'_center']
                amplitude = out.best_values['l'+str(jj)+'_amplitude']
                legend_labels.append('l'+str(jj))
                my_plot.plot(X[x_range], lorentzian(X[x_range], center=center, sigma=sigma, amplitude=amplitude), '--')
            else:
                sigma = out.best_values['gamma_0']
                center = out.best_values['omega_0']
                amplitude = out.best_values['amplitude']
                #d = out.best_values['d']
                q = out.best_values['q']
                #N = out.best_values['N']
                legend_labels.append('fc'+str(jj))
                #my_plot.plot(X[x_range], total_intensity_fit_func(X[x_range], d, N, center, sigma, amplitude), '--')
                my_plot.plot(X[x_range], Fano_intensity(X[x_range], q,  center, sigma, amplitude), '--')
    my_plot.legend(legend_labels)
    return nanograin_size, single_lorentz_parameters[0], abs(single_lorentz_parameters[1]), single_lorentz_parameters[2], my_figure

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
    result = load_files.load_mapping_file("../content/2015-04-10Si-mapySingleFiles/SemiAmorphS1/25.5.txt")
    X_ex = result[2]
    Y_ex = fit_lorentz.perform_fitting(X_ex, result[3], 522, 4, 0, 3)
    d = find_grain_diameter(X_ex, Y_ex, 1., 20., 0.1, 522, 4, 2.13334158e+04, 0, 3)
    print (d[0])

    x_constraint = (X_ex < 600) * (X_ex > 450)
    X_fit = X_ex[x_constraint]
    Y_fit = Y_ex[x_constraint]


    # print "Y_fit",Y_fit
    # print "X_fit",X_fit
    # dla malych grubosci ziaren(<10nm) widac znaczace roznice, dla innych- widma sa takie same
    list = arange(1, 10, 0.1)  # lista grubosci nanoziaren
    suma = []
    omega_exp = arange(450, 600, 2)

    # for diameter in list:
    #     delta = 0
    #     delta = np.sum(np.abs(Y_ex - (
    #         height_from_fit * total_intensity(X_ex, diameter, 100, omega0, gamma0) / find_max_intensity(diameter,
    #                                                                                                     omega0,
    #                                                                                                     gamma0) + offset_from_fit)))
    #     suma.append(delta)

    # print (min(suma), suma.index(min(suma)) * 0.1 + 1)
    # nanograin = suma.index(min(suma)) * 0.1 + 1
    # plt.plot(omega_exp, [
    #     height_from_fit * total_intensity(omega, nanograin, 100, omega0, gamma0) / find_max_intensity(nanograin, omega0,
    #                                                                                                   gamma0) + offset_from_fit
    #     for omega in omega_exp], label="fit for %snm" % (nanograin))
    # plt.plot(X_fit, Y_fit, "o", label="data", color="c")
    # plt.legend(loc=1)
    plt.show()


#########################################################################
if __name__ == "__main__":
    main()
