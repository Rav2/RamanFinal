from pylab import *
import numpy as np
from math import pi,sqrt
import matplotlib.pyplot as plt
from lmfit import Model, CompositeModel
from lmfit.parameter import Parameters
from lmfit.models import LorentzianModel, LinearModel, ConstantModel
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


def w(q, center):
    return center - 120 * (q / q0) ** 2


def c(q, d):
    return math.exp(((-1) * q ** 2 * d ** 2) / (2 * alfa(d)))


def intensity(omega, q, d, center, sigma):
    return ((4 * pi ** 2 * q ** 2 * c(q, d)) / ((omega - w(q, center)) ** 2 + (sigma) ** 2))

def normalized_intensity(x, d, N, center, sigma):
    return total_intensity(x, d, N, center, sigma) / find_max_intensity(d, center, sigma)

def total_intensity(omega, d, N, center, sigma):
    # calkujemy po dq,N-wybieramy poczatkowaintensity liczbe krokow
    # metoda trapezow
    a = -pi / a0  # 0.0
    b = pi / a0
    h = (b - a) / N
    inten = 0.5 * intensity(omega, a, d, center, sigma) + 0.5 * intensity(omega, b, d, center, sigma)
    for k in range(1, N):  # w nawiasach dolna i gorna granica
        inten += intensity(omega, a + k * h, d, center, sigma)
    return h * inten


def total_intensity_fit_func(x, d, N, center, sigma, amplitude):
    return amplitude*normalized_intensity(x, d, N, center, sigma)


def find_max_intensity(d, center, sigma):
    omega1 = arange(500, 530, 0.01)
    return max(total_intensity(omega1, d, 100, center, sigma))

################################################################################################

def fano_function(x, q, center, sigma, amplitude ):
    epsilon = 2 * (x - center)/sigma
    return amplitude * ((q + epsilon)**2)/(1 + epsilon**2)

def fano_int_func(x, k, diam,  q, center, sigma):
    epsilon = (center - w(x, center))/ sigma
    return ((epsilon + q)**2 )/(1+epsilon**2)* c(k,diam)

def fano(x, diam, No, q, center, sigma):
    a = -pi / a0  # 0.0
    b = pi / a0
    h = (b - a) / No
    inten =  0.5 * fano_int_func(x, a, diam, q, center, sigma) + 0.5 * fano_int_func(x, b, diam, q, center, sigma)
    for k in range(1, No):  # w nawiasach dolna i gorna granica
        inten += fano_int_func(x, a + k * h, diam, q, center, sigma)
    return h*inten

def fano_total(x, diam, No, q, center, sigma, amplitude):
    return amplitude * fano(x, diam, No, q, center, sigma) / fano_find_max_intensity(diam,  q, No, center, sigma)

def fano_find_max_intensity(diam, q, No, center, sigma):
    omega1 = arange(500, 530, 0.01)
    return max(fano(omega1, diam, No, q, center, sigma))



def find_grain_diameter(X, Y, min_size, max_size, size_step, return_figure,  peaks_params, is_bg_active):
    """
    New parameters are stored in peaks_params. The old ones are still present to provide compatibility with old code.
    :param X:
    :param Y:
    :param min_size:
    :param max_size:
    :param size_step:
    :param center:
    :param gamma_0:
    :param peak_height:
    :param return_figure:
    :param sample_type:
    :param peaks_params: [[center, gamma_0, amplitude, from, to, is_visible, is_lorentz, is fano], ...]
    :return:
    """

    best_d = -1
    best_sigma = -1
    best_center = -1
    best_amplitude = -1
    best_q = 5

    X_ex = X
    max_y = max(Y)
    x_bg_range = (X>570)*(X<600)

    lin_mod = LinearModel(prefix="background_", min=0, max=0.5*max_y)
    pars = lin_mod.guess(Y, x=X)
    final_mod = lin_mod

    for ii in range(0, len(peaks_params)):
        if peaks_params[ii][6]==True:
            mod = LorentzianModel(prefix="l" + str(ii) + '_')
            pars.update(mod.make_params())
            pars['l'+str(ii)+'_center'].set(peaks_params[ii][0], min=peaks_params[ii][3], max=peaks_params[ii][4])
            pars['l'+str(ii)+'_sigma'].set(peaks_params[ii][1])
            pars['l'+str(ii)+'_amplitude'].set(peaks_params[ii][2], min=0)#, max=max_y)
        elif peaks_params[ii][7]==True:
            mod = Model(fano_function)
            pars.update(mod.make_params())
            pars['q'].set(-6, min=-6, max=0) #q may be negative as well as positive, the most common values for q: [-1,5]
            pars['center'].set(peaks_params[ii][0], vary=True)
            pars['sigma'].set(peaks_params[ii][1], min=0.1, vary=True)
            pars['amplitude'].set(peaks_params[ii][2], min=0)
        else:
            mod = Model(total_intensity_fit_func)#total_intensity_fit_func(omega, d, N, center, gamma_0, amplitude)
            pars.update(mod.make_params())
            pars['d'].set((max_size-min_size), min=min_size, max=max_size)
            pars['N'].set(100, vary=False)
            pars['center'].set(peaks_params[ii][0], vary=False)
            pars['sigma'].set(peaks_params[ii][1], vary=False)
            pars['amplitude'].set(peaks_params[ii][2], min=0, max=max_y)

        final_mod = final_mod + mod

###########PLOT & RETURN######################################
    my_figure = plt.Figure()
    x_range = (X>400)*(X<600)
    my_plot = my_figure.add_subplot(111)
    legend_labels = []
    for ax in my_figure.axes:
        ax.set_xlabel(u'\u03C9'+' [cm^-1]')
        ax.set_ylabel('Intensity [arb.]')

    my_plot.plot(X[x_range], Y[x_range], "ro")
    legend_labels.append('Data')

    out = final_mod.fit(Y, pars, x=X)
    if not out.success:
        print("!!!!!!!!!!!!!FIT HAVE FAILED!!!!!!!!!!!")
        legend_labels = ['FIT FAILED!']
        my_plot.legend(legend_labels)
        if return_figure:
            return best_d, best_sigma, best_center, best_amplitude, my_figure
        else:
            return best_d, best_sigma, best_center, best_amplitude
    elif not return_figure:
        if(peaks_params[0][6]==True):
            best_sigma = out.best_values['l'+str(0)+'_sigma']
            best_center = out.best_values['l'+str(0)+'_center']
            best_amplitude = out.best_values['l'+str(0)+'_amplitude']
        elif(peaks_params[0][7]==True):
            best_sigma = out.best_values['sigma']
            best_center = out.best_values['center']
            best_amplitude = out.best_values['amplitude']
            best_q = out.best_values['q']
            print(best_q)
        else:
            best_sigma = out.best_values['sigma']
            best_center = out.best_values['center']
            best_amplitude = out.best_values['amplitude']
            best_d = out.best_values['d']
            print (best_d)
        return best_d, best_sigma, best_center, best_amplitude
    else:
        print(out.fit_report(min_correl=0.5))
        legend_labels.append('Best fit')

        my_plot.plot(X[x_range], out.best_fit[x_range], "b-")
        my_plot.plot(X[x_range], out.best_values['background_slope']*X[x_range]+out.best_values['background_intercept'], '--')
        legend_labels.append('background')
    
        for jj in range(0, len(peaks_params)):
            if(peaks_params[jj][6]==True):
                    sigma = out.best_values['l'+str(jj)+'_sigma']
                    center = out.best_values['l'+str(jj)+'_center']
                    amplitude = out.best_values['l'+str(jj)+'_amplitude']
                    legend_labels.append('l'+str(jj))
                    my_plot.plot(X[x_range], lorentzian(X[x_range], center=center, sigma=sigma, amplitude=amplitude), '--')
            elif (peaks_params[jj][7]==True):
                    sigma = out.best_values['sigma']
                    center = out.best_values['center']
                    amplitude = out.best_values['amplitude']
                    q = out.best_values['q']
                    legend_labels.append('fano; q='+str(q))
                    my_plot.plot(X[x_range], fano_function(X[x_range], q,  center, sigma, amplitude), '--')
            else:
                sigma = out.best_values['sigma']
                center = out.best_values['center']
                amplitude = out.best_values['amplitude']
                d = out.best_values['d']
                N = out.best_values['N']
                legend_labels.append('fc; d='+str(d))
                my_plot.plot(X[x_range], total_intensity_fit_func(X[x_range], d, N, center, sigma, amplitude), '--')

        my_plot.legend(legend_labels)
        if(peaks_params[jj][7]==True):
            my_figure.suptitle("q:{0:.1f}[nm]   omega:{1:.4f}[cm^-1]   sigma:{2:.4f}[cm^-1]    inten:{3:.4}[arb]".format(float(best_q), float(best_center), float(best_sigma), float(best_amplitude)))
        else:
            my_figure.suptitle("d:{0:.1f}[nm]   omega:{1:.4f}[cm^-1]   sigma:{2:.4f}[cm^-1]    inten:{3:.4}[arb]".format(float(best_d), float(best_center), float(best_sigma), float(best_amplitude)))
        my_plot.legend(loc=1)
        return best_d, best_sigma, best_center, best_amplitude, my_figure
    
         
#########################################################################    
def main():
  print('Hello!')

#########################################################################
if __name__ == "__main__":
    main()
