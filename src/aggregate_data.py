import numpy as np
import nanograin_size
from numpy import round
# def createMap(name, minimum, maximum):


def agg_data(data, min_size, max_size, size_step, params, sample_type):
    X = data[0]
    Y = data[1]
    raman_shift = data[2]
    intensity = data[3]

    current_Y = Y[0]
    current_X = X[0]
    result_diameter = []
    column_diameter = []
    result_sigma = []
    column_sigma = []
    result_intensity = []
    column_intensity = []
    result_omega_0 = []
    column_omega_0 = []

    for ii in range(0, len(X)):

        if (Y[ii] != current_Y):

            choose_condition = (Y == current_Y)*(X == current_X)
            func_result = nanograin_size.find_grain_diameter(raman_shift[choose_condition], intensity[choose_condition],
                                                             min_size, max_size, size_step, False, params, True)
            column_diameter.append(func_result[0])  # function returns (diameter, sigma, omega_0, peak_high)
            column_sigma.append(func_result[1])
            column_omega_0.append(func_result[2])
            column_intensity.append(func_result[3])
            current_Y = Y[ii]
        if (X[ii] != current_X):
            print("y")
            print(current_Y)
            print("x")
            print(current_X)
            progress = round(ii/len(X) * 100,1)
            msg = """#####################################################################################"""+\
                    "   PROGRESS:  " + str(progress) + \
                    """    ###################################################################################"""
            print(msg)

            result_diameter.append(column_diameter)
            result_sigma.append(column_sigma)
            result_omega_0.append(column_omega_0)
            result_intensity.append(column_intensity)
            column_diameter = []
            column_sigma = []
            column_omega_0 = []
            column_intensity = []
            current_X = X[ii]
        if (ii == (len(X) - 1)):
            choose_condition = (Y == current_Y)*(X == current_X)
            func_result = nanograin_size.find_grain_diameter(raman_shift[choose_condition], intensity[choose_condition],
                                                             min_size, max_size, size_step, False, params,
                                                             True)
            column_diameter.append(func_result[0])  # function returns (diameter, sigma, omega_0, peak_high)
            column_sigma.append(func_result[1])
            column_omega_0.append(func_result[2])
            column_intensity.append(func_result[3])
            result_diameter.append(column_diameter)
            result_sigma.append(column_sigma)
            result_omega_0.append(column_omega_0)
            result_intensity.append(column_intensity)
            print("y")
            print(current_Y)
            print("x")
            print(current_X)


    #return result_diameter, result_sigma, result_omega_0, result_intensity
    return  result_diameter, result_sigma, result_omega_0, result_intensity


