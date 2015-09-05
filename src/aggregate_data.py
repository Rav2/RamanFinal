import numpy as np
import nanograin_size
# def createMap(name, minimum, maximum):
def agg_data(data, min_size, max_size, size_step, omega_0, gamma_0, amplitude, offset):
    X = data[0]
    Y = data[1]
    raman_shift = data[2]
    intensity = data[3]

    current_Y = Y[0]
    current_X = X[0]
    result = []
    column = []


    for ii in range(0, len(X)):        
        if (Y[ii] != current_Y):

            choose_condition = (Y == current_Y)*(X == current_X)
            func_result = nanograin_size.find_grain_diameter(raman_shift[choose_condition], intensity[choose_condition],
                                                             min_size, max_size, size_step, omega_0, gamma_0, amplitude,
                                                             offset)
            column.append(func_result[0])  # function returns (diameter, figure with plot)
            current_Y = Y[ii]
        if (X[ii] != current_X):
            print("y")
            print(current_Y)
            print("x")
            print(current_X)
            result.append(column)
            column = []
            current_X = X[ii]
        if (ii == (len(X) - 1)):
            choose_condition = (Y == current_Y)*(X == current_X)
            func_result = nanograin_size.find_grain_diameter(raman_shift[choose_condition], intensity[choose_condition],
                                                             min_size, max_size, size_step, omega_0, gamma_0, amplitude,
                                                             offset)
            column.append(func_result[0])  # function returns (diameter, figure with plot)
            result.append(column)
            print("y")
            print(current_Y)
            print("x")
            print(current_X)

    return result



