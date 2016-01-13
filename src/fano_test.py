__author__ = 'rafal'


import numpy as np
import matplotlib.pyplot as plt

def fano_function(x, q, center, sigma, amplitude ):
    epsilon = 2 * (x - center)/sigma
    return amplitude * ((q + epsilon)**2)/(1 + epsilon**2)

q_val = np.linspace(-5,0,5)
center = 521
sigma = 5
amplitude = 90000

X = np.linspace(500,550,50)

legend_vals = []

for q in q_val:
    plt.plot(X, fano_function(X, q, center, sigma, amplitude))
    legend_vals.append(str(q))

plt.legend(legend_vals)
plt.show()
