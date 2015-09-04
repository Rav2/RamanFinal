__author__ = 'rafal'

import numpy as np

def load_two_column_file(filepath):
    content = np.loadtxt(filepath)
    x = content[:,0]
    y = content[:,1]
    return [[0], [0], x, y]

def load_mapping_file(filepath):
    content = np.loadtxt(filepath)
    x_pos = content[:, 0]
    y_pos = content[:, 1]
    x = content[:, 2]
    y = content[:, 3]
    return [x_pos, y_pos, x, y]