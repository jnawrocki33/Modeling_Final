#turn on and off but use numerical sensetivity analysis

import math
import numpy as np
from numpy import linalg as ln
from scipy import optimize as opt
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def print_matrix(m):
    dim = len(m)
    for i in range(dim):
        for j in range(dim):
            print(m[i, j], end='  ')
        print("")

#numbers in years, multiply by 4 since timestep = 3 months
egg_stage = 3 #time egg to smolt
ocean1 = 3 #time smolt to mature in ocean
ocean2 = 0.25 #time mature smolt to get up river
spawn = 0 #irrelevant, spawning fish just die
stock = 1

#overall survival rates
egg_to_ocean1 = .01
ocean1_to_ocean2 = 0.1
ocean2_to_spawn = .1

dim = int(egg_stage * 4 + ocean1 * 4 + ocean2 * 4 + 1 + stock * 4)
#dim = 12 + 12 + 1 + 1 + 4
matrix = np.zeros((dim, dim))
#print(matrix)

egg_per_stage_survival = round(egg_to_ocean1 ** (1/float(12)), 2)
for i in range(12):
    matrix[i+1, i] = egg_per_stage_survival

ocean1_per_stage_survival = round(ocean1_to_ocean2 ** (1/float(12)), 2)
for i in range(12):
    matrix[12 + i + 1, 12 + i] = ocean1_per_stage_survival


matrix[25, 24] = 1
matrix[26,25] = 1
matrix[27,26] = 1
matrix[28,27] = ocean2_to_spawn

matrix[0,12+12+4] = 8000 #fecundity


print_matrix(matrix)
