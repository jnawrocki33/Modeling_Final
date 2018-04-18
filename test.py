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
egg_stage = 1 #time egg to smolt
ocean1 = 0.5 #time smolt to mature in ocean
ocean2 = 0.25 #time mature smolt to get up river
spawn = 1 #irrelevant, spawning fish just die. for dimension
stock = 1 #same as above

#overall survival rates
egg_to_ocean1 = .01
ocean1_to_ocean2 = 0.1
ocean2_to_spawn = .1

dim = int(egg_stage * 4 + ocean1 * 4 + ocean2 * 4 + spawn + stock)
matrix = np.zeros((dim, dim))

#xth root of egg_to_ocean1, to account for being applied multiple times
egg_per_stage_survival = round(egg_to_ocean1 ** (1/float(4)), 2)
curr_index = 0
for i in range(egg_stage * 4):
    matrix[curr_index+1, curr_index] = egg_per_stage_survival
    curr_index += 1

#xth root of ocean1_to_ocean2, to account for being applied multiple times
ocean1_per_stage_survival = round(ocean1_to_ocean2 ** (1/float(2)), 2)

for i in range(int(ocean1 * 4)):
    matrix[curr_index + 1, curr_index] = ocean1_per_stage_survival
    curr_index += 1


matrix[curr_index + 1, curr_index] = ocean2_to_spawn

matrix[0,egg_stage*4 +int(ocean1*4) +1 ] = 8000 #fecundity

matrix[0, egg_stage*4 + int(ocean1*4) + 1 + 1] = 5000 #"stock"


print_matrix(matrix)

vector_i = np.array([1000,0,0,0,   1000,0,    0, 0, 0])

print("t = 0")
print(vector_i)

print("--------------------------------------------")

for i in range(1,300):
    print("t = " + str(i))
    vector_i = matrix.dot(vector_i)
    if (i+1) % 4 == 0:
        vector_i[len(vector_i) - 1] = 1
    print(vector_i)
    print("\n--------------------------------------------\n")
