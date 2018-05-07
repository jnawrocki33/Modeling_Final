import math
import numpy as np
from numpy import linalg as ln
from scipy import optimize as opt
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random


def print_matrix(m):
    dim = len(m)
    for i in range(dim):
        for j in range(dim):
            print(m[i, j], end='  ')
        print("")

#numbers in years. note, timestep = 3months
egg_stage = 3 #time egg to smolt
ocean1 = 2.5 #time smolt to mature in ocean. should end in .5
ocean2 = 0.25 #time mature smolt to get up river
spawn = 1 #irrelevant, spawning fish just die. for dimension
stock = 1 #same as above

#overall survival rates
egg_to_ocean1 = .011
ocean1_to_ocean2 = .102
ocean2_to_spawn = .09 #can also be determined by linear fitting function

num_fecundity = 8000 #num eggs produced

hatchery_disadvantage = 0.33
num_stocked_eggs = 475000 * hatchery_disadvantage
num_stocked_fry = 1000000 * hatchery_disadvantage
num_stocked_parr = 260000 * hatchery_disadvantage


dim = int(egg_stage * 4 + ocean1 * 4 + ocean2 * 4 + spawn + stock)
matrix = np.zeros((dim, dim))

##
# Call to reset the values of all the parameters to normal
# Used after iterating for sensetivity analysis
# MUST UPDATE VALUES HERE WHENEVER WE CHANGE THEM GLOBALLY
def set_defaults():
    global vector_i, egg_stage, ocean1, ocean2, spawn, stock, egg_to_ocean1, ocean1_to_ocean2, ocean2_to_spawn, num_fecundity, num_stock

    #numbers in years, multiply by 4 since timestep = 3 months
    egg_stage = 3 #time for egg to smolt
    ocean1 = 2.5 #time for smolt to mature in ocean
    ocean2 = 0.25 #time for mature smolt to get up river
    spawn = 1 #irrelevant, spawning fish just die. for dimension
    stock = 1 #same as above

    #overall survival rates
    egg_to_ocean1 = .011
    ocean1_to_ocean2 = .102
    ocean2_to_spawn = 0.09

    num_fecundity = 8000 #num eggs
    vector_i = np.array([1000,0,0,0,1000,0,0,0,1000,0,0,0,   1000,0,0,0,1000,0,0,0,1000,0,    0, 0, 0])


##
# Populate the matrix with the global parameter values
def populate_matrix(matrix):
    #xth root of egg_to_ocean1, to account for being applied multiple times
    egg_per_stage_survival = round(egg_to_ocean1 ** (1/float(egg_stage * 4)), 2)
    curr_index = 0
    for i in range(egg_stage * 4):
        matrix[curr_index+1, curr_index] = egg_per_stage_survival
        curr_index += 1

    #xth root of ocean1_to_ocean2, to account for being applied multiple times
    ocean1_per_stage_survival = round(ocean1_to_ocean2 ** (1/float(ocean1 * 4)), 2)

    for i in range(int(ocean1 * 4)):
        matrix[curr_index + 1, curr_index] = ocean1_per_stage_survival
        curr_index += 1


    matrix[curr_index + 1, curr_index] = ocean2_to_spawn

    matrix[0,egg_stage*4 +int(ocean1*4) + 1] = num_fecundity

    #stock eggs
    matrix[0, egg_stage*4 + int(ocean1*4) + 1 + 1] = num_stocked_eggs
    #stock fry and parr together to keep 12 month cycles intact. just adjust for survival rates
    matrix[4, egg_stage*4 + int(ocean1*4) + 1 + 1] = num_stocked_fry *egg_per_stage_survival + num_stocked_parr


def iterate(matrix, vector_i, should_print=False):
    for i in range(1,248):
        if should_print == True:
            print("t = " + str(i))
        vector_i = matrix.dot(vector_i)
        #turn on stocking every 4th iteration (annually)
        if (i+1) % 4 == 0:
            vector_i[len(vector_i) - 1] = 1
        vector_i = vector_i.astype(int)
        skip_eggs(vector_i)
        if should_print == True:
            print(vector_i)
            print("\n--------------------------------------------\n")
    return vector_i

# get survival rate going up river as a function of number of dams
# linear fit at this point
# could try quadratic fit.. but simple is better
def get_spawn_survival(num_dams):
    min_dams = 0
    max_dams = 5
    min_survival = 0.06
    max_survival = 0.085
    slope = float(max_survival - min_survival)/(min_dams-max_dams)
    return (slope * (num_dams - max_dams) + min_survival)

##
# Jump some eggs ahead various amounts of time, to disperse populations
# TODO can add some randomness to how far we jump
def skip_eggs(vector_i):
    egg_per_stage_survival = round(egg_to_ocean1 ** (1/float(egg_stage * 4)), 2)
    pct_to_skip_12_months = 0.25
    num_to_skip = int(vector_i[1] * pct_to_skip_12_months)
    vector_i[1] = vector_i[1] - num_to_skip
    vector_i[5] = vector_i[5] + num_to_skip * (egg_per_stage_survival ** 4)


populate_matrix(matrix)
print_matrix(matrix)
vector_i = np.array([1000,0,0,0,1000,0,0,0,1000,0,0,0,   1000,0,0,0,1000,0,0,0,1000,0,    0, 0, 0])
#iterate(matrix, vector_i)


## this loop is how we'll do one at a time sensetivity analysis
variable_range = [get_spawn_survival(0),get_spawn_survival(1),get_spawn_survival(2),get_spawn_survival(3),get_spawn_survival(4),get_spawn_survival(5)]
for i in variable_range:
    ocean2_to_spawn = i
    populate_matrix(matrix)
    final_vector_i = iterate(matrix, vector_i)
    print(final_vector_i)
    print("------------------------------------------")
    set_defaults()
