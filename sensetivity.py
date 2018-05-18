import math
import numpy as np
from numpy import linalg as ln
from scipy import optimize as opt
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
from scipy.linalg import eig

##
# This is the file we used to get sensetivity/elasticity values for the matrix
# It's essentially the same as main.py, but we just turn off stocking and
# calculate the elasticity matrix using the formula from the book.

np.set_printoptions(suppress=True)
def print_matrix(m):
    dim = len(m)
    for i in range(dim):
        for j in range(dim):
            print(m[i, j], end='  ')
        print("")

def print_wolfram(m):
    dim = len(m)
    print("[", end = '')
    for i in range(dim):
        print("[", end = '')
        for j in range(dim):
            if j != dim - 1:
                print(m[i,j], end = ', ')
            else:
                print(m[i,j], end = '')
        print("]", end = ',')
    print("]", end = '')

#numbers in years. note, timestep = 3months
egg_stage = 3 #time egg to smolt
ocean1 = 2.5 #time smolt to mature in ocean. should end in .5
ocean2 = 0.25 #time mature smolt to get up river
spawn = 1 #irrelevant, spawning fish just die. for dimension
stock = 1 #same as above

#overall survival rates
egg_to_ocean1 = .011
ocean1_to_ocean2 = .079
ocean2_to_spawn = .09 #can also be determined by linear fitting function

num_fecundity = 4000 #num eggs produced

# hatchery_disadvantage = 0.33
# num_stocked_eggs = 475000 * hatchery_disadvantage
# num_stocked_fry = 1000000 * hatchery_disadvantage
# num_stocked_parr = 260000 * hatchery_disadvantage


dim = int(egg_stage * 4 + ocean1 * 4 + ocean2 * 4 + spawn)# + stock)
matrix = np.zeros((dim, dim))

##
# Call to reset the values of all the parameters to normal
# Used after iterating for sensetivity analysis
# MUST UPDATE VALUES HERE WHENEVER WE CHANGE THEM GLOBALLY
# def set_defaults():
#     global vector_i, egg_stage, ocean1, ocean2, spawn, stock, egg_to_ocean1, ocean1_to_ocean2, ocean2_to_spawn, num_fecundity, num_stock
#
#     #numbers in years, multiply by 4 since timestep = 3 months
#     egg_stage = 3 #time for egg to smolt
#     ocean1 = 2.5 #time for smolt to mature in ocean
#     ocean2 = 0.25 #time for mature smolt to get up river
#     spawn = 1 #irrelevant, spawning fish just die. for dimension
#     stock = 1 #same as above
#
#     #overall survival rates
#     egg_to_ocean1 = .011
#     ocean1_to_ocean2 = .102
#     ocean2_to_spawn = 0.09
#
#     num_fecundity = 4000 #num eggs
#     vector_i = np.array([1000,0,0,0,1000,0,0,0,1000,0,0,0,   1000,0,0,0,1000,0,0,0,1000,0,    0, 0])


##
# Populate the matrix with the global parameter values
def populate_matrix(matrix):
    #xth root of egg_to_ocean1, to account for being applied multiple times
    egg_per_stage_survival = round(egg_to_ocean1 ** (1/float(egg_stage * 4)), 4)
    curr_index = 0
    for i in range(egg_stage * 4):
        matrix[curr_index+1, curr_index] = egg_per_stage_survival * 1.1
        curr_index += 1

        #increasing by 1.1 for testing eigenvalue behavior
    # matrix[1,0] = egg_per_stage_survival + .1 * egg_per_stage_survival
    # matrix[2,1] = egg_per_stage_survival + .1 * egg_per_stage_survival
    # matrix[3,2] = egg_per_stage_survival + .1 * egg_per_stage_survival
    # matrix[4,3] = egg_per_stage_survival + .1 * egg_per_stage_survival
    # matrix[5,4] = egg_per_stage_survival + .1 * egg_per_stage_survival
    # matrix[6,5] = egg_per_stage_survival + .1 * egg_per_stage_survival
    # matrix[7,6] = egg_per_stage_survival + .1 * egg_per_stage_survival
    # matrix[8,7] = egg_per_stage_survival + .1 * egg_per_stage_survival

    #original lambda: (0.9527051969088003+0j)

    #changing one by .1: (0.9564961570166772+0j)
    #diff is: .00379, approx .1 * orig_lambda * elasticity (=.00397)

    #changing by .2: (0.9599702027305008+0j)
    #diff is: .00726, approx .1 * orig_lambda * elasticity

    #changing 2 by .1: (0.9603022019362952+0j)
    #diff is: .007597, approx 2 * .1 * orig_lambda * elasticity (=.00794)


    #note, need .05 more in lambda to get to 1. .05/.00397 = 12.5
    #increase all 12 by .1, we get: (0.9992056402156342+0j).. ALMOST 1!
    #increase by .105 => (1.0014739873805125+0j), we did it!



    '''
    elasticity is (delta lambda / lambda) / (delta a_ij / aij)
    increase aij by say .1 (mult by 1.1) -- expect (delta lambda) = (lambda)(elasticity)(.1)

    its multiplicitavely the same, not the same in absolute change
    '''




    #xth root of ocean1_to_ocean2, to account for being applied multiple times
    ocean1_per_stage_survival = round(ocean1_to_ocean2 ** (1/float(ocean1 * 4)), 4)

    for i in range(int(ocean1 * 4)):
        if i == 0:
            matrix[curr_index + 1, curr_index] = ocean1_per_stage_survival
        else:
            matrix[curr_index + 1, curr_index] = ocean1_per_stage_survival
        curr_index += 1


    matrix[curr_index + 1, curr_index] = ocean2_to_spawn

    matrix[0,egg_stage*4 +int(ocean1*4) + 1] = num_fecundity


def iterate(matrix, vector_i, should_print=False):
    for i in range(1,500):
        if should_print == True:
            print("t = " + str(i))
        vector_i = matrix.dot(vector_i)
        vector_i = vector_i.astype(int)
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


populate_matrix(matrix)
print_matrix(matrix)
vector_i = np.array([1000,0,0,0,1000,0,0,0,1000,0,0,0,   1000,0,0,0,1000,0,0,0,1000,0,    0, 0])# 0])

#elasticity stuff
vals, left_eigenvectors, right_eigenvectors = eig(matrix, left=True)
max = 0
max_val = 0
for i in range(len(vals)):
    if vals[i] > max_val:
        max_val = vals[i]
        max = i
print(vals)
print(max_val)
print(max)
left = left_eigenvectors[:,max]
right = right_eigenvectors[:,max]
print(left)
print(right)
E = np.zeros((dim, dim))
dot = np.dot(left, right)
print(dot)
for i in range(dim):
    for j in range(dim):
        E[i,j] = (right[j] * left[i]) / dot * matrix[i,j] / max_val
print_matrix(E)
print('')
print_matrix(matrix)
exit()


# iterate(matrix, vector_i, should_print=True)
# exit()
#iterate(matrix, vector_i)


## this loop is how we'll do one at a time sensetivity analysis
# variable_range = [get_spawn_survival(0),get_spawn_survival(1),get_spawn_survival(2),get_spawn_survival(3),get_spawn_survival(4),get_spawn_survival(5)]
# for i in variable_range:
#     ocean2_to_spawn = i
#     populate_matrix(matrix)
#     final_vector_i = iterate(matrix, vector_i)
#     print(final_vector_i)
#     print("------------------------------------------")
#     set_defaults()
