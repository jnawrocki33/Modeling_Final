import math
import numpy as np
from numpy import linalg as ln
from scipy import optimize as opt
from scipy.integrate import odeint
import matplotlib.pyplot as plt



J = 128
L = 128
N = (L + 1) * (J + 1)

A = np.zeros((N,N))
b = np.zeros((N))
T_0 = 1

for j in range(J+1):
    for l in range(L+1):
        I = j*(L+1)+l
        if (l == 0) or (j == 0 and l != J) or (j == J and l != L):
            A[I][I] = 1
            b[I] = 0
        elif (l == L):
            #print(j)

            A[I][I] = 1

            if j <= ((J)/2):
                b[I] = (2*T_0) * (j/J)
            if j > ((J)/2):
                b[I] = (2*T_0) * (1 - (j/(J)))

        else:
            A[I][I] = -4.0 
            A[I][I+1] = 1.0
            A[I][I-1] = 1.0
            A[I][I+L+1] = 1.0
            A[I][I-(L+1)] = 1.0
            b[I] = 0

print('A =', A, '\n')
print('b =', b, '\n')
T = ln.solve(A, b)
print('T =', T)



def plotter(T):
    
    length = len(T) 
    temps = []
    x_list = []
    for I in range(length):
        j = I // (L+1)
        l = I % (L+1)
        
        if l == L/2:
            #print('(j,l) =', (j,l))
            #print('I =', I)
            #print('T =', T[I])
            x_list.append(j/ J)
            temps.append(T[I])
                     
    x_data = []
    temps_1 = []
    x = 0
    
    while x < 1:
        n = 1
        summation = 0
        x_data.append(x)
        while n < 50:
            summation += ((math.sin(n*math.pi/2))/n**2)*math.sin(math.pi*n*x)*((math.sinh(math.pi*n/2))/math.sinh(n*math.pi))
            n += 1
        temps_1.append((8*T_0/(math.pi)**2)*summation)
        x += 0.01

    plt.plot(x_data, temps_1, label = 'Analytical')     
    plt.plot(x_list, temps, label ='Numerical')
    plt.xlabel('x / a')
    plt.ylabel('T')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            
    #print(x_list)
    #print(temps)
    
    plt.show()
            


plotter(T)

        
        


    
    
    
    
