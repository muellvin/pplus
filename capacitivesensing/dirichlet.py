#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 13:01:22 2024

@author: vinzenz
"""

import numpy as np
import math
import matplotlib.pyplot as plt



def set_params(R1value, svalue, dvalue):
    global epsilon
    global R1
    global s
    global d
    global Rout
    global voltage_sensor
    global Nr
    global Nz
    global n
    global m
    global dr
    global dz
    global r
    global z
    global nR1
    global nR1s
    
    epsilon = 8.854e-12 #farads per meter
    voltage_sensor = 1

    R1 = R1value
    s = svalue
    d = dvalue
    Rout = 3*R1 

    #spacing of the finite elements
    Nr = 100
    Nz = 100

    #nxm matrix
    n = Nr+1
    m = Nz+1

    #spacing
    dr = Rout/Nr
    dz = d/Nz

    #vectors containing the coordinates
    r = np.arange(0,n)*dr
    z = np.arange(0,m)*dz

    #special indices:
    nR1 = math.floor(R1/Rout*Nr)
    nR1s = math.ceil((R1+s)/Rout*Nr)
    



#make matrix into vector
def flatten(matrix):
    return matrix.ravel()

def mapindex(i, j):
    """
    Maps the indices (i, j) of a 2D array to the corresponding index in the flattened array.

    Parameters:
    i (int): Row index starting at 0
    j (int): Column index starting at 0
    m (int): Number of columns in the 2D array.

    Returns:
    int: The corresponding index in the flattened array.
    """
    return i * m + j


def voltage_matrix():
    global voltage

    voltage = np.empty((n,m), dtype=float)
    voltage_vec = flatten(voltage)
    b = np.zeros(n*m, dtype=float)
    
    
    #laplacian operator
    L = np.zeros((n*m, n*m), dtype=float)
    
    #fill the operator,
    for i in range(n):
        for j in range(m):
            
            #for boundaries
    
            #preset voltage
            if j == Nz:
                L[mapindex(i, j), mapindex(i, j)] = 1
                b[mapindex(i, j)] = 0
                
                
            elif j == 0:
                
                #preset voltage
                if i <= nR1:
                    L[mapindex(i, j), mapindex(i, j)] = 1
                    b[mapindex(i, j)] = voltage_sensor
                elif i >= nR1s:
                    L[mapindex(i, j), mapindex(i, j)] = 1
                    b[mapindex(i, j)] = 0
                    
                    
                else:
                    L[mapindex(i, j), mapindex(i+1, j)] += 1/(r[i]*dr) + 1/dr**2
                    L[mapindex(i, j), mapindex(i, j)] += -1/(r[i]*dr) -2/dr**2
                    L[mapindex(i, j), mapindex(i-1, j)] += 1/dr**2
                    L[mapindex(i, j), mapindex(i, j+1)] += 1/dz**2
                    L[mapindex(i, j), mapindex(i, j)] += -1/dz**2 #here only -1
                    #L[mapindex(i, j), mapindex(i, j-1)] += 1/dz**2
    
                     
                    
            #to the very left with j inbetween     
            #nothing in radial direction
            elif i == 0:
                L[mapindex(i, j), mapindex(i+1, j)] +=  1/dr**2 #1/(r[i]*dr)
                L[mapindex(i, j), mapindex(i, j)] += -1/dr**2 #-1/(r[i]*dr)
                #L[mapindex(i, j), mapindex(i-1, j)] += 1/dr**2
                L[mapindex(i, j), mapindex(i, j+1)] += 1/dz**2
                L[mapindex(i, j), mapindex(i, j)] += -2/dz**2
                L[mapindex(i, j), mapindex(i, j-1)] += 1/dz**2
            
            #to the very right
            elif i == Nr:
                #L[mapindex(i, j), mapindex(i+1, j)] += 1/(r[i]*dr) + 1/dr**2
                L[mapindex(i, j), mapindex(i, j)] += -1/(r[i]*dr) -1/dr**2
                L[mapindex(i, j), mapindex(i-1, j)] += 1/dr**2
                L[mapindex(i, j), mapindex(i, j+1)] += 1/dz**2
                L[mapindex(i, j), mapindex(i, j)] += -2/dz**2
                L[mapindex(i, j), mapindex(i, j-1)] += 1/dz**2
                
                
            #the usual laplacian
            else:
                L[mapindex(i, j), mapindex(i+1, j)] += 1/(r[i]*dr) + 1/dr**2
                L[mapindex(i, j), mapindex(i, j)] += -1/(r[i]*dr) -2/dr**2
                L[mapindex(i, j), mapindex(i-1, j)] +=  1/dr**2
                L[mapindex(i, j), mapindex(i, j+1)] += 1/dz**2
                L[mapindex(i, j), mapindex(i, j)] += -2/dz**2
                L[mapindex(i, j), mapindex(i, j-1)] += 1/dz**2
    
    
    voltage_vec = np.linalg.solve(L, b)
    
    voltage = np.reshape(voltage_vec, voltage.shape)



def printvoltage():

    # Create the heatmap
    plt.pcolormesh(r, z, np.transpose(voltage), cmap='hot')
    plt.colorbar()  # Add a color bar indicating the scale
    plt.xlabel('r Coordinate')
    plt.ylabel('z Coordinate')
    plt.title('Heatmap of Voltage')
    plt.gca().set_aspect('equal')  # Set aspect ratio to 'auto'
    
    plt.show()


def capacitance():

    """compute the electric field and the energy stored in it"""
    global Er
    global Ez
    global W
    global C
    
    Er = np.zeros((n,m), dtype=float)
    Ez = np.zeros((n,m), dtype=float)
    
    W=0
    
    #computation of the del operator
    for i in range(n-1):
        for j in range(m-1):
            Er[i,j] = - (voltage[i+1,j]-voltage[i,j])/(dr)
            Ez[i,j] = - (voltage[i,j+1]-voltage[i,j])/(dz)
            
            W += 1/2 * epsilon * (Er[i,j]**2 + Ez[i,j]**2) * 2*np.pi * r[i] * dr * dz
    
    """capacitance"""
    #W = 1/2 * C *V**2
    C = 2*W/voltage_sensor**2
    print('capacitance: {:.2f} pF'.format(C*1e12))
    return C

def drawEfield():
    # Create a grid of coordinates for the vectors
    rr, zz = np.meshgrid(r,z)
    
    magnitude = np.sqrt(Er**2+Ez**2)
    Er_normed = Er/magnitude
    Ez_normed = Ez/magnitude
    
    # Plot the vector field
    plt.figure()
    plt.quiver(rr[::10,::10], zz[::10,::10], np.transpose(Er_normed[::10,::10]), np.transpose(Ez_normed[::10,::10]))
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Vector Field')
    plt.show()

def computation(R1value, svalue, dvalue):
    set_params(R1value, svalue, dvalue)
    voltage_matrix()
    #printvoltage()
    c = capacitance()
    #drawEfield()
    return c





dvalues = np.array([0.5,1,2,3,5,10])*1e-3

mil = 0.0254e-3

R1values = np.array([125,225,425,625])*mil
s = 195*mil

A_C = np.empty(len(dvalues), dtype=float)

for i, dvalue in enumerate(dvalues):
    A_C[i] = computation(20e-3, s, dvalue)

print(A_C)



"""
# Create the heatmap
plt.pcolormesh(r, z, np.transpose(Er), cmap='hot')
plt.colorbar()  # Add a color bar indicating the scale
plt.xlabel('r Coordinate')
plt.ylabel('z Coordinate')
plt.title('Heatmap of Er')
plt.gca().set_aspect('equal')  # Set aspect ratio to 'auto'

plt.show()

# Create the heatmap
plt.pcolormesh(r, z, np.transpose(Ez), cmap='hot')
plt.colorbar()  # Add a color bar indicating the scale
plt.xlabel('r Coordinate')
plt.ylabel('z Coordinate')
plt.title('Heatmap of Ez')
plt.gca().set_aspect('equal')  # Set aspect ratio to 'auto'

plt.show()
"""


