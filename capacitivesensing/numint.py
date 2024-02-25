#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 14:35:05 2024

@author: vinzenz
"""

import numpy as np
import matplotlib.pyplot as plt






#Gauss_quadrature integrating f(R,r,z,x) with respect to x on [a,b]
def gauss_quadrature(f, R, r, z, a, b):
    # Gauss two-point quadrature rule
    # Nodes and weights
    nodes = np.array([-np.sqrt(1/3), np.sqrt(1/3)])
    weights = np.array([1, 1])
    
    # Map nodes from [-1, 1] to [a, b]
    mapped_nodes = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    
    # Calculate integral approximation
    integral = 0
    for i in range(len(nodes)):
        integral += weights[i] * f(R, r, z, mapped_nodes[i])
        
    integral *= 0.5 * (b - a)
    return integral

#Integration of f(R,r,z,x) with respect to x on [a,b]
#using N subintervals and performing gauss_quadrature on it
def integrate(f, R, r, z, a, b, N):
    integral = 0
    x = np.linspace(a,b, N+1)
    for i in range(N):
        integral += gauss_quadrature(f, R, r, z, x[i], x[i+1])
    return integral

# integrand; radius R, position (r,z), integration variable phi dash = x
def integrandEr(R, r, z, x):
    numerator = r-R*np.cos(x)
    d = np.sqrt(r**2 + R**2 + z**2 - 2*r*R*np.cos(x))
    denominator = d**(3)
    return numerator/denominator

#only for sanity check, should be zero
def integrandEy(R, r, z, x):
    numerator = R*np.sin(x)
    d = np.sqrt(r**2 + R**2 + z**2 - 2*r*R*np.cos(x))
    denominator = d**(3)
    return numerator/denominator

def integrandEz(R, r, z, x):
    numerator = z
    d = np.sqrt(r**2 + R**2 + z**2 - 2*r*R*np.cos(x))
    denominator = d**(3)
    return numerator/denominator



def Efield(R,r,z, a, b, N):
    Er = integrate(integrandEr, R, r, z, a, b, N)
    Ez = integrate(integrandEz, R, r, z, a, b, N)
    return np.array([Er, Ez])



# Define integration limits
a = 0
b = 2*np.pi
N = 1000
R1 = 0.7
R2 = 3.2
constant1 = 1
constant2 = -1

#z = 10
#r = 40
#print(integrate(integrandEy, R, r, z, a, b, N))

# Define the grid
r = np.linspace(0, 4, 20)
z = np.linspace(-4, 4, 20)
rr, zz = np.meshgrid(r, z)

# Define the vector field
def vectorfield(rr,zz, a,b,N):
    E1 = constant1* Efield(R1, rr,zz, a,b,N)
    E2 = constant2* Efield(R2, rr,zz, a,b,N)
    return E1+E2

# Compute the vector field values
U, V = vectorfield(rr, zz, a,b,N)

# Plot the vector field
plt.figure(figsize=(6, 6))
plt.quiver(rr, zz, U, V, color='b', angles='xy', scale_units='xy', scale=10)
plt.xlabel('R')
plt.ylabel('Z')
plt.title('Electric field in the (r-z)-plane')
plt.grid()
plt.show()
