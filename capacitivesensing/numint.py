#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 14:35:05 2024

@author: vinzenz
"""

import numpy as np
import matplotlib.pyplot as plt

"""global constants"""
#charge and charge density
#no radius because cancels with r' from r'dphi' 
q = 1
qdens1 = q/(2*np.pi)
qdens2 = 0#-q/(2*np.pi)

#constants
epsilon = 8.854e-12 #farads per meter
constant1 = qdens1/(4*np.pi*epsilon)
constant2 = qdens2/(4*np.pi*epsilon)





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
def integrate(f, R, r, z):
    # Define integration limits
    a = 0
    b = 2*np.pi
    #number of subintervals
    N = 1000
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


def Efield(R,r,z, d = 0):
    Er = integrate(integrandEr, R, r, z)
    Ez = integrate(integrandEz, R, r, z)
    if d!=0:
        Er -= integrate(integrandEr, R, r, z-2*d)
        Ez -= integrate(integrandEz, R, r, z-2*d)

    return np.array([Er, Ez])


# Define the vector field
def vectorfield(R1, R2, rr, zz, d=0):
    E1 = constant1* Efield(R1, rr, zz, d)
    E2 = constant2* Efield(R2, rr, zz, d)
    return E1+E2


def plotfield(R1, R2, d):
    #plotting window
    r_min = 0
    r_max = 0.04
    z_min = -d
    z_max = 3*d

    if d==0:
        z_min = -0.04
        z_max = 0.04

    ndots = 5*4#multiples of 4!!!
    h = (z_max-z_min)/ndots

    # Define the grid
    r = np.linspace(r_min, r_max, ndots)
    z = np.array([i * h for i in range(int(-ndots/4), int(3/4*ndots+2))])-h/2
    rr, zz = np.meshgrid(r, z)
    
    # Compute the vector field values
    U, V = vectorfield(R1, R2, rr, zz, d)
    
    # Plot the vector field
    plt.figure(figsize=(6, 6))
    plt.quiver(rr, zz, U, V, color='b', angles='xy', scale_units='xy', scale=1e17)
    if d!= 0:
        plt.hlines(d, 0, 0.04, color='red')
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title('Electric field in the (r-z)-plane')
    plt.grid()
    plt.show()


#compute voltage difference
def capacitance(R1, R2, d=0):
    N = 1000
    x = np.linspace(R1, R2, N+1)
    dx = (R2-R1)/N
    E1 = constant1*Efield(R1,x,0, d)
    E2 = constant2*Efield(R2,x,0, d)
    Er = E1[0] + E2[0]
    voltage = np.sum(Er*dx)
    return q/voltage


#radii of the two rings with R1<R2
R1 = 10e-3#m
s = 10e-3#m
R2 = R1+s#m

#distance of the conducting plane
#choose d=0 for no conductor
d= 13e-3#m

plotfield(R1, R2, d)
print('capacitance: ', capacitance(R1, R2, d),'farads')


d = np.arange(1e-3,40e-3, 1e-3)
c = np.zeros_like(d, dtype=float)
for i, di in enumerate(d):
    c[i] = capacitance(R1, R2, di)

plt.figure()
plt.plot(d*1e3,c*1e12)
plt.title(f'R = {R1 * 1e3} mm, s = {s * 1e3} mm')
plt.xlabel('distance conductor [mm]')
plt.ylabel('Capacitance [pF]')
plt.savefig('fig.pdf')  # Save the figure as a png
plt.show()
plt.close()

