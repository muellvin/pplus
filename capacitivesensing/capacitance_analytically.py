#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 22:05:53 2024

@author: vinzenz
"""

import numpy as np
from scipy.special import ellipk
import matplotlib.pyplot as plt


epsilon = 8.854e-12 #farads per meter

def p(r,R) -> float:
    return (r**2 + R**2)/(2*r*R)-1

def prefactor(r,R) -> float:   
    return 1/(2*np.pi**2 *epsilon* np.sqrt(2*r*R*p(r,R)))

def m(r,R) -> float:
    return -2/p(r,R)

def voltage (r,R):
    return prefactor(r,R) * ellipk(m(r,R))

def voltageofboth(r, R1, R2):
    return voltage(r,R1) - voltage(r,R2)

def capacitance(R1, R2):
    return 1/voltage(R1,R2)*1e12 #in picofarad


R1 = 10*1e-3
R2 = 20*1e-3

r = np.linspace(10.1e-3 , 20e-3 , 199)



plt.figure()
plt.plot(R1+r,voltageofboth(r, R1, R2))
