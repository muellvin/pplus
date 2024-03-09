#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 22:05:53 2024

@author: vinzenz
"""

import numpy as np
from scipy.special import ellipk


epsilon = 8.854e-12 #farads per meter

def p(r,R) -> float:
    return (r**2 + R**2)/(2*r*R)-R 

def prefactor(r,R) -> float:   
    return 1/(2*np.pi**2 *epsilon* np.sqrt(2*r*R*p(r,R)))

def m(r,R) -> float:
    return -2/p(r,R)

def voltage (r,R):
    return prefactor(r,R) * ellipk(m(r,R))

def voltagediff(R1, R2):
    voltageatR1 = voltage(R1,R1) - voltage(R1,R2)
    voltageatR2 = voltage(R2,R1) - voltage(R2,R2)
    return voltageatR1 - voltageatR2

def capacitance(R1, R2):
    return 1/voltagediff(R1,R2)*1e12 #in picofarad


R1 = np.array([5, 10, 20, 30])*1e-3

s = 5e-3

R2 = R1+s
capacitances=capacitance(R1,R2)

print(capacitances)