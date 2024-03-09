#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 13:40:56 2024

@author: vinzenz
"""
import numpy as np
import matplotlib.pyplot as plt

"""
time to load the capacitor t = 0.69CR
hence the period is T = 2*t
measuring the frequency f = 1/T = 1/2t
we get t=1/2f
hence the capacitance is C = 1/(2*0.69Rf)
"""

#correcting for capacitance of the circuit itself

# Ohne Sensor, 0pF, 100 kO
#517.8, 97.99, 103
C_circuit_100k = 1/(2*0.69*100e3*517.8e3)*1e12 #pico Farad

# Ohne Sensor, 100 pF, 100 kO
#61.35, 55.97, 119
C_100pF_100k = 1/(2*0.69*100e3*61.35e3)*1e12

#with capacitance of the circuit being in series
C_100pF_100k_corrected = C_100pF_100k-C_circuit_100k
print('corrected measured capacitance, 100pF, 100kOhm: {}'.format(C_100pF_100k_corrected))



# Ohne Sensor, 100pF, 6.8kO
#833.5, 0, 90
C_100pF_6d8k = 1/(2*0.69*6.8e3*883.5e3)*1e12

# Ohne Sensor, 0pF, 6.8kO, Square wave destroyed (stil visible), bad sampling rate
#4500, 4349, 104
C_0pF_6d8k = 1/(2*0.69*6.8e3*4500e3)*1e12

#correcting the capacitance
C_100pF_6d8k_corrected = C_100pF_6d8k - C_0pF_6d8k
print('corrected measured capacitance, 100pF, 6.8kOhm: {}'.format(C_100pF_6d8k_corrected))



capacitance_correction = C_circuit_100k

#sensor A
# Abstand [mm], Mittelwert [kHz], std [Hz], n
A = np.array([
[10.0, 343.2, 0, 107],
[5.0, 333.0, 16.30, 104],
[3.0, 319.1, 0, 104],
[2.0, 303.3, 0, 105],
[1.0, 265.6, 0, 101],
[0.5, 216.3, 105.3, 104]])
A_open = np.array([1000, 363.7, 317.6, 107])
A_C_open = 1/(2*0.69*100e3*A_open[1]*1e3)*1e12 - capacitance_correction
print('corrected measured capacitance, A open, 100kOhm: {}'.format(A_C_open))

A_d = A[:,0]
A_C = 1/(2*0.69*100e3*A[:,1]*1e3)*1e12 - capacitance_correction

B = np.array([
[10.0,400.7, 0, 103],
[5.0, 394.0, 488.8, 105],
[3.0, 383.5, 0, 103],
[2.0, 371.9, 283.3, 106],
[1.0, 341.7, 225.0, 100],
[0.5, 308.9, 221.9, 102]])
B_open = np.array([1000, 423.3, 658.3, 113])
B_C_open = 1/(2*0.69*100e3*B_open[1]*1e3)*1e12 - capacitance_correction
print('corrected measured capacitance, B open, 100kOhm: {}'.format(B_C_open))


B_d = B[:,0]
B_C = 1/(2*0.69*100e3*B[:,1]*1e3)*1e12 - capacitance_correction


C = np.array([
[10.0, 423.4, 0, 104],
[5.0, 422.3, 120.4, 97],
[3.0, 417.7, 394.8, 104],
[2.0, 414.0, 0, 102],
[1.0, 400.9, 430.5, 103],
[0.5, 386.1, 267.6, 105]])
C_open = np.array([1000, 444.8, 431.2, 105])
C_C_open = 1/(2*0.69*100e3*C_open[1]*1e3)*1e12 - capacitance_correction
print('corrected measured capacitance, C open, 100kOhm: {}'.format(C_C_open))


C_d = C[:,0]
C_C = 1/(2*0.69*100e3*C[:,1]*1e3)*1e12 - capacitance_correction

D = np.array([
[10.0, 440.1, 28.44, 97],
[5.0, 438.9, 419.9, 105],
[3.0, 437.2, 0, 107],
[2.0, 435.5, 0, 105],
[1.0, 428.8, 152.3, 108],
[0.5, 422.6, 283.1, 98]])
D_open = np.array([1000, 459.9, 232.8, 105])
D_C_open = 1/(2*0.69*100e3*D_open[1]*1e3)*1e12 - capacitance_correction
print('corrected measured capacitance, D open, 100kOhm: {}'.format(D_C_open))


D_d = D[:,0]
D_C = 1/(2*0.69*100e3*D[:,1]*1e3)*1e12 - capacitance_correction

fig = plt.figure()
plt.plot(A_d, A_C)
plt.plot(B_d, B_C)
plt.plot(C_d, C_C)
plt.plot(D_d, D_C)

"""
dist = B_C
dist = dist[::-1]
print(np.array2string(dist, precision=2, separator=', '))
"""
print(D_C_open)