#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create Poincare plot
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


file_name = "../output.csv"
start = 0 #offset for plotting; choose poloidal plane.
stride = 1000 #if this is N, you plot every N points, so these lie in one poloidal plane.
num_tracers = 10

positions = []

with open(file_name, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in reader:
        positions.append(float(row[0]))
        
positions = np.array(positions)

major_radius = positions[::3]
phi = positions[1::3]
vertical_position = positions[2::3]
size = len(major_radius)
new_size = size//num_tracers

#unpack data
major_radius = np.reshape(major_radius, (num_tracers,new_size))
phi = np.reshape(phi, (num_tracers,new_size))
vertical_position = np.reshape(vertical_position, (num_tracers,new_size))

#poincare plot
plt.figure()
plt.scatter(major_radius[:,start::stride],vertical_position[:,start::stride],c='k')
plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.axis('equal')

# #3D plot of a field line, to show the helical nature of it
# ax = plt.figure().add_subplot(projection='3d')
# X = major_radius*np.cos(phi)
# Y = major_radius*np.sin(phi)
# Z = vertical_position
# ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z))) 
# ax.plot(X,Y,Z)
# plt.show()