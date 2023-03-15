#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create Poincare plot
"""

import csv
import numpy as np
import matplotlib.pyplot as plt


file_name = "../output.csv"
start = 0 #offset for plotting; choose poloidal plane.
stride = 1000 #if this is N, you plot every N points, so these lie in one poloidal plane.

positions = []

with open(file_name, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in reader:
        positions.append(float(row[0]))
        
positions = np.array(positions)

R = positions[::2]
Z = positions[1::2]

plt.figure()
plt.scatter(R[start::stride],Z[start::stride],c='k')
plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.axis('equal')
plt.show()