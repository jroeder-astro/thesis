import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

x1 = []
y1 = []

x2 = []
y2 = []
    
with open('tov_rk4.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

with open('tov_euler.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))

plt.plot(x1, y1, label='4th order\nRunge-Kutta')
plt.plot(x2, y2, label='1st order\nEuler')

plt.title('M-R-Relation\n1st vs. 4th order')
#plt.axis([0, 0.07, 0, 11])
plt.ylabel('M')
plt.xlabel('R')
plt.legend()
plt.show()
