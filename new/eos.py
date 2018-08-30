import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

K = np.power(10., -5.)

def eos(p):
    return np.power(p/10., 3./5.)

p1 = np.arange(0., 0.00032, 0.00001)
plt.plot(p1, eos(p1), label='known')

x1 = []
y1 = []

with open('log', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[3]))
        y1.append(float(row[2]))

plt.plot(x1, y1, 'bo', label='recon')

plt.title('Equation of State\nReconstruction Algorithms')
plt.ylabel('e(p)')
plt.xlabel('p')
plt.legend(loc=4)
plt.show()
