import matplotlib.pyplot as plt
import numpy as np
import csv

#from numpy import *

K = np.power(10., -5.)
P = 1.324*np.power(10., -6.)
Q = np.power(P, -2./3.)

def eos(p):
    return np.power(p*P/(Q*10.), 3./5.)

x1 = []
y1 = []

with open('log_save', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[3]))
        y1.append(float(row[2]))

x2 = [a*P for a in x1]
y2 = [b*P for b in y1]

plt.plot(x2, y2, 'bo', label='recon')

p1 = np.arange(0., 0.0005, 0.00000001)
plt.plot(p1*P, eos(p1), label='known')

plt.title('Equation of State\nReconstruction Algorithms')
plt.ylabel('$e(p)$ /MeV/fm$^3$')
plt.xlabel('$p$ /MeV/fm$^3$')
plt.legend(loc=4)
plt.show()
