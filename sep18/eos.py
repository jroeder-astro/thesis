import matplotlib.pyplot as plt
import numpy as np
import csv

#from numpy import *

K = np.power(10., -5.)
P = 7.55616208*np.power(10., 5.)
Q = np.power(P, -2./3.)

def eos(p):
    return np.power(p/(10.), 3./5.)

x1 = []
y1 = []

with open('logn', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[3]))
        y1.append(float(row[2]))

x2 = [a*P for a in x1]
y2 = [b*P for b in y1]

plt.plot(x2, y2, 'bo', label='recon')

p1 = np.arange(0., 0.00025, 0.00000001)
plt.plot(p1*P, eos(p1)*P, label='known')

plt.title('Equation of State\nReconstruction Algorithms')
plt.ylabel('$\epsilon(p)$ /MeV/fm$^3$')
plt.xlabel('$p$ /MeV/fm$^3$')
plt.legend(loc=4)
plt.show()
