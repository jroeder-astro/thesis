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

x2 = []
y2 = []

with open('out.dat', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))
        y2.append(float(row[2]))
        x2.append(float(row[3]))

x3 = [] 
y3 = []

x3 = [a*P for a in x2]
y3 = [b*P for b in y2]

#plt.plot(x2, y2, 'bo', label='recon')

#p1 = np.arange(0., 0.0005, 0.00000001)
#plt.plot(p1*P, eos(p1), label='known')

plt.figure(figsize=(11,5))

plt.subplot(2,2,(1,3))
plt.plot(x1, y1, label='MRR')
plt.ylabel('$M(R)$ /M$_\odot$', fontsize=20)
plt.xlabel('$R$ /km', fontsize=20)
plt.legend()

plt.subplot(2,2,(2,4))
plt.plot(x3, y3, label='EOS')
plt.ylabel('$\epsilon(p)$ /MeV/fm$^3$', fontsize=20)
plt.xlabel('$p$ /MeV/fm$^3$', fontsize=20)
plt.legend(loc=4)

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
plt.show()
