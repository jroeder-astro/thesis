import matplotlib.pyplot as plt
import numpy as np
import csv

#from numpy import *

P = 7.55616208*np.power(10., 5.) 

def eos(p):
    return np.power(p/(10.), 3./5.)

p1 = np.arange(0., 25*np.power(10.,-5.), 0.000001)

x1 = [] # Radius calc
y1 = [] # Mass calc

x2 = [] # p wrong
y2 = [] # e wrong

x4 = [] # Radius orig
y4 = [] # Mass orig

with open('results', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=' ')
    for row in plots:
        x1.append(float(row[1]))
        y1.append(float(row[3]))
        x2.append(float(row[4]))
        y2.append(float(row[5]))
        x4.append(float(row[0]))
        y4.append(float(row[2]))

x3 = [a*P for a in x2] # p right
y3 = [b*P for b in y2] # e right

#plt.plot(x2, y2, 'bo', label='recon')

#p1 = np.arange(0., 0.0005, 0.00000001)
#plt.plot(p1*P, eos(p1), label='known')

plt.figure(figsize=(11,5))

plt.subplot(2,2,(1,3))
plt.plot(x1, y1, label='MRR reconstructed')
plt.plot(x4, y4, label='MRR input')
plt.ylabel('$M(R)$ /M$_\odot$', fontsize=20)
plt.xlabel('$R$ /km', fontsize=20)
plt.legend()

plt.subplot(2,2,(2,4))
plt.plot(x3, y3, label='EOS reconstructed')
plt.plot(p1*P, eos(p1)*P, label='EOS used for MRR')
plt.ylabel('$\epsilon(p)$ /MeV/fm$^3$', fontsize=20)
plt.xlabel('$p$ /MeV/fm$^3$', fontsize=20)
plt.legend(loc=2)

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
plt.show()
