import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

def eos(p):
    return np.power(p/10., 3./5.)

def inv(e):
    return 10*np.power(e, 5./3.)

p1 = np.arange(0., 0.001, 0.0000001)
e1 = np.arange(0., 0.004, 0.0000001)

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

plt.figure(figsize=(8,4))

plt.subplot(2,2,(1,3))
plt.plot(x1, y1, label='MRR')
plt.plot([10.68],[1.21468], 'go')
plt.ylabel('$M(R)$',fontsize = 20)
plt.xlabel('$R$',fontsize = 20)
plt.yticks([])
plt.xticks([])
plt.legend()

plt.subplot(2,2,(2,4))
plt.plot(x2, y2, label='EOS')
plt.plot([0.0004], [0.0022974], 'go')
plt.ylabel('$e(p)$',fontsize = 20)
plt.xlabel('$p$',fontsize = 20)
plt.yticks([])
plt.xticks([])
plt.legend(loc=4)

plt.show()
