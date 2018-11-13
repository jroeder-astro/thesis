import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

def eos(p):
    return np.power(p/(10.), 3./5.)

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

x3 = [a*P for a in x2]
y3 = [b*P for b in y2]

plt.figure(figsize=(8,4))

plt.subplot(2,2,(1,3))
plt.plot(x1, y1, label='MRR')
plt.ylabel('M(R)/M$_\odot$', fontsize=15)
plt.xlabel('R/km', fontsize=15)
plt.legend(prop={'size':14})

plt.subplot(2,2,(2,4))
plt.plot(x3, y3, label='EOS')
plt.ylabel('$\epsilon$(p)/MeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2, prop={'size':14})

plt.tight_layout()
plt.show()
