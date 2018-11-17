import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

def eos(p):
    return np.power(p/(10.), 3./5.)

x1 = []
y1 = []

with open('recon_me_fail.dat', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[3]))
        y1.append(float(row[2]))

x2 = [a*P for a in x1]
y2 = [b*P for b in y1]

p1 = np.arange(0., 0.00025, 0.00000001)
plt.plot(p1*P, eos(p1)*P, label='known EOS')
plt.plot(x2,y2,'bo',label='reconstruction')

plt.ylabel('$\epsilon$(p)/MeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2, prop={'size':14})
plt.show()
