import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

R = []
M = []
e = []
p = []

with open('mr.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        R.append(float(row[0]))
        M.append(float(row[1]))
        e.append(float(row[2]))
        p.append(float(row[3]))

plt.figure(figsize=(8,4))

en = [a*P for a in e]
pn = [a*P for a in p]

eG = [a/1000 for a in en]
pG = [a/1000 for a in pn]

plt.subplot(2,2,(1,3))
plt.plot(R, M, label = 'MRR')
plt.ylabel('M/M$_\odot$', fontsize=15)
plt.xlabel('R/km', fontsize=15)
plt.legend(loc=2,prop={'size':14})

plt.subplot(2,2,(2,4))
plt.plot(pG, eG, label ='EOS')
plt.ylabel('$\epsilon$(p)/GeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/GeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2,prop={'size':14})

plt.tight_layout()
plt.show()

