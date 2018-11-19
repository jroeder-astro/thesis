import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

def eos(p):
    return np.power(p/(10.), 3./5.)

Rl2 = []
Ml2 = []
el2 = []
pl2 = []

Rl3 = []
Ml3 = []
el3 = []
pl3 = []

Rd = []
Md = []
ed = []
pd = []

with open('results_l3', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl3.append(float(row[1]))
        Ml3.append(float(row[3]))
        pl3.append(float(row[4]))
        el3.append(float(row[5]))

with open('results_l2', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl2.append(float(row[1]))
        Ml2.append(float(row[3]))
        pl2.append(float(row[4]))
        el2.append(float(row[5]))

with open('mr.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rd.append(float(row[0]))
        Md.append(float(row[1]))
        pd.append(float(row[2]))
        ed.append(float(row[3]))

plt.figure(figsize=(8,4))

el3n = [a*P for a in el3]
pl3n = [a*P for a in pl3]
el2n = [a*P for a in el2]
pl2n = [a*P for a in pl2]

p1 = np.arange(0., 0.00018, 0.00000001)

plt.subplot(2,2,(1,3))
plt.plot(Rd, Md, label = 'MRR input')
plt.plot(Rl2, Ml2,label = 'MRR $\lambda = 0.01$')
plt.plot(Rl3, Ml3,label = 'MRR $\lambda = 0.001$')
plt.ylabel('M/M$_\odot$', fontsize=15)
plt.xlabel('R/km', fontsize=15)
plt.legend(loc=3,prop={'size':14})

plt.subplot(2,2,(2,4))
plt.plot(p1*P, eos(p1)*P, label = 'EOS input')
plt.plot(pl2n, el2n, label ='EOS $\lambda = 0.01$')
plt.plot(pl3n, el3n, label ='EOS $\lambda = 0.001$')
plt.ylabel('$\epsilon$(p)/GeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/GeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2,prop={'size':14})

plt.tight_layout()
plt.show()

