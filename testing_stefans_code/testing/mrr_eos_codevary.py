import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

def eos(p):
    return np.power(p/(10.), 3./5.)

Rl26 = []
Ml26 = []
el26 = []
pl26 = []

Rl26w1 = []
Ml26w1 = []
el26w1 = []
pl26w1 = []

Rl26w2 = []
Ml26w2 = []
el26w2 = []
pl26w2 = []

Rl26wb = []
Ml26wb = []
el26wb = []
pl26wb = []

Rd = []
Md = []
ed = []
pd = []

with open('results_l2_7', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl26.append(float(row[1]))
        Ml26.append(float(row[3]))
        pl26.append(float(row[4]))
        el26.append(float(row[5]))

with open('results_l2_7_wo1', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl26w1.append(float(row[1]))
        Ml26w1.append(float(row[3]))
        pl26w1.append(float(row[4]))
        el26w1.append(float(row[5]))

with open('results_l2_7_wo2', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl26w2.append(float(row[1]))
        Ml26w2.append(float(row[3]))
        pl26w2.append(float(row[4]))
        el26w2.append(float(row[5]))

with open('results_l2_7_wob', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl26wb.append(float(row[1]))
        Ml26wb.append(float(row[3]))
        pl26wb.append(float(row[4]))
        el26wb.append(float(row[5]))

with open('mr.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rd.append(float(row[0]))
        Md.append(float(row[1]))
        pd.append(float(row[2]))
        ed.append(float(row[3]))

plt.figure(figsize=(14,7))

el26n = [a*P for a in el26]
pl26n = [a*P for a in pl26]

el26w1n = [a*P for a in el26w1]
pl26w1n = [a*P for a in pl26w1]

el26w2n = [a*P for a in el26w2]
pl26w2n = [a*P for a in pl26w2]

el26wbn = [a*P for a in el26wb]
pl26wbn = [a*P for a in pl26wb]

p1 = np.arange(0., 0.00025, 0.00000001)

plt.subplot(2,2,(1,3))
plt.xlim(11.5, 14)
plt.ylim(1.1, 1.2)
plt.plot(Rd, Md, label = 'MRR input')
plt.plot(Rl26, Ml26,label = 'MRR $\lambda = 0.07$')
plt.plot(Rl26w1, Ml26w1,label = 'MRR, w/o p1')
plt.plot(Rl26w2, Ml26w2,label = 'MRR, w/o p2')
plt.plot(Rl26wb, Ml26wb,label = 'MRR, w/o b')
plt.ylabel('M/M$_\odot$', fontsize=15)
plt.xlabel('R/km', fontsize=15)
plt.legend(loc=1,prop={'size':12})

plt.subplot(2,2,(2,4))
plt.plot(p1*P, eos(p1)*P, label = 'EOS input')
plt.plot(pl26n, el26n, label ='EOS $\lambda = 0.07$')
plt.plot(pl26w1n, el26w1n, label ='EOS, w/o p1')
plt.plot(pl26w2n, el26w2n, label ='EOS, w/o p2')
plt.plot(pl26wbn, el26wbn, label ='EOS, w/o b')
plt.ylabel('$\epsilon$(p)/MeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2,prop={'size':12})

plt.tight_layout()
plt.show()

