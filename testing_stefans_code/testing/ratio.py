import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

def eos(p):
    return np.power(p/(10.), 3./5.)

def inv(e):
    return 10.*np.power(e, 5./3.)

Rl2 = []
Ml2 = []
el2 = []
pl2 = []

Rl25 = []
Ml25 = []
el25 = []
pl25 = []

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

with open('results_l2_5', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl25.append(float(row[1]))
        Ml25.append(float(row[3]))
        pl25.append(float(row[4]))
        el25.append(float(row[5]))

with open('mr.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rd.append(float(row[0]))
        Md.append(float(row[1]))
        pd.append(float(row[2]))
        ed.append(float(row[3]))

plt.figure(figsize=(10,5))

el3n = [a*P for a in el3]
pl3n = [a*P for a in pl3]
el2n = [a*P for a in el2]
pl2n = [a*P for a in pl2]
el25n = [a*P for a in el25]
pl25n = [a*P for a in pl25]

edn =  [P*eos(a) for a in pl25]
eRatL25 = []
endind = len(edn)
for i in range(endind):
    eRatL25.append(el25n[i]/edn[i])


#eRatL3 = [a/eos(b) for a in el3n and b in pl3n]
#epInvL3 = [inv(a) for a in el3n]
#pRatL3 = [a/b for a in pl3n and b in epInvL3]


p1 = np.arange(0., 0.00025, 0.00000001)
oneX = [0, 180]
oneY = [1, 1]

#plt.subplot(2,2,(1,3))
#plt.plot(Rd, Md, label = 'MRR input')
#plt.plot(Rl2, Ml2,label = 'MRR $\lambda = 0.01$')
#plt.plot(Rl25, Ml25,label = 'MRR $\lambda = 0.05$')
#plt.plot(Rl3, Ml3,label = 'MRR $\lambda = 0.001$')
#plt.ylabel('M/M$_\odot$', fontsize=15)
#plt.xlabel('R/km', fontsize=15)
#plt.legend(loc=3,prop={'size':12})

#plt.subplot(2,2,(2,4))
#plt.plot(p1*P, eos(p1)*P, label = 'EOS input')
#plt.plot(pl2n, el2n, label ='EOS $\lambda = 0.01$')
#plt.plot(pl25n, el25n, label ='EOS $\lambda = 0.05$')

plt.plot(pl25n, eRatL25, label ='EOS ratio $\lambda = 0.05$')
plt.plot(oneX, oneY, 'b--')
plt.ylabel('$\epsilon$(p)/$\epsilon_{default}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2,prop={'size':12})

plt.tight_layout()
plt.show()

