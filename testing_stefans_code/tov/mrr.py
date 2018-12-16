import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

Rd = []
Md = []
ed = []
pd = []

with open('mr_eos_com.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        Rd.append(float(row[0]))
        Md.append(float(row[1]))
        ed.append(float(row[2]))
        pd.append(float(row[3]))

edn = [a*P for a in ed]
pdn = [a*P for a in pd]

R25 = []
M25 = []
e25 = []
p25 = []

with open('results_el1_1', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R25.append(float(row[1]))
        M25.append(float(row[3]))
        e25.append(float(row[5]))
        p25.append(float(row[4]))

e25n = [a*P for a in e25]
p25n = [a*P for a in p25]

R255 = []
M255 = []
e255 = []
p255 = []

with open('results_el2_55', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R255.append(float(row[1]))
        M255.append(float(row[3]))
        e255.append(float(row[5]))
        p255.append(float(row[4]))

e255n = [a*P for a in e255]
p255n = [a*P for a in p255]

R257 = []
M257 = []
e257 = []
p257 = []

with open('results_el2_57', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R257.append(float(row[1]))
        M257.append(float(row[3]))
        e257.append(float(row[5]))
        p257.append(float(row[4]))

e257n = [a*P for a in e257]
p257n = [a*P for a in p257]


plt.figure(figsize=(8,4))

plt.subplot(2,2,(1,3))
#plt.xlim(10, 13)
#plt.ylim(1.1, 1.2)
plt.plot(Rd, Md, label = 'MRR')
#plt.plot(R25, M25, label = '0.05')
#plt.plot(R255, M255, label = '0.055')
#plt.plot(R257, M257, label = '0.057')
plt.ylabel('M/M$_\odot$', fontsize=15)
plt.xlabel('R/km', fontsize=15)
plt.legend(loc=2,prop={'size':14})

plt.subplot(2,2,(2,4))
#plt.xlim(0, 200)
plt.plot(pdn, edn, label ='EOS')
#plt.plot(p25n, e25n, label ='0.05')
#plt.plot(p255n, e255n, label ='0.055')
#plt.plot(p257n, e257n, label ='0.057')
plt.ylabel('$\epsilon$(p)/MeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2,prop={'size':14})

plt.tight_layout()
plt.show()

