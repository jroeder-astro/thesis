import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

R25 = []
M25 = []
e25 = []
p25 = []

with open('mr_eos_com.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        R25.append(float(row[0]))
        M25.append(float(row[1]))
        e25.append(float(row[2]))
        p25.append(float(row[3]))

e25n = [a*P for a in e25]
p25n = [a*P for a in p25]

R28 = []
M28 = []
e28 = []
p28 = []

with open('results', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R28.append(float(row[1]))
        M28.append(float(row[3]))
        e28.append(float(row[5]))
        p28.append(float(row[4]))

e28n = [a*P for a in e28]
p28n = [a*P for a in p28]

plt.figure(figsize=(8,4))

plt.subplot(2,2,(1,3))
plt.xlim(10, 13)
#plt.ylim(1.1, 1.2)
plt.plot(R25, M25, label = 'MRR')
#plt.plot(R28, M28, label = 'recon')
plt.ylabel('M/M$_\odot$', fontsize=15)
plt.xlabel('R/km', fontsize=15)
plt.legend(prop={'size':14})

plt.subplot(2,2,(2,4))
#plt.xlim(0, 200)
plt.plot(p25n, e25n, label ='EOS')
#plt.plot(p28n, e28n, label ='recon')
plt.ylabel('$\epsilon$(p)/MeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=15)
plt.legend(loc=2,prop={'size':14})

plt.tight_layout()
plt.show()

