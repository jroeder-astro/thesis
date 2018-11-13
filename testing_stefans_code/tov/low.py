import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

e = []
p = []

with open('lowdensity.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        e.append(float(row[1]))
        p.append(float(row[0]))

pn = [a*P for a in p]
en = [a*P for a in e]

plt.plot(pn, en, label ='EOS')
plt.ylabel('$\epsilon$(p)/MeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=15)
plt.legend(prop={'size':14})

plt.tight_layout()
plt.show()

