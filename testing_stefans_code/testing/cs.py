import matplotlib.pyplot as plt
import numpy as np
import csv

cs = []
ld = []

with open('cs.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter=' ')
    for row in plots:
        cs.append(float(row[0]))
        ld.append(float(row[1]))

plt.plot(ld, cs, 'bo')

plt.ylabel('$\chi$$^2$', fontsize=15)
plt.xlabel('$\epsilon$', fontsize=15)
plt.legend(loc=2, prop={'size':14})
plt.tight_layout()
plt.show()
