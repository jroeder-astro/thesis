import matplotlib.pyplot as plt
import numpy as np
import csv

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

plt.subplot(2,2,(1,3))
plt.plot(R, M, label = 'MRR')
plt.ylabel('Mass', fontsize=20)
plt.xlabel('Radius', fontsize=20)
plt.legend()

plt.subplot(2,2,(2,4))
plt.plot(p, e, label ='EOS')
plt.ylabel('e(p)', fontsize=20)
plt.xlabel('p', fontsize=20)
plt.legend()

plt.tight_layout()
plt.show()

