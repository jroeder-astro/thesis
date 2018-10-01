import matplotlib.pyplot as plt
import numpy as np
import csv
def eos(p):
    return np.power(p/(10.), 3./5.)

x1 = []
y1 = []

with open('twin.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

#x2 = [a*P for a in x1]
#y2 = [b*P for b in y1]

plt.plot(x1, y1, 'bo', label='MRR')

#p1 = np.arange(0., 0.0005, 0.00000001)
#plt.plot(p1*P, eos(p1)*P, label='known')

plt.title('Twin Star MRR/EOS')
plt.ylabel('$M/M_\odot$')
plt.xlabel('$R/km$')
plt.legend()
plt.show()
