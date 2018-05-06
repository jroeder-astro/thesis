import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

x = []
y = []

#p = np.arange(0., 50., 0.01)

with open('one.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))

plt.plot(x, y, label='4th order\nRunge-Kutta')
#plt.plot(p, (p/10.)**(3./5.), label='EoS')

plt.title('TOV equation')
#plt.axis([0, 0.07, 0, 11])
plt.ylabel('E(r)')
plt.xlabel('r')
plt.legend()
plt.show()
