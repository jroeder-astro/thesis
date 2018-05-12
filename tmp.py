import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

x = []
y = []

p = np.arange(0., 0.1, 0.0001)

with open('inv.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))

plt.subplot(3,1,1)
plt.plot(x, y, label='RK4 reconstr.')
plt.axis([0,0.1,0,0.07])
plt.title('TOV equation\nReconstructed EoS')
plt.ylabel('e(p)')
plt.legend()

plt.subplot(3,1,2)
plt.plot(p, (p/10.)**(3./5.), label='Initial EoS')
plt.axis([0,0.1,0,0.07])
plt.ylabel('e(p)')
plt.legend()

plt.subplot(3,1,3)
plt.plot(p, (p/10.)**(3./5.), label='Initial EoS')
plt.plot(x, y, label='RK4 reconstr.')
plt.axis([0,0.1,0,0.07])
plt.ylabel('e(p)')
plt.legend()


#plt.title('TOV equation\nReconstructed EoS')
#plt.axis([0, 0.07, 0, 11])
#plt.ylabel('e(p)')
plt.xlabel('p')
#plt.legend()
plt.show()
