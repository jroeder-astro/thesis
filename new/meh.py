import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

x1 = []
y1 = []

with open('mr.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

plt.plot(x1, y1, label='Euler')

plt.plot([15.5],[1.0005], 'ro', label='Stefan')
plt.plot([15.62],[0.993025], 'go', label='Jan')
#with open('RK4th.dat', 'r') as csvfile:
#    plots = csv.reader(csvfile, delimiter=' ')
#    for row in plots:
#        x2.append(float(row[0]))
#        y2.append(float(row[1]))

#plt.subplot(3,1,1)
#plt.plot(x1, y1)
#plt.plot(x2, y2, label='RK4')
#plt.axis([0.00007,0.000076,0.0006,0.001])
#plt.title('TOV equation')
plt.ylabel('M')
#plt.legend()

#plt.subplot(3,1,2)
#plt.plot(p, (p/10.)**(3./5.), label='Initial EoS')
#plt.axis([0,0.1,0,0.07])
#plt.ylabel('e(p)')
#plt.legend()

#plt.subplot(3,1,3)
#plt.plot(p, (p/10.)**(3./5.), label='Initial EoS')
#plt.plot(x, y, label='RK4 reconstr.')
#plt.axis([0,0.1,0,0.07])
#plt.ylabel('e(p)')
#plt.legend()


#plt.title('TOV equation\nReconstructed EoS')
#plt.axis([0, 0.07, 0, 11])
#plt.ylabel('e(p)')
plt.xlabel('R')
plt.legend()
plt.show()