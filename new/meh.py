import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

x1 = []
y1 = []

x2 = []
y2 = []

x3 = []
y3 = []

x4 = []
y4 = []

with open('mr.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

plt.plot(x1, y1, label='Euler')

#plt.plot([15.5, 14.99],[1.0005, 1.03217], 'ro', label='Stefan loop')
#plt.plot([15.59],[0.994556], 'bo', label='Jan loop')

with open('rec_stefan.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))

plt.plot(x2, y2, 'ro', label='Rec_Stefan')

with open('rec_jan.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x3.append(float(row[0]))
        y3.append(float(row[1]))

plt.plot(x3, y3, 'go', label='Rec_Jan')

with open('rec_stefan_rk4.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x4.append(float(row[0]))
        y4.append(float(row[1]))

plt.plot(x2, y2, 'co', label='Rec_Stefan_RK4')


#plt.plot([15.5, 14.99, 14.59],[1.00014, 1.03397, 1.06196], 'go', label='Rec_Jan')



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
