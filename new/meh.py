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

with open('rec_stefan.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))

#plt.plot(x2, y2, 'ro', label='Rec_Stefan')

with open('log_new', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x3.append(float(row[0]))
        y3.append(float(row[1]))

plt.plot(x3, y3, 'go', label='Rec_Stefan')

with open('rec_stefan_rk4.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x4.append(float(row[0]))
        y4.append(float(row[1]))

#plt.plot(x4, y4, 'co', label='Rec_Stefan_RK4')


#plt.subplot(2,1,1)
#plt.plot(x1, y1, label='Euler')
#plt.plot(x2, y2, 'ro', label='Rec_Stefan')
#plt.plot(x3, y3, 'go', label='Rec_Jan')
#plt.ylabel('M')
#plt.legend()

#plt.subplot(2,1,2)
#plt.plot(x1, y1, label='Euler')
#plt.plot(x4, y4, 'ro', label='Rec_Stefan_RK4')
#plt.plot(x3, y3, 'go', label='Rec_Jan')
#plt.ylabel('M')
#plt.legend()


plt.title('Mass-Radius-Relation\nReconstruction Algorithms')
plt.ylabel('M(R)')
plt.xlabel('R')
plt.legend()
plt.show()
