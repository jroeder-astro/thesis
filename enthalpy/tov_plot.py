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

x2 = []
y2 = []

with open('mr2.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))

plt.plot(x2, y2, label='add')


plt.title('TOV equation')
plt.ylabel('M')
plt.xlabel('R')
plt.legend()
plt.show()
