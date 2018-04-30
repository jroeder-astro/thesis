import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

x = []
y = []

with open('tov.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))

#p1 = polyfit(x, y, 1)
#print(p1)

plt.plot(x, y, label='4th order\nRunge-Kutta')
#plt.plot(x, polyval(p1,x))

plt.title('TOV equation\nAttempt for numerical solution')
#plt.axis([0, 0.07, 0, 11])
plt.ylabel('P(r)')
plt.xlabel('r')
plt.legend()
plt.show()
