import matplotlib.pyplot as plt
import numpy as np
import csv

x1 = []
y1 = []

with open('stuff.dat', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

plt.plot(x1, y1)

plt.show()
