import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

def eos(p):
    return np.power(p/10., 3./5.)

p1 = np.arange(0., 0.00012, 0.0000001)

plt.plot(p1, eos(p1))

K = np.power(10., -5.)

# Me

x1 = [3.*K, 3.1*K, 4.2*K, 5.3*K, 6.4*K, 7.5*K, 8.6*K, 9.7*K, 0.000108, 0.000119]
y1 = [eos(3.*np.power(10.,-5.)), 0.00048661, 0.000619172, 0.000639204, 0.00065824, 0.000767281, 0.000693327, 0.000709379, 0.00072444, 0.00073851]
plt.plot(x1, y1, 'r-', label='Me')

# Stefan

x2 = [4.3184*K, 4.9191*K, 5.8289*K]
y2 = [0.000202575, 0.000210364, 0.000217736]
plt.plot(x2, y2, 'g-', label='Stefan')

# Me.2

x3 = [2.1*K, 3.2*K, 4.3*K, 5.4*K]
y3 = [0.000381743, 0.00047331, 0.000493342, 0.000512378]
plt.plot(x3, y3, 'b-', label='MeAsw')


#

plt.plot([3.8538*K],[0.000455623],'ro')

#


#plt.plot([3.*np.power(10., -5.)],[0.000485593],'go')

#plt.plot([3.1*np.power(10., -5.)],[0.000495502],'bo')

















#with open('output.dat', 'r') as csvfile:
#    plots = csv.reader(csvfile, delimiter=',')
#    for row in plots:
#        x1.append(float(row[0]))
#        y1.append(float(row[1]))

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
plt.ylabel('e')
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
plt.xlabel('p')
plt.legend()
plt.show()
