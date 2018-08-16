import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *

x1 = []
y1 = []
<<<<<<< HEAD
x2 = []
y2 = []
x3 = []
y3 = []
=======
#x2 = []
#y2 = []
>>>>>>> 226f81eceedb97273e306e32a84fd58efc62e375

def eos(p):
    return np.power(p/10., 3./5.)

def line(p):
    return m*p+b;
<<<<<<< HEAD
=======

p1 = np.arange(0., 0.0004, 0.0000001)
>>>>>>> 226f81eceedb97273e306e32a84fd58efc62e375

p1 = np.arange(0., 0.0004, 0.0000001)

#plt.plot(p1, eos(p1))

K = np.power(10., -5.)

#x = [3.*K, 3.1*K, 4.2*K, 5.3*K, 6.4*K, 7.5*K, 8.6*K, 9.7*K, 0.000108, 0.000119]
#y = [eos(3.*np.power(10.,-5.)), 0.00048661, 0.000619172, 0.000639204, 0.00065824, 0.000767281, 0.000693327, 0.000709379, 0.00072444, 0.00073851]
#plt.plot(x, y, 'r--')


## new ##

#m = [2.1*K, 3.2*K, 4.3*K, 5.4*K, 0.000321]
#n = [0.000381743, 0.000473831, 0.000484831, 0.000495831, 0.00244383]
#plt.plot(m,n,'g--')

## well then ##

<<<<<<< HEAD
=======
## well then ##

>>>>>>> 226f81eceedb97273e306e32a84fd58efc62e375
#p = [3.1*K, 4.2*K, 5.3*K, 6.4*K, ]
#q = []
#plt.plot(p,q,'gr--')




#plt.plot([3.*np.power(10., -5.)],[0.000485593],'go')

#plt.plot([3.1*np.power(10., -5.)],[0.000495502],'bo')

















<<<<<<< HEAD
with open('MR.dat', 'r') as csvfile:
=======
with open('log', 'r') as csvfile:
>>>>>>> 226f81eceedb97273e306e32a84fd58efc62e375
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

<<<<<<< HEAD
plt.plot(x1,y1,'r', label='MR.dat')

with open('mrcurve', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))

plt.plot(x2,y2,'g', label='mrcurve')

with open('mr_rel_rk4.out', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x3.append(float(row[0]))
        y3.append(float(row[1]))

plt.plot(x3,y3,'b', label='mr_rel_rk4')






plt.plot([3.8538*K],[0.000455623],'ro')
=======
plt.plot(x1,y1,'c--')

>>>>>>> 226f81eceedb97273e306e32a84fd58efc62e375

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
plt.ylabel('M(R)')
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
