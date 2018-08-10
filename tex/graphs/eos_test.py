import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *


def eos(p):
    return np.power(p/10., 3./5.)

def deriv(p):
    return 3./5. * np.power(10., -3./5.) * np.power(p, -2./5.)

def cause(p):
    return p


p1 = np.arange(0., 1.*np.power(10., -5.), 0.0000001)
p2 = np.arange(0., 2.*np.power(10., -5.), 0.0000001)


plt.subplot(1,3,1)
plt.plot(p1, eos(p1), label='given EOS')
plt.plot(p2, cause(p2), 'r--')
plt.ylabel('$\epsilon(p)$')
plt.xlabel('$p$')
plt.yticks([])
plt.xticks([])
plt.title('a)')
plt.axis([0,2.*np.power(10.,-5.), 0,  eos(2.*np.power(10.,-5.))+np.power(10.,-6.) ] )
plt.legend(loc = 2)


plt.subplot(1,3,2)
plt.plot(p1, eos(p1), label='given EOS')
plt.plot([1.*np.power(10.,-5.), 2.*np.power(10.,-5.)],[eos(1.*np.power(10.,-5.)), eos(2.*np.power(10.,-5.))-0.00005], label='added line')
plt.yticks([])
plt.xticks([])
plt.axis([0,2.*np.power(10.,-5.), 0,  eos(2.*np.power(10.,-5.))+np.power(10.,-6.) ] )
plt.title('b)')
plt.legend(loc=2)
plt.xlabel('$p$')
plt.ylabel('$\epsilon(p)$')


plt.subplot(1,3,3)
plt.plot(p1, eos(p1), label='given EOS')
plt.plot([1.*np.power(10.,-5.), 1.5*np.power(10.,-5.)],[eos(1.*np.power(10.,-5.)), eos(1.5*np.power(10.,-5.))-0.00005], label='1st line')
plt.plot([1.5*np.power(10.,-5.), 2.*np.power(10.,-5.)],[eos(1.5*np.power(10.,-5.))-0.00005, eos(2.*np.power(10.,-5.))-0.00003], label='2nd line')
plt.yticks([])
plt.xticks([])
plt.axis([0,2.*np.power(10.,-5.), 0,  eos(2.*np.power(10.,-5.))+np.power(10.,-6.) ] )
plt.title('c)')
plt.legend(loc=2)
plt.xlabel('$p$')
plt.ylabel('$\epsilon(p)$')


plt.show()