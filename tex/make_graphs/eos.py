import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import *


def eos(p):
    return np.power(p/10., 3./5.)

def deriv(p):
    return 3./5. * np.power(10., -3./5.) * np.power(p, -2./5.)

p1 = np.arange(0., 1.*np.power(10., -5.), 0.0000001)

plt.figure(figsize=(12, 5))

plt.subplot(1,3,1)
plt.plot(p1, eos(p1), label='given EOS')
plt.ylabel('$\epsilon$(p)', fontsize = 15)
plt.xlabel('p', fontsize = 15)
plt.yticks([])
plt.xticks([])
plt.title('a)')
plt.axis([0,2.*np.power(10.,-5.), 0,  eos(2.*np.power(10.,-5.))+np.power(10.,-6.) ] )
plt.legend(loc=2, prop={'size': 14})


plt.subplot(1,3,2)
plt.plot(p1, eos(p1), label='given EOS')
#plt.plot([1.*np.power(10.,-5.), 2.*np.power(10.,-5.)],[eos(1.*np.power(10.,-5.)), eos(2.*np.power(10.,-5.))-0.00005], label='added line')

plt.plot([1.*np.power(10.,-5.), 1.5*np.power(10.,-5.)],[eos(1.*np.power(10.,-5.)), eos(1.5*np.power(10.,-5.))-0.00005], label='added line')

plt.yticks([])
plt.xticks([])
plt.axis([0,2.*np.power(10.,-5.), 0,  eos(2.*np.power(10.,-5.))+np.power(10.,-6.) ] )
plt.title('b)')
plt.legend(loc=2, prop={'size': 14})
plt.xlabel('p', fontsize = 15)
plt.ylabel('$\epsilon$(p)', fontsize = 15)


plt.subplot(1,3,3)
plt.plot(p1, eos(p1), label='given EOS')
plt.plot([1.*np.power(10.,-5.), 1.5*np.power(10.,-5.)],[eos(1.*np.power(10.,-5.)), eos(1.5*np.power(10.,-5.))-0.00005], label='1st line')
plt.plot([1.5*np.power(10.,-5.), 2.*np.power(10.,-5.)],[eos(1.5*np.power(10.,-5.))-0.00005, eos(2.*np.power(10.,-5.))-0.00003], label='2nd line')
plt.yticks([])
plt.xticks([])
plt.axis([0,2.*np.power(10.,-5.), 0,  eos(2.*np.power(10.,-5.))+np.power(10.,-6.) ] )
plt.title('c)')
plt.legend(loc=2, prop={'size': 14})
plt.xlabel('p',fontsize = 15)
plt.ylabel('$\epsilon$(p)',fontsize = 15)

plt.tight_layout()
plt.show()
