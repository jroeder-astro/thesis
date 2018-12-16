import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)


def eos(p,i):
    return m[i-1] + (p-m[i-2]) * (m[i+1]-m[i-1]) / (m[i]-m[i-2]) 

Rd = []
Md = []
ed = []
pd = []

with open('mr_eos_com.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    for row in plots:
        Rd.append(float(row[0]))
        Md.append(float(row[1]))
        ed.append(float(row[2]))
        pd.append(float(row[3]))

edn = [a*P for a in ed]
pdn = [a*P for a in pd]

m = []
for i in range(len(ed)):
  #  print pd[i]
  #  print ed[i]
    m.append(pd[i])
    m.append(ed[i])
lena = len(m)

R11 = []
M11 = []
e11 = []
p11 = []

with open('results_el1_1', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R11.append(float(row[1]))
        M11.append(float(row[3]))
        e11.append(float(row[5]))
        p11.append(float(row[4]))

e11n = [a*P for a in e11]
p11n = [a*P for a in p11]
eRatL11 = []
eDn = []
for x in p11:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p11)):
    eRatL11.append(e11n[i]/eDn[i])

R21 = []
M21 = []
e21 = []
p21 = []

with open('results_el2_1', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R21.append(float(row[1]))
        M21.append(float(row[3]))
        e21.append(float(row[5]))
        p21.append(float(row[4]))

e21n = [a*P for a in e21]
p21n = [a*P for a in p21]
eRatL21 = []
eDn = []
for x in p21:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p21)):
    eRatL21.append(e21n[i]/eDn[i])

R24 = []
M24 = []
e24 = []
p24 = []

with open('results_el2_4', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R24.append(float(row[1]))
        M24.append(float(row[3]))
        e24.append(float(row[5]))
        p24.append(float(row[4]))

e24n = [a*P for a in e24]
p24n = [a*P for a in p24]
eRatL24 = []
eDn = []
for x in p24:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p24)):
    eRatL24.append(e24n[i]/eDn[i])

R25 = []
M25 = []
e25 = []
p25 = []

with open('results_el2_5', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R25.append(float(row[1]))
        M25.append(float(row[3]))
        e25.append(float(row[5]))
        p25.append(float(row[4]))

e25n = [a*P for a in e25]
p25n = [a*P for a in p25]
eRatL25 = []
eDn = []
for x in p25:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p25)):
    eRatL25.append(e25n[i]/eDn[i])

R255 = []
M255 = []
e255 = []
p255 = []

with open('results_el2_55', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R255.append(float(row[1]))
        M255.append(float(row[3]))
        e255.append(float(row[5]))
        p255.append(float(row[4]))

e255n = [a*P for a in e255]
p255n = [a*P for a in p255]
eRatL255 = []
eDn = []
for x in p255:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p255)):
    eRatL255.append(e255n[i]/eDn[i])

R257 = []
M257 = []
e257 = []
p257 = []

with open('results_el2_57', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R257.append(float(row[1]))
        M257.append(float(row[3]))
        e257.append(float(row[5]))
        p257.append(float(row[4]))

e257n = [a*P for a in e257]
p257n = [a*P for a in p257]
eRatL257 = []
eDn = []
for x in p257:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p257)):
    eRatL257.append(e25n[i]/eDn[i])

R26 = []
M26 = []
e26 = []
p26 = []

with open('results_el2_6', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R26.append(float(row[1]))
        M26.append(float(row[3]))
        e26.append(float(row[5]))
        p26.append(float(row[4]))

e26n = [a*P for a in e26]
p26n = [a*P for a in p26]
eRatL26 = []
eDn = []
for x in p26:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p26)):
    eRatL26.append(e26n[i]/eDn[i])

R27 = []
M27 = []
e27 = []
p27 = []

with open('results_el2_7', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R27.append(float(row[1]))
        M27.append(float(row[3]))
        e27.append(float(row[5]))
        p27.append(float(row[4]))

e27n = [a*P for a in e27]
p27n = [a*P for a in p27]
eRatL27 = []
eDn = []
for x in p27:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p27)):
    eRatL27.append(e27n[i]/eDn[i])

R28 = []
M28 = []
e28 = []
p28 = []

with open('results_el2_8', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R28.append(float(row[1]))
        M28.append(float(row[3]))
        e28.append(float(row[5]))
        p28.append(float(row[4]))

e28n = [a*P for a in e28]
p28n = [a*P for a in p28]
eRatL28 = []
eDn = []
for x in p28:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p28)):
    eRatL28.append(e28n[i]/eDn[i])

R31 = []
M31 = []
e31 = []
p31 = []

with open('results_el3_1', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        R31.append(float(row[1]))
        M31.append(float(row[3]))
        e31.append(float(row[5]))
        p31.append(float(row[4]))

e31n = [a*P for a in e31]
p31n = [a*P for a in p31]
eRatL31 = []
eDn = []
for x in p31:
    for i in range(0, lena, 2):
        if (x < m[i]):
            break
    eDn.append(P*eos(x,i))
for i in range(len(p31)):
    eRatL31.append(e31n[i]/eDn[i])

oneX = [0,700]
oneY = [1,1]
plt.figure(figsize=(12,7))


#plt.plot(pdn, edn, label ='EOS input')
#plt.plot(p11n, e11n, label ='EOS $\lambda = 0.1$')
#plt.plot(p21n, e21n, label ='EOS $\lambda = 0.01$')
#plt.plot(p24n, e24n, label ='EOS $\lambda = 0.04$')
#plt.plot(p25n, e25n, label ='EOS $\lambda = 0.05$')
#plt.plot(p255n, e255n, label ='EOS $\lambda = 0.055$')
#plt.plot(p257n, e257n, label ='EOS $\lambda = 0.057$')
#plt.plot(p26n, e26n, label ='EOS $\lambda = 0.06$')
#plt.plot(p27n, e27n, label ='EOS $\lambda = 0.07$')
#plt.plot(p28n, e28n, label ='EOS $\lambda = 0.08$')
#plt.plot(p31n, e31n, label ='EOS $\lambda = 0.001$')

plt.xlim(70,550)
plt.ylim(0.80, 1.1)
plt.plot(p11n,  eRatL11,  label = 'EOS $\lambda = 0.1$')
plt.plot(p21n,  eRatL21,  label = 'EOS $\lambda = 0.01$')
plt.plot(p24n,  eRatL24,  label = 'EOS $\lambda = 0.04$')
plt.plot(p25n,  eRatL25,  label = 'EOS $\lambda = 0.05$')
plt.plot(p255n, eRatL255, label = 'EOS $\lambda = 0.055$')
plt.plot(p257n, eRatL257, label = 'EOS $\lambda = 0.057$')
plt.plot(p26n,  eRatL26,  label = 'EOS $\lambda = 0.06$')
plt.plot(p27n,  eRatL27,  label = 'EOS $\lambda = 0.07$')
plt.plot(p28n,  eRatL28,  label = 'EOS $\lambda = 0.08$')
#plt.plot(p31n,  eRatL31,  label = 'EOS $\lambda = 0.001$')


plt.plot(oneX,oneY, 'b--')
plt.ylabel('$\epsilon_R$/$\epsilon_D$', fontsize=20)
#plt.ylabel('$\epsilon$(p)/MeVfm$^{-3}$', fontsize=15)
plt.xlabel('p/MeVfm$^{-3}$', fontsize=20)
plt.legend(loc=1,prop={'size':20})

#plt.subplots_adjust(hspace = 0.5)
plt.tight_layout()
plt.show()
#plt.savefig('mrr_eos_ext_upright.pdf')
