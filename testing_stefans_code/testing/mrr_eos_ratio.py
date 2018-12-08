import matplotlib.pyplot as plt
import numpy as np
import csv

P = 7.55616208*np.power(10., 5.)

def eos(p):
    return np.power(p/(10.), 3./5.)

Rl2 = []
Ml2 = []
el2 = []
pl2 = []

Rl245 = []
Ml245 = []
el245 = []
pl245 = []

Rl25 = []
Ml25 = []
el25 = []
pl25 = []

Rl255 = []
Ml255 = []
el255 = []
pl255 = []

Rl257 = []
Ml257 = []
el257 = []
pl257 = []

Rl26 = []
Ml26 = []
el26 = []
pl26 = []

Rl263 = []
Ml263 = []
el263 = []
pl263 = []

Rl265 = []
Ml265 = []
el265 = []
pl265 = []

Rl27 = []
Ml27 = []
el27 = []
pl27 = []

Rl275 = []
Ml275 = []
el275 = []
pl275 = []

Rl28 = []
Ml28 = []
el28 = []
pl28 = []

Rd = []
Md = []
ed = []
pd = []

with open('results_l2', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl2.append(float(row[1]))
        Ml2.append(float(row[3]))
        pl2.append(float(row[4]))
        el2.append(float(row[5]))

with open('results_l2_45', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl245.append(float(row[1]))
        Ml245.append(float(row[3]))
        pl245.append(float(row[4]))
        el245.append(float(row[5]))

with open('results_l2_5', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl25.append(float(row[1]))
        Ml25.append(float(row[3]))
        pl25.append(float(row[4]))
        el25.append(float(row[5]))

with open('results_l2_55', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl255.append(float(row[1]))
        Ml255.append(float(row[3]))
        pl255.append(float(row[4]))
        el255.append(float(row[5]))

with open('results_l2_57', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl257.append(float(row[1]))
        Ml257.append(float(row[3]))
        pl257.append(float(row[4]))
        el257.append(float(row[5]))

with open('results_l2_6', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl26.append(float(row[1]))
        Ml26.append(float(row[3]))
        pl26.append(float(row[4]))
        el26.append(float(row[5]))

with open('results_l2_63', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl263.append(float(row[1]))
        Ml263.append(float(row[3]))
        pl263.append(float(row[4]))
        el263.append(float(row[5]))

with open('results_l2_65', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl265.append(float(row[1]))
        Ml265.append(float(row[3]))
        pl265.append(float(row[4]))
        el265.append(float(row[5]))

with open('results_l2_7', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl27.append(float(row[1]))
        Ml27.append(float(row[3]))
        pl27.append(float(row[4]))
        el27.append(float(row[5]))

with open('results_l2_75', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl275.append(float(row[1]))
        Ml275.append(float(row[3]))
        pl275.append(float(row[4]))
        el275.append(float(row[5]))

with open('results_l2_8', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rl28.append(float(row[1]))
        Ml28.append(float(row[3]))
        pl28.append(float(row[4]))
        el28.append(float(row[5]))

with open('mr.out', 'rb') as csvfile:
    plots = csv.reader(csvfile, delimiter = ' ')
    for row in plots:
        Rd.append(float(row[0]))
        Md.append(float(row[1]))
        pd.append(float(row[2]))
        ed.append(float(row[3]))

plt.figure(figsize=(14,7))

el2n = [a*P for a in el2]
pl2n = [a*P for a in pl2]
edn = [P*eos(a) for a in pl2]
eRatL2 = []
endind = len(edn)
for i in range(endind):
    eRatL2.append(el2n[i]/edn[i])

el245n = [a*P for a in el245]
pl245n = [a*P for a in pl245]
edn = [P*eos(a) for a in pl245]
eRatL245 = []
endind = len(edn)
for i in range(endind):
    eRatL245.append(el245n[i]/edn[i])

el25n = [a*P for a in el25]
pl25n = [a*P for a in pl25]
edn = [P*eos(a) for a in pl25]
eRatL25 = []
endind = len(edn)
for i in range(endind):
    eRatL25.append(el25n[i]/edn[i])

el255n = [a*P for a in el255]
pl255n = [a*P for a in pl255]
edn = [P*eos(a) for a in pl255]
eRatL255 = []
endind = len(edn)
for i in range(endind):
    eRatL255.append(el255n[i]/edn[i])

el257n = [a*P for a in el257]
pl257n = [a*P for a in pl257]
edn = [P*eos(a) for a in pl257]
eRatL257 = []
endind = len(edn)
for i in range(endind):
    eRatL257.append(el257n[i]/edn[i])

el26n = [a*P for a in el26]
pl26n = [a*P for a in pl26]
edn = [P*eos(a) for a in pl26]
eRatL26 = []
endind = len(edn)
for i in range(endind):
    eRatL26.append(el26n[i]/edn[i])

el263n = [a*P for a in el263]
pl263n = [a*P for a in pl263]
edn = [P*eos(a) for a in pl263]
eRatL263 = []
endind = len(edn)
for i in range(endind):
    eRatL263.append(el263n[i]/edn[i])

el265n = [a*P for a in el265]
pl265n = [a*P for a in pl265]
edn = [P*eos(a) for a in pl265]
eRatL265 = []
endind = len(edn)
for i in range(endind):
    eRatL265.append(el265n[i]/edn[i])

el27n = [a*P for a in el27]
pl27n = [a*P for a in pl27]
edn = [P*eos(a) for a in pl27]
eRatL27 = []
endind = len(edn)
for i in range(endind):
    eRatL27.append(el27n[i]/edn[i])

el275n = [a*P for a in el275]
pl275n = [a*P for a in pl275]
edn = [P*eos(a) for a in pl275]
eRatL275 = []
endind = len(edn)
for i in range(endind):
    eRatL275.append(el275n[i]/edn[i])

el28n = [a*P for a in el28]
pl28n = [a*P for a in pl28]
edn = [P*eos(a) for a in pl28]
eRatL28 = []
endind = len(edn)
for i in range(endind):
    eRatL28.append(el28n[i]/edn[i])

p1 = np.arange(0., 0.00025, 0.00000001)
oneX = [25, 160]
oneY = [1, 1]

#plt.subplot(2,2,(1,3))
#plt.xlim(11.5, 14)
#plt.ylim(1.1, 1.2)
#plt.plot(Rd, Md, label = 'MRR input')
#plt.plot(Rl2, Ml2,label = 'MRR $\lambda = 0.01$')
#plt.plot(Rl245, Ml245,label = 'MRR $\lambda = 0.045$')
#plt.plot(Rl25, Ml25,label = 'MRR $\lambda = 0.05$')
#plt.plot(Rl255, Ml255,label = 'MRR $\lambda = 0.055$')
#plt.plot(Rl257, Ml257,label = 'MRR $\lambda = 0.057$')
#plt.plot(Rl26, Ml26,label = 'MRR $\lambda = 0.06$')
#plt.plot(Rl263, Ml263,label = 'MRR $\lambda = 0.063$')
#plt.plot(Rl265, Ml265,label = 'MRR $\lambda = 0.065$')
#plt.plot(Rl27, Ml27,label = 'MRR $\lambda = 0.07$')
#plt.plot(Rl275, Ml275,label = 'MRR $\lambda = 0.075$')
#plt.plot(Rl28, Ml28,label = 'MRR $\lambda = 0.08$')
#plt.ylabel('M/M$_\odot$', fontsize=15)
#plt.xlabel('R/km', fontsize=15)
#plt.legend(loc=1,prop={'size':12})

#plt.subplot(2,2,(2,4))
#plt.plot(p1*P, eos(p1)*P, label = 'EOS input')
#plt.plot(pl2n,   eRatL2,   label ='EOS $\lambda = 0.01$')
#plt.plot(pl245n, eRatL245, label ='EOS $\lambda = 0.045$')
plt.plot(pl25n,  eRatL25,  label ='EOS $\lambda = 0.05$')
#plt.plot(pl255n, eRatL255, label ='EOS $\lambda = 0.055$')
#plt.plot(pl257n, eRatL257, label ='EOS $\lambda = 0.057$')
#plt.plot(pl26n,  eRatL26,  label ='EOS $\lambda = 0.06$')
#plt.plot(pl263n, eRatL263, label ='EOS $\lambda = 0.063$')
plt.plot(pl265n, eRatL265, label ='EOS $\lambda = 0.065$')
plt.plot(pl27n,  eRatL27,  label ='EOS $\lambda = 0.07$')
plt.plot(pl275n, eRatL275, label ='EOS $\lambda = 0.075$')
#plt.plot(pl28n,  eRatL28,  label ='EOS $\lambda = 0.08$')
plt.plot(oneX, oneY, 'b--')

plt.ylabel('$\epsilon_R$/$\epsilon_D$', fontsize=20)

plt.xlabel('p/MeVfm$^{-3}$', fontsize=20)
plt.legend(loc=2,prop={'size':20})

plt.tight_layout()
plt.show()

