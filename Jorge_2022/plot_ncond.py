import numpy as np
import matplotlib.pyplot as plt

file = 'Static_Conc.dat'

data   = open(file)
dummy  = data.readline()
dimens = data.readline()
dimens = np.array(dimens.split())
NELEM  = int(dimens[0])
NMOLE  = int(dimens[1])
NDUST  = int(dimens[2])
NPOINT = int(dimens[3])
header = data.readline()
data.close()

keyword = np.array(header.split())

dat = np.loadtxt(file,skiprows=3)

ncond = dat[:, 6]
NPOINT = len(dat[0:])
T = dat[:,0]

index = np.where(keyword=='SMg2SiO4')[0][0]
Mg2SiO4 = dat[:,index]
print(Mg2SiO4)

plt.plot(T, Mg2SiO4)
plt.xlabel("Temperature T (K)")
plt.ylabel("log10 n_solid/ <nH>")
plt.title("God Please: Mg2SiO4")
plt.savefig("SMg2SiO4_vs_T.png")
plt.show()
