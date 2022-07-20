import numpy as np
import matplotlib.pyplot as plt
import pyvalem.formula as fo
from molmass import Formula
import astropy.units as u
from astropy.constants import astropyconst20 as const 
from HD142527.properties import *

# Checking Gauss for 
r = np.linspace(0.1, 375, 3000)
sigma = 3000 * r**(-0.75)
print(np.round(r, 1))

rgap = 0.7
wgap = 0.7
sgap = 10**-2
# sgap = ((wgap/2.0) * np.sqrt(2 * np.pi))**-1

rgap_ind = np.where(np.round(r, 1) == rgap)[0][0]
print(rgap_ind)

print(rgap - wgap/2.0)
wgap_ind1 = np.where(np.round(r, 1) == np.round(rgap - wgap/2.0, 1))[0][0]
wgap_ind2 = np.where(np.round(r, 1) == np.round(rgap + wgap/2.0, 1))[0][0]
print(wgap_ind1, wgap_ind2)

for i in range(wgap_ind1, wgap_ind2 + 1):
	
	if i <= rgap_ind:
		
		m = (sgap * sigma[rgap_ind] - sigma[wgap_ind1])/(r[rgap_ind] - r[wgap_ind1])
		sigma[i] = sigma[wgap_ind1] + m * (r[i] - r[wgap_ind1])
	
	else:
		
		m = (sigma[rgap_ind] - sigma[wgap_ind2])/(r[rgap_ind] - r[wgap_ind2])
		sigma[i] = sigma[rgap_ind] + m * (r[i] - r[rgap_ind])
		 
	
	
	# sigma[i] = 1 - (sgap * np.exp(-(r[i] - rgap)**2 / (2 * (wgap/2.0)**2)))
	# sigma[i] = 1 - (np.exp(-(r[i] - rgap)**2 / (2 * (wgap/2.0)**2)))	
	# sigma[i] = sigma[i] * sgap

plt.plot(r, sigma)
plt.show()

printff

disk = 'HD144432' 										# The name of the disk
Folder = disk + '/'  								    # Path where output files are saved
file = '{0}{1}_Static_Conc.dat'.format(Folder, disk)    # Simulation output file
print(file)

print(p)
print(q)

Na = const.N_A.cgs                    		# Avogadro's number in /mol
solid = 'Mg0.7Fe0.3SiO3'
f = Formula(solid)
molwt = f.mass * u.g / u.mol / Na
print(molwt)

f = fo.Formula(solid)
fancy = "$" + f.latex + "$"
raw_solid = r"{}".format(fancy)
print(raw_solid)
