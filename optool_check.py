# Will I finish optool or will optool finish me?

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import pandas as pd
import glob

files = 'Optool_Check/*.dat'
a = 0.1 * u.micron
rho = 3.27 * u.g / u.cm**3

for opfile in glob.glob(files):
	
	filename = opfile.split('/')[1]
	f = filename.split('_')[3][:-4]
	Q_curve = pd.read_csv(opfile, delimiter='\s+', skiprows = 23, names = ['wavelength', 'Kabs', 'Ksca', 'g_asymm'])
	lamda = u.Quantity(Q_curve['wavelength'].to_numpy(), u.micron)
	kappa = Q_curve['Kabs'].to_numpy() * u.cm**2 / u.g
	Q = 4.0 * a.to(u.cm) * rho * kappa / 3.0
	plt.plot(lamda, Q, label=f)

plt.xlabel(r'$\lambda$ ($\mu$m)')
plt.ylabel(r'Absorption efficiency (Q)')
plt.title("Optool Qcurves Forsterite a=0.1 wl = 7-14 $\mu$m")
plt.legend()
plt.savefig("Optool_Check/optool_for_a0.1_multf.png")
plt.show()
