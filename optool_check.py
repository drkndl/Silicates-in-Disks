# Will I finish optool or will optool finish me?

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import pandas as pd

filename = 'optool_for_a0.1_f0.999.dat'
a = 0.1 * u.micron
rho = 3.27 * u.g / u.cm**3
Q_curve = pd.read_csv(filename, delimiter='\s+', skiprows = 23, names = ['wavelength', 'Kabs', 'Ksca', 'g_asymm'])
lamda = u.Quantity(Q_curve['wavelength'].to_numpy(), u.micron)
kappa = Q_curve['Kabs'].to_numpy() * u.cm**2 / u.g
Q = 4.0 * a.to(u.cm) * rho * kappa / 3.0
print(Q.unit)

plt.plot(lamda, Q)
plt.xlabel(r'$\lambda$ ($\mu$m)')
plt.ylabel(r'Absorption efficiency (Q)')
plt.title("Optool Qcurve for Forsterite a=0.1 f=0.999 wl = 7-14 $\mu$m")
plt.savefig("optool_for_a0.1_f0.999.png")
plt.show()
