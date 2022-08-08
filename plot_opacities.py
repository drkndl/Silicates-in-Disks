import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from astropy import units as u
from astropy.constants import astropyconst20 as const
import cmasher as cmr
from fancy_name import latex_name
from spectra import get_l_and_k

plt.rcParams["font.family"] = "serif"
plt.rcParams['font.serif'] = ['TeX Gyre Schola']
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['figure.titlesize'] = 20

opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_Pyroxene_rv0.1_fmax0.7.dat': 3.01, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_Olivine_rv0.1_fmax0.7.dat': 3.71, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64}

lmin = 8.0 * u.micron 						  			# Lower limit of wavelength (microns)
lmax = 13.0 * u.micron						  			# Upper limit of wavelength (microns)
lsize = 450 								  			# Number of wavelength (and kappa) points 

fig, ax = plt.subplots(1,1, figsize=(6, 13))
n = len(opfile_dens.keys())
scale = np.linspace(0, 0.9, n)
colours = cmr.take_cmap_colors('cmr.guppy', n, cmap_range=(0.1, 0.9))
i = 0

for opfile, density in opfile_dens.items():
	mineral, rv, fmax, lamda, kappa = get_l_and_k(opfile, density, lmin, lmax, lsize)
	kappa = preprocessing.normalize([kappa])
	ax.plot(lamda, kappa[0,:] + scale[i], color = colours[i], linewidth = 2)
	ax.text(0.75, scale[i]+0.07, latex_name(mineral), transform=ax.transAxes, color=colours[i], fontsize = 18)
	i += 1

ax.set_xlabel(r"Wavelength $\lambda$ $(\mu$m)")
ax.set_ylabel(r"Normalized Opacities $\kappa$")
ax.set_title("Opacity curves")
ax.get_yaxis().set_ticks([])
# ~ ax.legend()
fig.tight_layout()
fig.savefig("opacities.png")
plt.show()
