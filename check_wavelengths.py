import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from spectra import get_l_and_k 
import pandas as pd


def get_req_lamda(opfile, lmin, lmax):
	
	Q_curve = pd.read_csv(opfile, delimiter='\s+', skiprows=1, names=['wavelength','Q'])	
	lamda = Q_curve['wavelength'].to_numpy() * u.micron
	indices = np.where(np.logical_and(lamda >= lmin, lamda <= lmax))
	lamda = lamda[indices]
	
	return lamda
	
	
top5_solids = ['CaMgSi2O6', 'Fe', 'Fe2SiO4', 'Fe3O4', 'FeS', 'Mg2SiO4', 'Mg3Si2O9H4', 'MgAl2O4', 'MgSiO3', 'Olivine', 'Pyroxene']
lmin = 5.0 * u.micron 						  			# Lower limit of wavelength (microns) 
lmax = 21.0 * u.micron						  			# Upper limit of wavelength (microns)
lsize = 450 								  			# Number of wavelength (and kappa) points 

# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv2.0.dat' : 3.2, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv2.0.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe3O4_rv2.0_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe2SiO4_rv2.0_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_Fe_met_rv2.0_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_FeS_rv2.0_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv2.0_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64, 'Qcurve_inputs/qval_MgAl2O4_rv2.0_fmax0.7.dat' : 3.64, 'Qcurve_inputs/Q_Olivine_rv0.1_fmax0.7.dat': 3.71, 'Qcurve_inputs/Q_Olivine_rv2.0_fmax0.7.dat' : 3.71, 'Qcurve_inputs/Q_Pyroxene_rv0.1_fmax0.7.dat': 3.01, 'Qcurve_inputs/Q_Pyroxene_fmax0.7_rv2.0.dat': 3.01}

# Adding units to the material densities using astropy
for key, value in opfile_dens.items():
	opfile_dens[key] = value * u.g / u.cm**3

req_lamda = np.round(get_req_lamda('Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat', lmin, lmax), 2)
print(req_lamda)

# Initializing dictionaries for the wavelengths, opacities, grain sizes and emptiness fractions used in the opfiles for each solid 	
gs_ranges = ['0.1', '2.0']
lamdas = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
kappas = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
fmaxs = {key: None for key in top5_solids}

for opfile, density in opfile_dens.items():
	mineral, rv, fmax, lamda, kappa = get_l_and_k(opfile, density, lmin, lmax, lsize, req_lamda)		
	# Qcurve_plotter(lamda, kappa, mineral, rv, fmax, folder)
	lamdas[mineral][rv] = lamda
	kappas[mineral][rv] = kappa
	fmaxs[mineral] = fmax


for solid in top5_solids:
	for size in gs_ranges:
		
		if solid=="CaMgSi2O6" and size=='2.0':
			print("Skipping: ", solid, size)
			continue
				
		x = range(len(lamdas[solid][size]))
		print(len(x))
		plt.plot(x, lamdas[solid][size], label="{0} {1}".format(solid, size))
		
plt.title("Checking wavelength distributions for all solids")
plt.legend()
plt.savefig("check_wavelengths.png")
plt.show()
