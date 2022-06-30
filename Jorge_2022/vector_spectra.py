# I live to make my own life harder every step of the way
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from molmass import Formula
from pyvalem.formula import Formula
from jorge_diskprop import inner_radius, r_from_T
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from opacities import latex_name, molecular_weight, surface_density, Qcurve_plotter, Plancks, tau_calc, flux_map, f
from scipy.interpolate import UnivariateSpline
plt.rcParams['axes.titlesize'] = 10


# Some constants in CGS
Na = 6.022E23                    # Avogadro's number in /mol
h = 6.6261e-27                   # Planck's constant in cm^2 g s-1
c = 2.99792458e10                # Speed of light in cm/s                
k = 1.3807e-16                   # Boltzmann constant in cm^2 g s^-2 K^-1
bar = 1.E+6                      # 1 bar in dyn/cm^2
AU = 1.496E13                    # 1 astronomical unit in cm
Jy = 1.E-23                      # Converting flux from CGS units to Jy



def vector_spectra(tau, I, R_arr, lamda, Rmin, Rmax):
	
	"""
	Plots the integrated flux vs wavelength
	
	Summ shape (lamda, 1)
	"""
	
	# Finding the indices of Rmin and Rmax by taking the first instance of where they are in the rounded R_arr 
	R_rounded = np.round(R_arr, 2)	
	Rmin_id = np.where(R_rounded == Rmin)[0][0]
	Rmax_id = np.where(R_rounded == Rmax)[0][0]
	
	R_arr = R_arr * AU
	print(R_arr)
	R_arr = R_arr[np.newaxis]
	R_arr = R_arr.T 
	I = I.T
	
	delr = (R_arr[Rmin_id+1: Rmax_id+1, :] - R_arr[Rmin_id: Rmax_id, :])
	print("R_arr: ", R_arr.shape) # 500,1
	print("delr: ", delr.shape)   # 499,1
	print(delr)
	print(tau.shape, I.shape)     # 500,450
	extra = np.zeros(len(lamda))
	
	f1 = R_arr[Rmin_id: Rmax_id, :] * tau[Rmin_id: Rmax_id, :] * I[Rmin_id: Rmax_id, :] * 2 * np.pi * delr * 0.5
	# f1 = np.vstack([f1, extra])
	f2 = R_arr[Rmin_id+1: Rmax_id+1, :] * tau[Rmin_id+1: Rmax_id+1, :] * I[Rmin_id+1: Rmax_id+1, :] * 2 * np.pi * delr * 0.5
	# f2 = np.vstack([extra, f2])
	
	print("f1: ", f1.shape)  # 499,450
	print("f2: ", f2.shape)  # 499,450
	
	temp = (f1 + f2)
	summ = temp.sum(axis=0)
	print(summ.shape)        # 450,1
	
	"""
	for r1 in range(Rmin_id, Rmax_id-1):
		for r2 in range(r1+1, Rmax_id):
	
			# Numerical integration using the trapezoidal rule
			fr1 = f(tau[r1, :], I[:, r1], R_arr[r1])
			fr2 = f(tau[r2, :], I[:, r2], R_arr[r2])
			summ += delr * 0.5 * (fr1 + fr2)
	"""
	return summ / Jy

    
def main():
	
	file   = 'Sun/sun_Static_Conc.dat'      # Simulation output file
	data   = open(file)
	dummy  = data.readline()                # Ignoring first line
	dimens = data.readline()                
	dimens = np.array(dimens.split())
	NELEM  = int(dimens[0])                 # Total number of elements used
	NMOLE  = int(dimens[1])                 # Number of molecules created
	NDUST  = int(dimens[2])                 # Number of condensates created
	NPOINT = int(dimens[3])                 # Number of points in simulation
	header = data.readline()                # Saves parameter names such as molecule names
	data.close()
	
	dat = np.loadtxt(file,skiprows=3)
	keyword = np.array(header.split())      # Array of parameter names such as molecule names 
	
	# Extracting data from the output file
	Tg    = dat[:,0]                        # T [K]
	nHtot = dat[:,1]                        # n<H> [cm^-3]
	lognH = np.log10(nHtot)          
	press = dat[:,2]                        # p [dyn/cm^2]
	Tmin  = np.min(Tg)                      # Minimum gas temperature
	Tmax  = np.max(Tg)                      # Maximum gas temperature
	
	# Converting temperatures to corresponding radii
	T0 = 1500                               # Sublimation temperature (K)
	Qr = 1                                  # Ratio of absorption efficiencies (assumed to be black body)
	R_sun = 0.00465047                      # Sun's radius (AU)
	T_sun = 5780                            # Effective temperature of the sun (K)
	
	R_in = inner_radius(Qr, T0, R_sun, T_sun)
	R_arr = r_from_T(R_in, Tg, T0)
	
	top = 5                                 # Top X condensates whose abundance is the highest	
	lmin = 4.0 								# Lower limit of wavelength (microns)
	lmax = 20.0 							# Upper limit of wavelength (microns)
	lsize = 450 							# Number of wavelength (and kappa) points 
	Rmin = 0.03 							# Minimum radius for spectrum plotting (AU) ENSURE IT IS ONLY 2 DECIMAL PLACES LONG
	Rmax = 1.0 							# Maximum radius for spectrum plotting (AU) ENSURE IT IS ONLY 2 DECIMAL PLACES LONG
	gs = 0.1E-4                             # Grain radius (cm)
	
	# All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
	minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])
	one = ['CaMgSi2O6']
	
	# Finding the most abundant condensates
	abundances, solid_names, abunds_dict = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
	top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names)
	top5_solids, topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds, abunds_dict)
	
	# Calculating the surface density
	molwt = molecular_weight(one)
	surf_dens = surface_density(one, molwt, topabunds_radii, nHtot)
	
	# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
	opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278}
	
	# Plotting the Qcurves	
	lamdas = {key: None for key in one}
	kappas = {key: None for key in one}
	rvs = {key: None for key in one}
	fmaxs = {key: None for key in one}
	
	for opfile, density in opfile_dens.items():
		mineral, rv, fmax, lamda, kappa = Qcurve_plotter(opfile, density, gs, lmin, lmax, lsize)
		lamdas[mineral] = lamda
		kappas[mineral] = kappa
		rvs[mineral] = rv
		fmaxs[mineral] = fmax
	
	print(lamdas)        
	# Since some opacity files are missing at the moment, so some lamda and kappa values are None
	
	# Plotting the flux map and calculating the integrated flux for each solid
	I = {key: None for key in one}
	tau = {key: None for key in one}
	F_map = np.zeros((NPOINT, lsize))
	intflux_sum = np.zeros(lsize)
	
	for solid in one:
							
		I[solid] = Plancks(T0, R_arr, R_in, lamdas[solid]) 
		tau[solid] = tau_calc(surf_dens[solid], kappas[solid])
		F_map += flux_map(solid, rvs[solid], fmaxs[solid], tau[solid], I[solid], lamdas[solid], R_arr)
		intflux_sum += vector_spectra(tau[solid], I[solid], R_arr, lamdas[solid], Rmin, Rmax)
	
	
	# Plotting the overall spectrum
	fig = plt.figure()
	# lamda_array = np.linspace(lmin, lmax, lsize)
	plt.plot(lamdas['CaMgSi2O6'], intflux_sum)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux [Jy]')
	plt.title(r'Overall Spectrum r=0.1 $\mu$m R={0}-{1} AU'.format(Rmin, Rmax))
	# plt.savefig("Temp/Overall_spectrum_r0.1_R{0}-{1}.png".format(Rmin, Rmax))
	plt.show()
	

if __name__ == "__main__":
    main()
