import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.constants import astropyconst20 as const
from astropy.modeling.models import BlackBody
from fancy_name import latex_name
from all_solids import get_all_solids
from molmass import Formula
from diskprop import inner_radius, r_from_T
from clean_plots import T_plot
from radial_plot import R_plot
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from spectra import molecular_weight, surface_density, r_to_rad, slice_lQ, get_l_and_k, Qcurve_plotter, Plancks, tau_calc, flux_map, plot_fluxmap, calculate_spectra, plot_spectra, hankel_transform
from scipy.interpolate import UnivariateSpline
from scipy.special import j0


# Some constants in CGS
Na = const.N_A.cgs                    		# Avogadro's number in /mol
h = const.h.cgs                   			# Planck's constant in cm^2 g s-1
c = const.c.cgs              				# Speed of light in cm/s              
k = const.k_B.cgs                   		# Boltzmann constant in cm^2 g s^-2 K^-1
bar = 1.E+6                      			# 1 bar in dyn/cm^2
AU = const.au.cgs                			# 1 astronomical unit in cm
mp = 1.6733E-24         					# Mass of proton (g)
mu = 2.3                					# Mean molecular weight
Rc = 8.314E7            					# Ideal gas constant (erg K^-1 mol^-1)
G = 6.6725E-8           					# Gravitational constant (cm3 g^-1 s^-2)

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
	Tg    = dat[:,0] * u.K                        # T [K]
	nHtot = dat[:,1] / u.cm**3                    # n<H> [cm^-3]          
	press = dat[:,2]                        	  # p [dyn/cm^2]
	Tmin  = np.min(Tg)                            # Minimum gas temperature
	Tmax  = np.max(Tg)                      	  # Maximum gas temperature
	
	# Converting temperatures to corresponding radii
	T0 = 1500.0 * u.K                          	# Dust sublimation temperature (K)
	Qr = 1                                  	# Ratio of absorption efficiencies (assumed to be black body)
	R_star = const.R_sun.to(u.AU)               # Star's radius (AU)
	T_star = 5780 * u.K                         # Effective temperature of the star (K)
	Sigma0 = 1700 * u.g / u.cm**2          		# Surface density with MMSN (g/cm^2)
	M_star = 1.99E33         					# Solar mass (g)
	
	R_in = inner_radius(Qr, T0, R_star, T_star)   # Inner-most radius beyond which the dust is sublimated (AU)
	R_arr = r_from_T(R_in, Tg, T0)                # 1D array of radii obtained from the power law disk model (AU)
	
	top = 5                                 	  			# Top X condensates whose abundance is the highest	
	lmin = 0.0 * u.micron 						  			# Lower limit of wavelength (microns)
	lmax = 20.0 * u.micron						  			# Upper limit of wavelength (microns)
	lsize = 450 								  			# Number of wavelength (and kappa) points 
	Rmin = np.round(np.min(R_arr), 3) 						# Minimum radius for spectrum plotting (AU) ENSURE IT IS ONLY 2 DECIMAL PLACES LONG
	Rmax = np.round(np.max(R_arr), 3)						# Maximum radius for spectrum plotting (AU) ENSURE IT IS ONLY 2 DECIMAL PLACES LONG
	dist_pc = 100 * u.pc  			            			# Assuming a distance to the Sun-like star in parsec
	gs = 0.1E-4 * u.cm                            			# Grain radius (cm)
	wl = 5.5 * u.micron                           			# Observing wavelength (microns)
	wl_list = [1.0, 2.0, 3.2, 5.5, 10.0, 12.0] * u.micron	# 1D list of wavelengths to plot correlated flux against baselines (microns)
	B = np.arange(0.0, 130.0, 2.0) * u.m          			# 1D array of baselines (m)
	B_small = np.linspace(0.0, 130.0, 5) * u.m    			# 1D array of a few baselines to plot correlated flux against wavelengths (m)
	folder = 'Temp/'                                        # Folder where all the results go
	
	minerals = get_all_solids(keyword, dat, NELEM, NMOLE, NDUST)
	
	# Plotting the abunances as a function of radius and temperature
	# R_plot(minerals, dat, keyword, R_arr, R_in, Rmin, Rmax, T0, folder, NELEM, NMOLE, NDUST)
	
	# Finding the most abundant condensates
	abundances, solid_names, abunds_dict = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
	top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names)
	top5_solids, topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds, abunds_dict)
	
	# Write down abundances and corresponding solids element by element in a file in a human readable format 
	filename = folder + 'most_abundant.dat'
	with open(filename, 'w') as f:
		
		f.write('{0:47} {1:47} {2:47} {3:47} {4:47} Radius \n'.format(str(1), str(2), str(3), str(4), str(5)))                  # Adding headers with some formatting to help readability
		for radius in range(len(top_abunds)):
			for element in range(len(top_abunds[radius])):
				f.write('{0:20} : {1:20} '.format(str(top_solids[radius][element]), str(top_abunds[radius][element])))          # Adding the top 5 most abundant solids and their corresponding abundances
				if (element+1) % top == 0:
					f.write(str(R_arr[radius]) + '\n')                                                                          # If all 5 solids are written, then add the radius in the last column and move to the next line
				else:
					f.write('\t ')
	
	# Removing the solids without opfiles from top5_solids
	not_there = ['SiO', 'Mg3Si4O12H2', 'Fe3Si2O9H4', 'Ni', 'NaAlSi3O8', 'NaMg3AlSi3O12H2', 'CaAl2Si2O8', 'H2O']
	top5_solids = np.setdiff1d(top5_solids, not_there)
	
	# Calculating the surface density
	molwt = molecular_weight(top5_solids)
	surf_dens = surface_density(top5_solids, molwt, topabunds_radii, nHtot)
	
	# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
	opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64}
	# , 'Qcurve_inputs/Q_H2O_rv0.1_fmax0.7.dat' : 0.92}
	
	# Adding units to the material densities using astropy
	for key, value in opfile_dens.items():
		opfile_dens[key] = value * u.g / u.cm**3
	
	# Initializing dictionaries for the wavelengths, opacities, grain sizes and emptiness fractions used in the opfiles for each solid 	
	lamdas = {key: None for key in top5_solids}
	kappas = {key: None for key in top5_solids}
	rvs = {key: None for key in top5_solids}
	fmaxs = {key: None for key in top5_solids}
	
	for opfile, density in opfile_dens.items():
		mineral, rv, fmax, lamda, kappa = get_l_and_k(opfile, density, gs, lmin, lmax, lsize)
		# Qcurve_plotter(lamda, kappa, mineral, rv, fmax, folder)
		lamdas[mineral] = lamda
		kappas[mineral] = kappa
		rvs[mineral] = rv
		fmaxs[mineral] = fmax		
	
	# Plotting the flux map and calculating the integrated flux for each solid
	I = {key: None for key in top5_solids}
	tau = {key: None for key in top5_solids}
	F_map_sum = np.zeros((NPOINT, lsize)) * u.erg / (u.s * u.Hz * u.sr * u.cm**2)
	intflux_sum = np.zeros(lsize) * u.Jy
	
	for solid in top5_solids:
			
			I[solid] = Plancks(T0, R_arr, R_in, lamdas[solid]) 
			tau[solid] = tau_calc(surf_dens[solid], kappas[solid])
			
			F_map = flux_map(tau[solid], I[solid])
			# plot_fluxmap(solid, rvs[solid], fmaxs[solid], F_map, lamdas[solid], R_arr, folder)
			F_map_sum += F_map
			
			intflux = calculate_spectra(tau[solid], I[solid], R_arr, Rmin, Rmax)
			# plot_spectra(lamdas[solid], intflux, solid, rvs[solid], fmaxs[solid], Rmin, Rmax, folder)
			intflux_sum += intflux
			
	
	# Plotting the overall flux map
	fig, ax = plt.subplots(1,1)
	img = ax.imshow(F_map_sum, cmap='plasma', interpolation='none')
	
	x_axis_locations = np.linspace(0, lsize-1, 8).astype(int)
	ax.set_xticks(x_axis_locations)
	x_axis_labels = np.round(lamda[x_axis_locations], 1)
	ax.set_xticklabels(x_axis_labels.value)
	ax.set_xlabel(r'$\lambda$ ($\mu$m)')
	
	y_axis_locations = np.linspace(0, len(R_arr)-1, 8).astype(int)
	ax.set_yticks(y_axis_locations)
	y_axis_labels = np.round(R_arr[y_axis_locations], 3)
	ax.set_yticklabels(y_axis_labels.value)
	ax.set_ylabel('R (AU)')
	
	ax.set_title(r"Overall Flux Map for r=0.1 microns")
	fig.colorbar(img)
	plt.savefig(folder + "overall_fluxmap.png", bbox_inches = 'tight')
	plt.show()
	
	# Plotting the overall spectrum
	fig = plt.figure()
	plt.plot(lamdas['Mg2SiO4'], intflux_sum)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title(r'Overall Spectrum r=0.1 $\mu$m R={0}-{1} AU'.format(Rmin.value, Rmax.value))
	plt.savefig(folder + "Overall_spectrum_r0.1_R{0}-{1}.png".format(Rmin.value, Rmax.value))
	plt.show()
	
	# Plotting the overall spectrum considering multiple radii together
	Rmin_list = [0.035, 0.035, 0.1, 0.1, 0.5, 0.75] * u.AU
	Rmax_list = [0.05, 0.1, 0.25, 0.5, 1.0, 1.277] * u.AU
	intflux_sum_mr = np.zeros(lsize) * u.Jy
	fig = plt.figure()
	
	for i in range(len(Rmax_list)):		
		for solid in top5_solids:		 			
			intflux_sum_mr += calculate_spectra(tau[solid], I[solid], R_arr, Rmin_list[i], Rmax_list[i])
			
		plt.plot(lamdas['Mg2SiO4'], intflux_sum_mr, label=r"($R_{{min}}$, $R_{{max}}$) = ({0},{1}) AU".format(Rmin_list[i].value, Rmax_list[i].value))
		intflux_sum_mr = np.zeros(lsize) * u.Jy
			
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title(r'Overall spectrum for multiple radii')
	plt.legend()	
	plt.savefig(folder + "Overall_spectrum_multiple_radii_limits.png")
	plt.show()
	
	# Plotting the correlated flux density for multiple baselines against wavelengths
	fig = plt.figure()
	
	for Bl in B_small:
		
		corr_flux_absB = hankel_transform(F_map_sum, R_arr, lamdas['Mg2SiO4'], wl, Bl, wl_array = True)
		plt.plot(lamdas['Mg2SiO4'], corr_flux_absB, label = 'B = {0} m'.format(Bl.value))
	
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Correlated flux (Jy)')             
	plt.title(r'Correlated flux for multiple baselines, rv = 0.1 $\mu$m')
	plt.legend()
	plt.savefig(folder + 'Correlated_flux_multB_rv0.1.png')
	plt.show()
	
	# Plotting the correlated flux density for a single wavelength against baselines
	inter_flux = np.zeros(len(B)) * u.Jy
	
	for wl in wl_list:
		for bl in range(len(B)):			
			inter_flux[bl] = hankel_transform(F_map_sum, R_arr, lamdas['Mg2SiO4'], wl, B[bl], wl_array = False)
		plt.plot(B, inter_flux, label=r"{0} $\mu$m".format(wl.value))	
		inter_flux = np.zeros(len(B)) * u.Jy		
	
	plt.xlabel(r'Baseline B (m)')
	plt.ylabel('Correlated flux (Jy)')
	plt.title(r'Correlated flux for multiple wavelengths, rv = 0.1 $\mu$m')
	plt.legend()
	plt.savefig(folder + "Correlated_flux_mulWL_rv0.1.png")
	plt.show()


if __name__ == "__main__":
    main()
