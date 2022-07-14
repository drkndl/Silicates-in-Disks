# Attempts to recreate disk spectra from van Boekel (2005)

import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.constants import astropyconst20 as const
from fancy_name import latex_name
from all_solids import get_all_solids
from molmass import Formula
from diskprop import inner_radius, r_from_T, scale_height, star_radius
from T_plot import T_plot
from radial_plot import R_plot
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from spectra import molecular_weight, surface_density, r_to_rad, slice_lQ, get_l_and_k, Plancks, tau_calc, tau_calc_amorphous, flux_map, calculate_spectra, hankel_transform
from no_thoughts_just_plots import Qcurve_plotter, plot_Bv, plot_tau, plot_fluxmap, plot_spectra
from Qcurve_compare import gimme_k


# ~ def get_paper_spectra(filename):
	
	# ~ return wl, flux
	
	
def gimme_l(filename):
	
	Q_curve = pd.read_csv(filename, delimiter='\s+', skiprows = 23, names = ['wavelength', 'Kabs', 'Ksca', 'g_asymm'])
	lamda = u.Quantity(Q_curve['wavelength'].to_numpy(), u.micron)
	
	return lamda
	
	
def main():
	
	file = 'HD144432/HD144432_Static_Conc.dat'      # Simulation output file
	disk = 'HD144432'
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
	T0 = 1500.0 * u.K                          		# Dust sublimation temperature (K)
	Qr = 3                                  		# Ratio of absorption efficiencies 
	L_star = 10**1.01 * const.L_sun.cgs         	# Stellar luminosity
	T_star = 7345.138 * u.K                         # Effective temperature of the star (K)
	R_star = star_radius(L_star, T_star).to(u.AU)   # Star's radius (AU)
	Sigma0 = 1700 * u.g / u.cm**2          			# Surface density with MMSN (g/cm^2)
	M_star = 1.8 * 1.99E33 * u.g         			# Solar mass (g)
	q = -0.5 										# Disk temperature gradient exponent
	e = -1.0 										# Disk surface density gradient exponent
	dist_pc = 145 * u.pc                            # Star distance in parsec
	
	R_in = inner_radius(Qr, T0, R_star, T_star)   # Inner-most radius beyond which the dust is sublimated (AU)
	R_arr = r_from_T(R_in, Tg, T0, q)             # 1D array of radii obtained from the power law disk model (AU)
	# H = scale_height(M_star, R_arr, Tg)
	H = 1.0 * u.cm 								  # Scale height (cm)
	
	top = 5                                 	  			# Top X condensates whose abundance is the highest	
	lmin = 0.0 * u.micron 						  			# Lower limit of wavelength (microns)
	lmax = 20.0 * u.micron						  			# Upper limit of wavelength (microns)
	lsize = 300 								  			# Number of wavelength (and kappa) points 
	Rmin = np.round(np.min(R_arr), 3) 						# Minimum radius for spectrum plotting (AU) ENSURE IT IS ONLY 2 DECIMAL PLACES LONG
	Rmax = np.round(np.max(R_arr), 3)						# Maximum radius for spectrum plotting (AU) ENSURE IT IS ONLY 2 DECIMAL PLACES LONG
	dist_pc = 100 * u.pc  			            			# Assuming a distance to the Sun-like star in parsec
	gs = 0.1E-4 * u.cm                            			# Grain radius (cm)
	wl = 5.5 * u.micron                           			# Observing wavelength (microns)
	wl_list = [1.0, 2.0, 3.2, 5.5, 10.0, 12.0] * u.micron	# 1D list of wavelengths to plot correlated flux against baselines (microns)
	B = np.arange(0.0, 130.0, 2.0) * u.m          			# 1D array of baselines (m)
	B_small = np.linspace(0.0, 130.0, 5) * u.m    			# 1D array of a few baselines to plot correlated flux against wavelengths (m)
	folder = 'HD144432/'                               			# Folder where all the results go
	
	minerals = get_all_solids(keyword, dat, NELEM, NMOLE, NDUST)
	paths = ['Qcurve_inputs_mult_GS/0.1_to_1.5/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_10/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_20/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_100/*.dat']
	gs_ranges = ['0.1_to_1.5', '0.1_to_10', '0.1_to_20', '0.1_to_100']
	
	# Plotting the abunances as a function of radius and temperature
	# R_plot(minerals, dat, keyword, R_arr, R_in, Rmin, Rmax, T0, q, folder, NELEM, NMOLE, NDUST)
	
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
	# not_there = ['SiO', 'Mg3Si4O12H2', 'Fe3Si2O9H4', 'Ni', 'NaAlSi3O8', 'NaMg3AlSi3O12H2', 'CaAl2Si2O8', 'H2O', 'Ca2MgSi2O7', 'NH3', 'Al2O3', 'Ca2Al2SiO7', 'Ca3Al2Si3O12', 'CaAl2Si2O8', 'ZrO2', 'Ti3O5', 'W', 'VO', 'CaTiO3', 'NaAlSiO4']
	# top5_solids = np.setdiff1d(top5_solids, not_there)
	opcurves_available = ['Fe2SiO4', 'Fe3O4', 'Fe', 'FeS', 'Mg2SiO4', 'MgSiO3', 'Ni', 'SiO']
	top5_solids = list(set(top5_solids).intersection(opcurves_available))
	
	# Calculating the surface density
	molwt = molecular_weight(top5_solids)
	surf_dens = surface_density(top5_solids, molwt, topabunds_radii, nHtot, H)
	
	# Defining a nested dictionary where each dictionary within is named after a mineral and each such dictionary consists of grain size ranges as keys with their opacities as values
	kdict = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	lamda = gimme_l('Qcurve_inputs_mult_GS/0.1_to_1.5/Q_Fe2SiO4_rv0.1-1.5_fmax0.7.dat')
	
	for path in paths:	
		
		# Obtaining the grain size range from the given path
		size = path.split('/')[1]
		for opfile in glob.glob(path):
			
			# Obtaining the mineral from the given path
			filename = opfile.split('/')[2]
			mineral = filename.split('_')[1]
			kdict[mineral][size] = gimme_k(opfile)
			
	# Calculating and plotting spectra
	intflux_sum = np.zeros(len(lamda)) * u.Jy
	
	for size in gs_ranges:
		for mineral in top5_solids:
			
			I = Plancks(lamda, Tg)
			tau = tau_calc(surf_dens[mineral], kdict[mineral][size])
			F_map = flux_map(tau, I)
			intflux = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)
			intflux_sum += intflux
		
		# Plot spectra between 8 to 13 microns as in the van Boekel paper
		first = np.where(lamda >= 8 * u.micron)[0][0]
		last = np.where(lamda <= 13 * u.micron)[0][-1]
		lamda_plot = lamda[first: last+1]
		intflux_plot = intflux_sum[first: last+1]
		plt.plot(lamda_plot, intflux_plot * 10**4, label = r'{0} $\mu$m'.format(size))
		
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title(r'{0} Spectrum multiple grain sizes $\mu$m R={1}-{2} AU'.format(disk, Rmin.value, Rmax.value))
	plt.legend()
	plt.savefig(folder + "{0}_spectrum_multgs_R{1}-{2}.png".format(disk, Rmin.value, Rmax.value))
	plt.show()
			
			
if __name__ == "__main__":
	main()
