# Attempts to recreate disk spectra from van Boekel (2005)

import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.constants import astropyconst20 as const
from scipy.special import j0
from fancy_name import latex_name
from all_solids import get_all_solids
from molmass import Formula
from diskprop import inner_radius, r_from_T, scale_height, star_radius
from T_plot import T_plot
from radial_plot import R_plot
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from spectra import molecular_weight, surface_density, r_to_rad, slice_lQ, get_l_and_k, Plancks, tau_calc, tau_calc_amorphous, flux_map, calculate_spectra, hankel_transform
from no_thoughts_just_plots import Qcurve_plotter, plot_Bv, plot_tau, plot_fluxmap, plot_spectra
from compare_Qcurve import gimme_k, gimme_l
from HD179218.properties import *


def get_paper_spectra(filename):
	
	"""
	Data points from the spectra in van Boekel (2005) are obtained using the WebPlotDigitizer software. The output of the software is a csv file with 2 columns: wavelength (microns) and flux (Jy). This function obtains the column data from the file. 
	
	Parameters:
	
	filename 	  : CSV file with data from spectra data from van Boekel (2005) with 2 columns: wavelength (microns) and flux (Jy)
	
	Returns:
	
	wl 			  : 1D array of wavelengths from the spectra images in van Boekel (2005) (microns)
	flux 		  : 1D array of flux from the spectra images in van Boekel (2005) (Jy)
	"""
	
	spectrum = pd.read_csv(filename, delimiter=',', names = ['wavelength', 'flux'])
	wl = u.Quantity(spectrum['wavelength'].to_numpy(), u.micron)
	flux = spectrum['flux'].to_numpy() * u.Jy
	
	return wl, flux

	
	
def main():
	
	subplot_size = 12                                        # Some font size definitions
	fullplot_size = 15
	labelsize = 11
	
	R_star = star_radius(L_star, T_star).to(u.AU)   		# Star's radius (AU)
	fmax = 0.7
	
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
	
	R_in = inner_radius(Qr, T0, R_star, T_star)   # Inner-most radius beyond which the dust is sublimated (AU)
	R_arr = r_from_T(R_in, Tg, T0, q)             # 1D array of radii obtained from the power law disk model (AU)
	Rmin = np.round(np.min(R_arr), 3) 						# Minimum radius for spectrum plotting (AU) ENSURE IT IS ONLY 3 DECIMAL PLACES LONG
	Rmax = np.round(np.max(R_arr), 3)						# Maximum radius for spectrum plotting (AU) ENSURE IT IS ONLY 3 DECIMAL PLACES LONG
	# H = scale_height(M_star, R_arr, Tg)
	
	minerals = get_all_solids(keyword, dat, NELEM, NMOLE, NDUST)
	paths = ['Qcurve_inputs_mult_GS/0.1_to_1.5/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_10/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_20/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_100/*.dat']
	gs_ranges = ['0.1_to_1.5', '0.1_to_10', '0.1_to_20', '0.1_to_100']
	
	# Plotting the abunances as a function of radius and temperature
	R_plot(minerals, dat, keyword, R_arr, R_in, Rmin, Rmax, T0, q, folder, disk, NELEM, NMOLE, NDUST)
	
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
	lamda = gimme_l('Qcurve_inputs_mult_GS/0.1_to_1.5/Q_Mg2SiO4_rv0.1-1.5_fmax0.7.dat')         # Since wavelength values for all solids are the same, just picking one file
	
	for path in paths:	
		
		# Obtaining the grain size range from the given path
		size = path.split('/')[1]
		for opfile in glob.glob(path):
			
			# Obtaining the mineral from the given path
			filename = opfile.split('/')[2]
			mineral = filename.split('_')[1]
			kdict[mineral][size] = gimme_k(opfile)
			# Qcurve_plotter(lamda, kdict[mineral][size], mineral, size, fmax, folder)
			
	#################################################################### Calculating and plotting spectra ###################################################################################################
	intflux = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	F_map_sum = {key: np.zeros((NPOINT, len(lamda))) * u.erg / (u.s * u.Hz * u.sr * u.cm**2) for key in gs_ranges}
	intflux_sum = {key: np.zeros(len(lamda)) * u.Jy for key in gs_ranges}
	
	for size in gs_ranges:
		for mineral in top5_solids:
			
			I = Plancks(lamda, Tg)
			tau = tau_calc(surf_dens[mineral], kdict[mineral][size])
			F_map = flux_map(tau, I)
			F_map_sum[size] += F_map
			intflux[mineral][size] = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)
			intflux_sum[size] += intflux[mineral][size]
		
		# Plot overall spectra between 8 to 13 microns as in the van Boekel paper for each grain size range
		first = np.where(lamda >= 8 * u.micron)[0][0]
		last = np.where(lamda <= 13 * u.micron)[0][-1]
		lamda_plot = lamda[first: last+1]
		intflux_plot = intflux_sum[size][first: last+1]
		plt.plot(lamda_plot, intflux_plot, label = r'{0} $\mu$m'.format(size))
		
	# Plotting paper spectrum for the given disk
	datfile = folder + 'van_Boekel_' + disk + '.dat'
	wl, flux = get_paper_spectra(datfile)
	plt.plot(wl, flux, label="van Boekel (2005)")
		
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title(r'{0} Spectrum multiple grain sizes R={1}-{2} AU'.format(disk, Rmin.value, Rmax.value))
	plt.legend()
	plt.savefig(folder + "{0}_spectrum_multgs_R{1}-{2}.png".format(disk, Rmin.value, Rmax.value))
	plt.show()
	
	#################################################### Plotting the spectra of each individual condensate for the various grain size ranges ############################################################
	
	# Obtaining indices for only 0-20 micron wavelength values
	ind_20 = np.where(lamda <= 20 * u.micron)[0][-1]
	
	for mineral in top5_solids:
		for size in gs_ranges:
			plt.plot(lamda[:ind_20+1], intflux[mineral][size][:ind_20+1], label=size)
		plt.xlabel(r'$\lambda$ ($\mu$m)')
		plt.ylabel('Flux (Jy)')
		plt.title(r'{0} {1} Overall Spectrum multiple grain sizes'.format(disk, latex_name(mineral)))
		plt.legend()
		plt.savefig(folder + "{0}_{1}_spectrum_multgs.png".format(disk, mineral))
		plt.show()
		
	##################################################################### Plotting the spectra for multiple radii limits ##################################################################################
	fig, axs = plt.subplots(2, 2, figsize=(20, 15))
	axes = [axs[0, 0], axs[0, 1], axs[1,0], axs[1, 1]]
	i = 0
	R_list = np.round(R_arr[np.linspace(0, len(R_arr)-1, 7).astype(int)], 3)
	Rmin_list = R_list[:-1]
	Rmax_list = R_list[1:]
	
	for size in gs_ranges:
		for j in range(len(Rmax_list)):		
				
			intflux_sum_mr = calculate_spectra(F_map_sum[size], R_arr, Rmin_list[j], Rmax_list[j], dist_pc)	
			axes[i].plot(lamda[:ind_20+1], intflux_sum_mr[:ind_20+1], label=r"($R_{{min}}$,$R_{{max}}$) = ({0},{1}) AU".format(Rmin_list[j].value, Rmax_list[j].value))
		
		axes[i].set_xlabel(r'$\lambda$ ($\mu$m)', fontsize = labelsize)
		axes[i].set_ylabel('Flux (Jy)', fontsize = labelsize)
		axes[i].set_title("gs={0} $\mu$m".format(size), fontsize = subplot_size) 
		axes[i].legend()
		i += 1
	
	axs[0,0].get_xaxis().set_visible(False)
	axs[0,1].get_xaxis().set_visible(False)
             
	fig.suptitle(r'{0} Spectra for multiple radial limits and grain sizes'.format(disk) , fontsize = fullplot_size)
	fig.tight_layout()
	plt.savefig(folder + 'spectra_multradii_multgs.png')
	plt.show()		
		
	##################################################### Plotting the correlated flux density for multiple baselines against wavelengths ##################################################################
	fig, axs = plt.subplots(2, 2, figsize=(20, 15))
	axes = [axs[0, 0], axs[0, 1], axs[1,0], axs[1, 1]]
	i = 0
	
	for size in gs_ranges:
		
		for Bl in B_small:
			
			corr_flux_absB = hankel_transform(F_map_sum[size], R_arr, lamda, wl, Bl, dist_pc, wl_array = True) 
			axes[i].plot(lamda[:ind_20+1], corr_flux_absB[:ind_20+1], label = r'B={0} m'.format(Bl.value))
		
		axes[i].set_xlabel(r'$\lambda$ ($\mu$m)', fontsize = labelsize)
		axes[i].set_ylabel('Correlated flux (Jy)', fontsize = labelsize)	
		axes[i].set_title("gs={0} $\mu$m".format(size), fontsize = subplot_size) 
		axes[i].legend()
		i += 1
	
	axs[0,0].get_xaxis().set_visible(False)
	axs[0,1].get_xaxis().set_visible(False)
             
	fig.suptitle(r'{0} correlated flux for multiple baselines and grain sizes'.format(disk), fontsize = fullplot_size)
	fig.tight_layout()
	plt.savefig(folder + 'Correlated_flux_multB_multgs.png')
	plt.show()
	
	############################################## Plotting correlated fluxes against baselines for multiple wavelengths and size ranges ############################################################3
	inter_flux = {key: np.zeros(len(B)) * u.Jy for key in gs_ranges}
	fig, axs = plt.subplots(2, 2, figsize=(20, 15))
	axes = [axs[0, 0], axs[0, 1], axs[1,0], axs[1, 1]]
	i = 0
	
	for size in gs_ranges:		
		for wl in wl_list:
			for bl in range(len(B)):			
				inter_flux[size][bl] = hankel_transform(F_map_sum[size], R_arr, lamda, wl, B[bl], dist_pc, wl_array = False)
			axes[i].plot(B, inter_flux[size], label=r"{0} $\mu$m".format(wl.value))	
		
		axes[i].set_xlabel(r'Baseline (m)', fontsize = labelsize)
		axes[i].set_ylabel('Correlated flux (Jy)', fontsize = labelsize)
		axes[i].set_title("gs={0} $\mu$m".format(size), fontsize = subplot_size) 
		axes[i].legend()
		i += 1
             
	fig.suptitle(r'{0} correlated flux for multiple wavelengths and grain sizes'.format(disk), fontsize = fullplot_size)
	fig.tight_layout()
	plt.savefig(folder + 'Correlated_flux_multwl_multgs.png')
	plt.show()
				
if __name__ == "__main__":
	main()
