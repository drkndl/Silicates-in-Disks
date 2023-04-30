import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cmasher as cmr
from astropy import units as u
from astropy.constants import astropyconst20 as const
from mpl_toolkits.axes_grid1 import make_axes_locatable
from fancy_name import latex_name
from all_solids import get_all_solids
from molmass import Formula
from diskprop import inner_radius, r_from_T, scale_height, star_radius, midplaneT_profile, surface_dens, scale_height, dens
from T_plot import T_plot
from radial_plot import R_plot
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from spectra import molecular_weight, surface_density, r_to_rad, slice_lQ, get_l_and_k, Plancks, tau_calc, tau_calc_amorphous, flux_map, calculate_spectra, hankel_transform
from no_thoughts_just_plots import add_textbox, Qcurve_plotter, plot_surf_dens_radial, plot_surf_dens_disk, plot_Bv, plot_tau, plot_fluxmap, plot_spectra
from compare_grain_sizes import get_paper_spectra
from HD144432_gapring.properties import *


# plt.rcParams["font.family"] = "serif"
# plt.rcParams['font.serif'] = ['TeX Gyre Schola']
# plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['legend.fontsize'] = 13
plt.rcParams['figure.titlesize'] = 16


# Some constants in CGS
Na = const.N_A.cgs                    		# Avogadro's number in /mol
h = const.h.cgs                   			# Planck's constant in cm^2 g s-1
c = const.c.cgs              				# Speed of light in cm/s              
k = const.k_B.cgs                   		# Boltzmann constant in cm^2 g s^-2 K^-1
AU = const.au.cgs                			# 1 astronomical unit in cm
mp = const.m_p.cgs         					# Mass of proton (g)
mu = 2.3                					# Mean molecular weight
Rc = const.R.cgs           					# Ideal gas constant (erg K^-1 mol^-1)
G = const.G.cgs           					# Gravitational constant (cm^3 g^-1 s^-2)	


def main():
	
	kwargs = {'q': q, 'e': e, 'Qr': Qr, 'Sigma0': Sigma0, 'amor_temp': amor_temp, 'add_gap': add_gap, 'rgap': rgap, 'wgap': wgap, 'sgap': sgap, 'add_ring': add_ring, 'rring': rring, 'wring': wring, 'sring': sring}
	
	R_star = star_radius(L_star, T_star).to(u.AU)   		# Star's radius (AU)
	
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
	Rmin = np.round(np.min(R_arr), 1) 						# Minimum radius for spectrum plotting (AU) ENSURE IT IS ONLY 3 DECIMAL PLACES LONG
	Rmax = np.round(np.max(R_arr), 1)						# Maximum radius for spectrum plotting (AU) ENSURE IT IS ONLY 3 DECIMAL PLACES LONG
	# H = scale_height(M_star, R_arr, Tg)

	# Calculating nH theoretical because nH from GGCHEM does not follow the prescribed disk shape; it is log linear
	# R_star = star_radius(L_star, T_star).to(u.AU)   		# Star's radius (AU) 
	# R_in = inner_radius(Qr, T0, R_star, T_star)
	# R_out = r_from_T(R_in, 100.0 * u.K, T0, q)
	# r_arr = np.linspace(R_in, R_out, 500)
	# T_arr = midplaneT_profile(R_in, T0, r_arr, q)	

	# Sigma = surface_dens(Sigma0, r_arr, e)
	# H = scale_height(M_star, r_arr, T_arr)
	# rho, nH = dens(Sigma, H)
	# nHtot = nH

	minerals = get_all_solids(keyword, dat, NELEM, NMOLE, NDUST)
	
	# Plotting the abunances as a function of radius and temperature
	R_plot(minerals, dat, keyword, R_arr, R_in, Rmin, Rmax, T0, q, folder, disk, NELEM, NMOLE, NDUST)
	frejfejrfbejfb
	
	# Finding the most abundant condensates
	abundances, solid_names, abunds_dict = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST) 
	top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names) 
	print(np.shape(abundances))
	print('\n')
	print(top_abunds, np.shape(top_abunds))
	print('\n')
	print(top_solids, np.shape(top_solids)) 
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
	not_there = ['SiO', 'Mg3Si4O12H2', 'Fe3Si2O9H4', 'Ni', 'NaAlSi3O8', 'NaMg3AlSi3O12H2', 'CaAl2Si2O8', 'H2O', 'Ca2MgSi2O7', 'NH3', 'Al2O3', 'Ca2Al2SiO7', 'Ca3Al2Si3O12', 'CaAl2Si2O8', 'ZrO2', 'Ti3O5', 'W', 'VO', 'CaTiO3', 'NaAlSiO4']
	top5_solids = np.setdiff1d(top5_solids, not_there)
	top5_solids = np.concatenate((top5_solids, ['Olivine', 'Pyroxene']))
	
	# Calculating the surface density
	molwt = molecular_weight(top5_solids)
	
	if add_gap and add_ring:
		surf_dens, rgap_ind, wgap_ind1, wgap_ind2, rring_ind, wring_ind1, wring_ind2 = surface_density(top5_solids, molwt, topabunds_radii, nHtot, H, add_gap, R_arr, rgap, wgap, sgap, add_ring, rring, wring, sring)
	elif add_gap:
		surf_dens, rgap_ind, wgap_ind1, wgap_ind2 = surface_density(top5_solids, molwt, topabunds_radii, nHtot, H, add_gap, R_arr, rgap, wgap, sgap, add_ring, rring, wring, sring)
	elif add_ring:
		surf_dens, rring_ind, wring_ind1, wring_ind2 = surface_density(top5_solids, molwt, topabunds_radii, nHtot, H, add_gap, R_arr, rgap, wgap, sgap, add_ring, rring, wring, sring)
	else:
		surf_dens = surface_density(top5_solids, molwt, topabunds_radii, nHtot, H, add_gap, R_arr, rgap, wgap, sgap, add_ring, rring, wring, sring)
		
	plot_surf_dens_radial(surf_dens, R_arr, folder, **kwargs)
	plot_surf_dens_disk(surf_dens, R_arr, folder, **kwargs)
	
	# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
	opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv2.0.dat' : 3.2, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv2.0.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe3O4_rv2.0_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe2SiO4_rv2.0_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_Fe_met_rv2.0_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_FeS_rv2.0_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv2.0_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64, 'Qcurve_inputs/qval_MgAl2O4_rv2.0_fmax0.7.dat' : 3.64, 'Qcurve_inputs/Q_Olivine_rv0.1_fmax0.7.dat': 3.71, 'Qcurve_inputs/Q_Olivine_rv2.0_fmax0.7.dat' : 3.71, 'Qcurve_inputs/Q_Pyroxene_rv0.1_fmax0.7.dat': 3.01, 'Qcurve_inputs/Q_Pyroxene_fmax0.7_rv2.0.dat': 3.01}
	
	# Adding units to the material densities using astropy
	for key, value in opfile_dens.items():
		opfile_dens[key] = value * u.g / u.cm**3
	
	# Initializing dictionaries for the wavelengths, opacities, grain sizes and emptiness fractions used in the opfiles for each solid 	
	gs_ranges = ['0.1', '2.0']
	lamdas = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	kappas = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	fmaxs = {key: None for key in top5_solids}
	
	for opfile, density in opfile_dens.items():
		mineral, rv, fmax, lamda, kappa = get_l_and_k(opfile, density, lmin, lmax, lsize)		
		# Qcurve_plotter(lamda, kappa, mineral, rv, fmax, folder)
		lamdas[mineral][rv] = lamda
		kappas[mineral][rv] = kappa
		fmaxs[mineral] = fmax			
	
	# Plotting the flux map and calculating the integrated flux for each solid
	I = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	tau = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	F_map = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	F_map_sum = np.zeros((NPOINT, lsize)) * u.erg / (u.s * u.Hz * u.sr * u.cm**2)
	intflux = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	intflux_sum = np.zeros(lsize) * u.Jy
	
	for solid in top5_solids:
		for size in gs_ranges:
			
			if lamdas[solid][size] == None:
				print("Skipping: ", solid, size)
				continue
		
			if solid == 'Olivine':
				continue
			
			elif solid == 'Pyroxene':
				continue
						
			elif solid == "Mg2SiO4":
					
				I[solid][size] = Plancks(lamdas[solid][size], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid][size] = tau_calc_amorphous(surf_dens[solid], surf_dens['Olivine'], kappas[solid][size], kappas['Olivine'][size], Tg, amor_temp, mass_fracs[solid][size])
				# plot_tau(tau[solid][size], solid, size, folder)
				
				F_map[solid][size] = flux_map(tau[solid][size], I[solid][size])
				# plot_fluxmap(solid, size, fmaxs[solid], F_map[solid][size], lamdas[solid][size], R_arr, folder, **kwargs) 
				F_map_sum += F_map[solid][size]
				
				intflux[solid][size] = calculate_spectra(F_map[solid][size], R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid][size], intflux[solid][size], solid, size, fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux[solid][size]
				
			elif solid == "MgSiO3":
					
				I[solid][size] = Plancks(lamdas[solid][size], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid][size] = tau_calc_amorphous(surf_dens[solid], surf_dens['Pyroxene'], kappas[solid][size], kappas['Pyroxene'][size], Tg, amor_temp, mass_fracs[solid][size])
				# plot_tau(tau[solid][size], solid, size, folder)
				
				F_map[solid][size] = flux_map(tau[solid][size], I[solid][size])
				# plot_fluxmap(solid, size, fmaxs[solid], F_map[solid][size], lamdas[solid][size], R_arr, folder, **kwargs) 
				F_map_sum += F_map[solid][size]
				
				intflux[solid][size] = calculate_spectra(F_map[solid][size], R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid][size], intflux[solid][size], solid, size, fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux[solid][size]
				
			else:
				
				I[solid][size] = Plancks(lamdas[solid][size], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid][size] = tau_calc(surf_dens[solid], kappas[solid][size], mass_fracs[solid][size])
				# plot_tau(tau[solid][size], solid, size, folder)
				
				F_map[solid][size] = flux_map(tau[solid][size], I[solid][size])
				# plot_fluxmap(solid, size, fmaxs[solid], F_map[solid][size], lamdas[solid][size], R_arr, folder, **kwargs)
				F_map_sum += F_map[solid][size]
				
				intflux[solid][size] = calculate_spectra(F_map[solid][size], R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid][size], intflux[solid][size], solid, size, fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux[solid][size]
			
	########################################################################## Plotting the overall flux map ################################################################################################
	fig, ax = plt.subplots(1,1)
	img = ax.imshow(F_map_sum, cmap='cmr.gem', interpolation='none')
	
	x_axis_locations = np.linspace(0, lsize-1, 8).astype(int)
	ax.set_xticks(x_axis_locations)
	x_axis_labels = np.round(lamda[x_axis_locations], 1)
	ax.set_xticklabels(x_axis_labels.value)
	ax.set_xlabel(r'$\lambda$ ($\mu$m)')
	
	y_axis_locations = np.linspace(0, len(R_arr)-1, 8).astype(int) 
	ax.set_yticks(y_axis_locations)
	y_axis_labels = np.round(R_arr[y_axis_locations], 1) 
	ax.set_yticklabels(y_axis_labels.value)
	ax.set_ylabel('R (AU)')
	
	ax.set_title("Combined Flux Map for all Solids")
	
	textstr = add_textbox(**kwargs)	
	plt.text(0.85, 0.25, textstr, transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 10, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	fig.colorbar(img, label=r'erg / (s Hz sr $cm^2$)', cax=cax)
	plt.tight_layout()
	plt.savefig(folder + "overall_fluxmap.png", bbox_inches = 'tight')
	plt.show()
	
	######################################################## Plotting the solid spectra, overall spectrum and the paper spectrum #############################################################################
	lamda = lamdas['Mg2SiO4']['0.1']
	diskname = disk.split('_')[0]
	# datfile = folder + 'van_Boekel_' + diskname + '.dat'
	datfile = folder + 'van_Boekel_HD144432.dat'
	wl, flux = get_paper_spectra(datfile)
	
	n = len(top5_solids)
	# ~ fig, axs = plt.subplots(1, 2, figsize=(14, 7))
	# ~ # colours = cmr.take_cmap_colors('cmr.wildfire', n, cmap_range=(0.1, 0.9))
	# ~ colours = plt.cm.jet(np.linspace(0,1,n))
	# ~ styles = np.tile(['solid', 'dashed'], n)
	# ~ i, j = 0, 0
	
	# ~ for solid in top5_solids:
		# ~ for size in gs_ranges:
			
			# ~ if lamdas[solid][size] == None or solid=="Olivine" or solid=="Pyroxene":
				# ~ print("Skipping: ", solid, size)
				# ~ continue
				
			# ~ axs[0].plot(lamdas[solid][size], intflux[solid][size], color=colours[i], linestyle=styles[j], label=r"{0} {1}$\mu$m".format(latex_name(solid), size))
			# ~ j += 1
		# ~ i += 1
		
	# ~ axs[0].plot(lamda, intflux_sum, color="black", label='Model')
	# ~ axs[0].plot(wl, flux, color="grey", label="Data")
	
	# ~ # 8 to 13 microns zoomed in version of the above plot 
	# ~ i, j = 0, 0

	# ~ for solid in top5_solids:
		# ~ for size in gs_ranges:
			
			# ~ if lamdas[solid][size] == None or solid=="Olivine" or solid=="Pyroxene":
				# ~ print("Skipping: ", solid, size)
				# ~ continue
			
			# ~ first = np.where(lamdas[solid][size] >= 8 * u.micron)[0][0]
			# ~ last = np.where(lamdas[solid][size] <= 13 * u.micron)[0][-1]		
			# ~ axs[1].plot(lamdas[solid][size][first: last+1], intflux[solid][size][first: last+1], color=colours[i], linestyle=styles[j], label=r"{0} {1}$\mu$m".format(latex_name(solid), size))
			# ~ j += 1
		# ~ i += 1
			
	# ~ axs[1].plot(lamda[first: last+1], intflux_sum[first: last+1], color="black", label='Model')
	# ~ axs[1].plot(wl, flux, color="grey", label="Data")
	
	# ~ textstr = add_textbox(**kwargs)
	
	# ~ for ax in axs:
		# ~ ax.set_xlabel(r'$\lambda$ ($\mu$m)')
		# ~ ax.set_ylabel('Flux (Jy)')
		# ~ ax.label_outer()
		
	# ~ axs[0].text(0.15, 0.7, textstr, transform=axs[0].transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 10, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))		
	# ~ axs[0].set_title(r"0 to 20 $\mu$m")
	# ~ axs[1].set_title(r"8 to 13 $\mu$m")
	# ~ fig.suptitle("Combined Spectrum for all Solids")
	# ~ axs[1].legend(loc='upper left', bbox_to_anchor=(1, 1))
	
	# ~ fig.tight_layout()
	# ~ plt.savefig(folder + "Overall_spectrum_multgs_R{0}-{1}.png".format(Rmin.value, Rmax.value))
	# ~ plt.show()	
	
	fig, axs = plt.subplots(1, 1, figsize=(12,10))
	# colours = cmr.take_cmap_colors('cmr.wildfire', n, cmap_range=(0.1, 0.9))
	# ~ colours = plt.cm.jet(np.linspace(0,1,n))
	# ~ colours = cmr.take_cmap_colors('cmr.tropical', n, cmap_range=(0.0, 1.0))
	colours = ['saddlebrown', 'red', 'orange', 'green', 'purple', 'blue', 'gold', 'darkorchid', 'limegreen']
	styles = np.tile(['solid', 'dashed'], n)
	i, j = 0, 0

	for solid in top5_solids:
		for size in gs_ranges:
			
			if lamdas[solid][size] == None or solid=="Olivine" or solid=="Pyroxene":
				print("Skipping: ", solid, size)
				continue
			
			first = np.where(lamdas[solid][size] >= 8 * u.micron)[0][0]
			last = np.where(lamdas[solid][size] <= 13 * u.micron)[0][-1]		
			axs.plot(lamdas[solid][size][first: last+1], intflux[solid][size][first: last+1], color=colours[i], linestyle=styles[j], label=r"{0} {1}$\mu$m".format(latex_name(solid), size))
			j += 1
		i += 1
			
	axs.plot(lamda[first: last+1], intflux_sum[first: last+1], color="black", label='Model')
	axs.plot(wl, flux, color="grey", label="MATISSE Data")
	
	textstr = add_textbox(**kwargs)
	
	axs.set_xlabel(r'$\lambda$ ($\mu$m)')
	axs.set_ylabel('Flux (Jy)')
		# ~ ax.label_outer()
		
	axs.text(0.15, 0.85, textstr, transform=axs.transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 13, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))		
	fig.suptitle("Single Dish Spectrum with Solid Contributions: Data vs Model")
	axs.legend(loc='upper left', bbox_to_anchor=(1, 1))
	
	fig.tight_layout()
	plt.savefig(folder + "Presentation_Overall_spectrum_multgs_R{0}-{1}.png".format(Rmin.value, Rmax.value))
	plt.show()

	
	################################################################### Plotting the overall spectrum before and after the gap/ring ###########################################################################
	
	if add_gap:
		
		R_gap1 = np.round(R_arr[wgap_ind1], 1)
		R_gap2 = np.round(R_arr[wgap_ind2], 1)
		
		# Calculating and plotting the spectrum before the gap
		intflux_beforegap = calculate_spectra(F_map_sum, R_arr, Rmin, R_gap1, dist_pc)	
		plt.plot(lamda, intflux_beforegap, color='blue', label=r"Before Gap = ({0},{1}) AU".format(Rmin.value, R_gap1.value))
		
		# Calculating and plotting the spectrum in the gap
		# ~ intflux_gap = calculate_spectra(F_map_sum, R_arr, R_gap1, R_gap2, dist_pc)	
		# ~ plt.plot(lamda, intflux_gap, label=r"Gap = ({0},{1}) AU".format(R_gap1.value, R_gap2.value))
		
		# Calculating and plotting the spectrum after the gap
		intflux_aftergap = calculate_spectra(F_map_sum, R_arr, R_gap2, Rmax, dist_pc)	
		plt.plot(lamda, intflux_aftergap, color='red', linestyle='dashed', label=r"After Gap = ({0},{1}) AU".format(R_gap2.value, Rmax.value))
		
		plt.xlabel(r'$\lambda$ ($\mu$m)')
		plt.ylabel('Flux (Jy)')
		
		plt.title(r'Single Dish Spectra Before and After Gap')
		plt.legend(loc='upper right')	
		plt.tight_layout()
		plt.savefig(folder + "Presentation_Gap_spectra.png")
		plt.show()	
		
	if add_ring:
		
		R_ring1 = np.round(R_arr[wring_ind1], 1)
		R_ring2 = np.round(R_arr[wring_ind2], 1)
		
		# Calculating and plotting the spectrum before the ring
		intflux_beforering = calculate_spectra(F_map_sum, R_arr, Rmin, R_ring1, dist_pc)	
		plt.plot(lamda, intflux_beforering, label=r"Before Ring = ({0},{1}) AU".format(Rmin.value, R_ring1.value))
		
		# Calculating and plotting the spectrum in the ring
		intflux_ring = calculate_spectra(F_map_sum, R_arr, R_ring1, R_ring2, dist_pc)	
		plt.plot(lamda, intflux_ring, label=r"Ring = ({0},{1}) AU".format(R_ring1.value, R_ring2.value))
		
		# Calculating and plotting the spectrum after the ring
		intflux_afterring = calculate_spectra(F_map_sum, R_arr, R_ring2, Rmax, dist_pc)	
		plt.plot(lamda, intflux_afterring, label=r"After Ring = ({0},{1}) AU".format(R_ring2.value, Rmax.value))
		
		plt.xlabel(r'$\lambda$ ($\mu$m)')
		plt.ylabel('Flux (Jy)')
		
		plt.title(r'Ring Spectra')
		plt.legend(loc='upper right')	
		plt.tight_layout()
		plt.savefig(folder + "Ring_spectra.png")
		plt.show()		
		
	############################################################## Plotting the overall spectrum considering multiple radii limits ###########################################################################
	R_list = np.round(R_arr[np.linspace(0, len(R_arr)-1, 7).astype(int)], 1)
	Rmin_list = R_list[:-1]
	Rmax_list = R_list[1:]
	fig, ax = plt.subplots(1, 1)
	
	n = len(Rmin_list)
	colours = cmr.take_cmap_colors('cmr.guppy', n, cmap_range=(0.1, 0.9))
	j = 0
	
	for i in range(len(Rmax_list)):			
		intflux_sum_mr = calculate_spectra(F_map_sum, R_arr, Rmin_list[i], Rmax_list[i], dist_pc)	
		plt.plot(lamda, intflux_sum_mr, color = colours[j], label=r"($R_{{min}}$,$R_{{max}}$) = ({0},{1}) AU".format(Rmin_list[i].value, Rmax_list[i].value))
		j += 1
		
	textstr = add_textbox(**kwargs)	
	plt.text(0.15, 0.7, textstr, transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 10, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))
	
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title(r'Spectra for Different Radii Limits')
	plt.legend(loc='upper right')	
	plt.savefig(folder + "Overall_spectrum_multiple_radii_limits.png")
	plt.show()
	
	##################################################### Plotting the correlated flux density for multiple baselines against wavelengths #####################################################################
	fig, axs = plt.subplots(1, 2, figsize=(14, 7))
	n = len(B_small)
	colours = plt.cm.viridis(np.linspace(0,1,n))
	i = 0
	
	for Bl in B_small:
		corr_flux = np.zeros(lsize) * u.Jy
		for solid in top5_solids:
			for size in gs_ranges:
				
				if solid == 'CaMgSi2O6' and size == '2.0':
					continue
				elif solid == 'Olivine' or solid == 'Pyroxene':
					continue
					
				corr_flux += np.abs(hankel_transform(F_map[solid][size], R_arr, lamdas[solid][size], 0.0, Bl, dist_pc, wl_array = True)) 
		axs[0].plot(lamda, corr_flux, color = colours[i], label = 'B = {0} m'.format(Bl.value))
		axs[1].plot(lamda[first: last+1], corr_flux[first: last+1], color = colours[i], label = 'B = {0} m'.format(Bl.value))	
		i += 1
	
	textstr = add_textbox(**kwargs)	
	
	for ax in axs:
		ax.set_xlabel(r'$\lambda$ ($\mu$m)')
		ax.set_ylabel('Correlated flux (Jy)')
		ax.label_outer()
		ax.legend(loc='upper right')
	
	axs[0].text(0.15, 0.7, textstr, transform=axs[0].transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 10, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))
	axs[0].set_title(r"0 to 20 $\mu$m")
	axs[1].set_title(r"8 to 13 $\mu$m")
		             
	fig.suptitle(r"Correlated Flux: Multiple Baselines")
	fig.tight_layout()
	plt.savefig(folder + 'Correlated_flux_multB_multgs_other.png')
	plt.show()
	
	##################################################### Plotting the correlated flux density for a single wavelength against baselines ######################################################################
	fig, axs = plt.subplots(2,1)
	inter_flux = np.zeros(len(B)) * u.Jy
	n = len(wl_list)
	colours = plt.cm.viridis(np.linspace(0,1,n))
	i = 0
	
	for wavel in wl_list:
		print(wavel)
		for bl in range(len(B)):			
			inter_flux[bl] = hankel_transform(F_map_sum, R_arr, lamda, wavel, B[bl], dist_pc, wl_array = False) 
		axs[0].plot(B, inter_flux, color = colours[i], label=r"{0} $\mu$m".format(wavel.value))
		if np.max(inter_flux) < 2.0 * u.Jy:
			axs[1].plot(B, inter_flux, color = colours[i], label=r"{0} $\mu$m".format(wavel.value))	
		inter_flux = np.zeros(len(B)) * u.Jy	
		i += 1
	
	for ax in axs:
		ax.set_xlabel(r'Baseline B (m)')
		ax.set_ylabel('Correlated flux (Jy)')
		ax.label_outer()
		ax.legend(loc='upper right')
		
	axs[0].text(0.7, 0.6, textstr, transform=axs[0].transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 10, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))
	axs[0].set_title(r"All wavelengths")
	axs[1].set_title(r"Zoomed-in version")
	
	fig.suptitle(r"Correlated Flux: Multiple Wavelengths")
	fig.tight_layout()
	plt.savefig(folder + "Correlated_flux_multWL_multgs_other.png")
	plt.show()


if __name__ == "__main__":
    main()
