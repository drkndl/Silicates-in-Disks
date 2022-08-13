import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cmasher as cmr
from astropy import units as u
from astropy.io import fits
from astropy.constants import astropyconst20 as const
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.special import j0
from fancy_name import latex_name
from all_solids import get_all_solids
from molmass import Formula
from diskprop import inner_radius, r_from_T, scale_height, star_radius
from T_plot import T_plot
from radial_plot import R_plot
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from spectra import molecular_weight, surface_density, r_to_rad, slice_lQ, get_l_and_k, Plancks, tau_calc, tau_calc_amorphous, flux_map, calculate_spectra, hankel_transform
from no_thoughts_just_plots import add_textbox, Qcurve_plotter, plot_surf_dens_radial, plot_surf_dens_disk, plot_Bv, plot_tau, plot_fluxmap, plot_spectra
from compare_grain_sizes import get_paper_spectra
from HD144432_q0373_p07_S4900_Tamf1350_mfchange4_gap.properties import *


plt.rcParams["font.family"] = "serif"
plt.rcParams['font.serif'] = ['TeX Gyre Schola']
plt.rcParams['mathtext.fontset'] = 'stix'
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


def baselines(fitsfile):
	
	hdul = fits.open('HD_144432_2022-03-22T09_11_24_N_TARGET_FINALCAL_INT_VISCORRFLUX.fits')
	data = hdul['OI_VIS'].data
	# print(data.columns)
	uc = data.field('UCOORD')
	vc = data.field('VCOORD')
	
	return uc, vc
	
	
def hankel_transform_individual(F_map, R_arr, lamda, B, dist_pc):
	
	"""
	Calculates the absolute correlated flux density (in Jy) using the combined flux map of all the solids. This involves finding the Fourier transform in cylindrical coordinates using the Bessel function of the first kind zeroth order (j0). It is assumed that the correlated fluxes (shape: (lsize,)) are to be plotted against the wavelengths. If, however, it is to be plotted against the baselines, then the correlated flux density (a single float value) corresponding to a specific wavelength (wl) is extracted and returned from the array (shape: (lsize,)). This can be specified with wl_array.
	
	Parameters: 
	
	F_map 		  : 2D array (shape: (NPOINT, lsize)) of the combined flux map for all solids in erg/(s Hz sr cm^2) i.e. CGS units (float)
	R_arr         : 1D array of radii in AU obtained from the temperature array in GGchem output based on the power law model (float)
	lamda         : 1D array of wavelengths (shape: (lsize,)) of the solid in microns (float)
	B 			  : The interferometric baseline value in m (float)
	
	Returns:
	
	inter_flux	  : The absolute correlated flux density in Jy. If wl_array = True, a float array of shape (lsize,) is returned, else a single float value is returned
	"""
	
	rad_arr = r_to_rad(R_arr, dist_pc) 						 											# Converting radius (AU) to arcseconds to radians	
	rad_arr = rad_arr[np.newaxis]
	rad_arr = rad_arr.T 																		# Shape: (NPOINT, 1)
	rad_arr = rad_arr.to('', equivalencies=u.dimensionless_angles()) 
	
	# Finding the Bessel function of the first kind, zeroth order
	bessel = j0(2.0 * np.pi * rad_arr * B / (lamda.to(u.m)))           							# Shape: (NPOINT, lsize)
	
	# Calculating the absolute interferometric flux
	inter_flux = 2.0 * np.pi * np.trapz(rad_arr * bessel * F_map, x = rad_arr, axis = 0) * u.rad**2 		# Shape: (lsize,)
	inter_flux = inter_flux.to(u.Jy, equivalencies = u.dimensionless_angles())
	
	return inter_flux * 10**4
		
		
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
	
	minerals = get_all_solids(keyword, dat, NELEM, NMOLE, NDUST)
	
	# Plotting the abunances as a function of radius and temperature
	# R_plot(minerals, dat, keyword, R_arr, R_in, Rmin, Rmax, T0, q, folder, disk, NELEM, NMOLE, NDUST)
	
	# Finding the most abundant condensates
	abundances, solid_names, abunds_dict = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST) 
	top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names) 
	top5_solids, topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds, abunds_dict)
	
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
		
	# plot_surf_dens_radial(surf_dens, R_arr, folder, **kwargs)
	# plot_surf_dens_disk(surf_dens, R_arr, folder, **kwargs)
	
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
	
	
	################################################################################# Loading the correlated flux data ########################################################################################
	
	hdul = fits.open('HD_144432_2022-03-22T09_11_24_N_TARGET_FINALCAL_INT_VISCORRFLUX.fits')
	data = hdul['OI_VIS'].data
	corrflux = data.field('VISAMP')
	plot_wl = hdul['OI_WAVELENGTH'].data['EFF_WAVE'] * 1e6 * u.micron
	
	uc, vc = baselines('HD_144432_2022-03-22T09_11_24_N_TARGET_FINALCAL_INT_VISCORRFLUX.fits')
	Blines = []
	for i in range(len(uc)):
		Blines.append(np.sqrt(uc[i]**2 + vc[i]**2))		
	Blines = np.round(Blines, 1) * u.m
	
	for i in range(len(corrflux)):	
		plt.plot(plot_wl, corrflux[i], label="%.1f m" % (Blines[i].value))
	
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Correlated flux (Jy)')
	plt.title("HD144432 Baseline Plot - MATISSE")
	plt.legend()
	plt.savefig(folder + "data_corrflux.png")	
	plt.show()
	
	################################################################################# Plotting the data spectrum ########################################################################################
	
	hdul = fits.open('HD_144432_2022-03-22T09_11_24_N_TARGET_FINALCAL_INT_VISCORRFLUX.fits')
	data = hdul['OI_FLUX'].data
	single_dish_spectrum = data.field('FLUXDATA')
	plot_wl = hdul['OI_WAVELENGTH'].data['EFF_WAVE'] * 1e6 * u.micron
	
	plt.plot(plot_wl, single_dish_spectrum[0,:])
	
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title("HD144432 Spectrum - MATISSE")
	plt.legend()
	plt.savefig(folder + "data_spectrum.png")	
	plt.show()
	
	##################################################### Plotting the correlated flux density for multiple baselines against wavelengths #####################################################################
	
	lamda = lamdas['Mg2SiO4']['0.1']
	
	fig, axs = plt.subplots(2, 3, figsize=(20, 10))
	# ~ plt.style.use('dark_background')
	
	axes = [axs[0, 0], axs[0, 1], axs[0, 2], axs[1,0], axs[1, 1], axs[1, 2]]
	n = len(top5_solids)
	# ~ colours = plt.cm.jet(np.linspace(0,1,n))
	# ~ colours = ['darkblue', 'blue',  'teal', 'darkgreen', 'green', 'goldenrod', 'orangered', 'red', 'orchid']
	colours = ['saddlebrown', 'red', 'orange', 'green', 'purple', 'blue', 'gold', 'darkorchid', 'limegreen']
	styles = np.tile(['solid', 'dashed'], n)
	
	for i in range(len(Blines)):
		j, k = 0, 0
		corr_flux_total = np.zeros(lsize) * u.Jy
		for solid in top5_solids:
			for size in gs_ranges:
				
				if solid == 'CaMgSi2O6' and size == '2.0':
					continue
				elif solid == 'Olivine' or solid == 'Pyroxene':
					continue
					
				first = np.where(lamdas[solid][size] >= 8 * u.micron)[0][0]
				last = np.where(lamdas[solid][size] <= 13 * u.micron)[0][-1]
				
				corr_flux = hankel_transform_individual(F_map[solid][size], R_arr, lamdas[solid][size], Blines[i], dist_pc) 
				corr_flux_total += corr_flux
				if Blines[i].value == 62.2:
					axs[0,2].plot(lamdas[solid][size][first: last+1], corr_flux[first: last+1], color=colours[j], linestyle=styles[k], label = '{0} {1}'.format(latex_name(solid), size))
				elif Blines[i].value == 46.3:
					axs[0,0].plot(lamdas[solid][size][first: last+1], corr_flux[first: last+1], color=colours[j], linestyle=styles[k], label = '{0} {1}'.format(latex_name(solid), size))	
				else:
					axes[i].plot(lamdas[solid][size][first: last+1], corr_flux[first: last+1], color=colours[j], linestyle=styles[k], label = '{0} {1}'.format(latex_name(solid), size))	
				k += 1
			j += 1
		
		# Plot the complete model	
		if Blines[i].value == 62.2:
			axs[0,2].plot(lamda[first: last+1], np.abs(corr_flux_total[first: last+1]), color = 'black', label = 'Model')
			axs[0,2].set_ylim([-1.0, 5.0])
			axs[0,2].set_title("B = {0} m".format(Blines[i].value))
		elif Blines[i].value == 46.3:
			axs[0,0].plot(lamda[first: last+1], np.abs(corr_flux_total[first: last+1]), color = 'black', label = 'Model')	
			axs[0,0].set_ylim([-1.0, 5.0])
			axs[0,0].set_title("B = {0} m".format(Blines[i].value))
		else:
			axes[i].plot(lamda[first: last+1], np.abs(corr_flux_total[first: last+1]), color = 'black', label = 'Model')
			axes[i].set_ylim([-1.0, 5.0])
			axes[i].set_title("B = {0} m".format(Blines[i].value))
			
		# ~ axes[i].plot(lamda[first: last+1], np.abs(corr_flux_total[first: last+1]), color = 'black', label = 'Model')
		# ~ axes[i].set_ylim([-1.0, 5.0])
		# ~ axes[i].set_title("B = {0} m".format(Blines[i].value))
		
		# Plot the data
		if Blines[i].value == 62.2:
			axs[0,2].plot(plot_wl, corrflux[i], color='grey', label='MATISSE Data')
		elif Blines[i].value == 46.3:
			axs[0,0].plot(plot_wl, corrflux[i], color='grey', label='MATISSE Data')	
		else:
			axes[i].plot(plot_wl, corrflux[i], color='grey', label='MATISSE Data')
			
		# ~ axes[i].plot(plot_wl, corrflux[i], color='grey', label='MATISSE Data')
		
	textstr = add_textbox(**kwargs)	
	
	for ax in axes:
		ax.set_xlabel(r'$\lambda$ ($\mu$m)')
		ax.set_ylabel('Correlated flux (Jy)')
		ax.label_outer()
		# ax.legend(loc='upper right')
	
	# ~ axes[0].text(0.15, 0.8, textstr, transform=axes[0].transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 13, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))
	
	handles, labels = axes[0].get_legend_handles_labels()
	lgd = fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1), bbox_transform=axes[2].transAxes)
	# fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))
	          
	fig.suptitle(r"Solid-wise Correlated Flux: Multiple Baselines")
	fig.tight_layout()
	plt.savefig(folder + 'Presentation_Correlated_flux_multsolid_multB_multgs.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
	plt.show()


if __name__ == "__main__":
    main()
