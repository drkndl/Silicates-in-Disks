import numpy as np
import matplotlib.pyplot as plt
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
from HD163296_q04_p075_S3400.properties import *

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
	Rmin = np.round(np.min(R_arr), 3) 						# Minimum radius for spectrum plotting (AU) ENSURE IT IS ONLY 3 DECIMAL PLACES LONG
	Rmax = np.round(np.max(R_arr), 3)						# Maximum radius for spectrum plotting (AU) ENSURE IT IS ONLY 3 DECIMAL PLACES LONG
	# H = scale_height(M_star, R_arr, Tg)
	
	minerals = get_all_solids(keyword, dat, NELEM, NMOLE, NDUST)
	
	# Plotting the abunances as a function of radius and temperature
	# R_plot(minerals, dat, keyword, R_arr, R_in, Rmin, Rmax, T0, q, folder, NELEM, NMOLE, NDUST)
	
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
	surf_dens = surface_density(top5_solids, molwt, topabunds_radii, nHtot, H, add_gap, R_arr, rgap, wgap, sgap) 
	
	# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
	opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64, 'Qcurve_inputs/Q_Olivine_rv0.1_fmax0.7.dat': 3.71, 'Qcurve_inputs/Q_Pyroxene_rv0.1_fmax0.7.dat': 3.01}
	
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
	
	for amor_temp in amor_temp_list:
		
		# Plotting the flux map and calculating the integrated flux for each solid
		I = {key: None for key in top5_solids}
		tau = {key: None for key in top5_solids}
		F_map_sum = np.zeros((NPOINT, lsize)) * u.erg / (u.s * u.Hz * u.sr * u.cm**2)
		intflux_sum = np.zeros(lsize) * u.Jy
		
		for solid in top5_solids:
			
			if solid == 'Olivine':
				continue
			
			elif solid == 'Pyroxene':
				continue
						
			elif solid == "Mg2SiO4":
				
				I[solid] = Plancks(lamdas[solid], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid] = tau_calc_amorphous(surf_dens[solid], surf_dens['Olivine'], kappas[solid], kappas['Olivine'], Tg, amor_temp)
				# plot_tau(tau[solid], solid, folder)
				
				F_map = flux_map(tau[solid], I[solid])
				# plot_fluxmap(solid, rvs[solid], fmaxs[solid], F_map, lamdas[solid], R_arr, folder) 
				F_map_sum += F_map
				
				intflux = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid], intflux, solid, rvs[solid], fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux
				
			elif solid == "MgSiO3":
				
				I[solid] = Plancks(lamdas[solid], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid] = tau_calc_amorphous(surf_dens[solid], surf_dens['Pyroxene'], kappas[solid], kappas['Pyroxene'], Tg, amor_temp)
				# plot_tau(tau[solid], solid, folder)
				
				F_map = flux_map(tau[solid], I[solid])
				# plot_fluxmap(solid, rvs[solid], fmaxs[solid], F_map, lamdas[solid], R_arr, folder) 
				F_map_sum += F_map
				
				intflux = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid], intflux, solid, rvs[solid], fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux
				
			else:
				
				I[solid] = Plancks(lamdas[solid], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid] = tau_calc(surf_dens[solid], kappas[solid])
				# plot_tau(tau[solid], solid, folder)
				
				F_map = flux_map(tau[solid], I[solid])
				# plot_fluxmap(solid, rvs[solid], fmaxs[solid], F_map, lamdas[solid], R_arr, folder) 
				F_map_sum += F_map
				
				intflux = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid], intflux, solid, rvs[solid], fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux
				
		# ~ # Plotting the overall flux map
		# ~ fig, ax = plt.subplots(1,1)
		# ~ img = ax.imshow(F_map_sum, cmap='plasma', interpolation='none')
		
		# ~ x_axis_locations = np.linspace(0, lsize-1, 8).astype(int)
		# ~ ax.set_xticks(x_axis_locations)
		# ~ x_axis_labels = np.round(lamda[x_axis_locations], 1)
		# ~ ax.set_xticklabels(x_axis_labels.value)
		# ~ ax.set_xlabel(r'$\lambda$ ($\mu$m)')
		
		# ~ y_axis_locations = np.linspace(0, len(R_arr)-1, 8).astype(int) 
		# ~ ax.set_yticks(y_axis_locations)
		# ~ y_axis_labels = np.round(R_arr[y_axis_locations], 3) 
		# ~ ax.set_yticklabels(y_axis_labels.value)
		# ~ ax.set_ylabel('R (AU)')
		
		# ~ ax.set_title(r"{0} Overall Flux Map for r=0.1 microns $T_{{am}}$ = {1}K".format(disk, amor_temp.value))
		# ~ fig.colorbar(img, label=r'{0}'.format(F_map_sum.unit))
		# ~ plt.savefig(folder + "overall_fluxmap_T{0}.png".format(amor_temp.value), bbox_inches = 'tight')
		# ~ plt.show()
		
		# Plotting the overall spectrum
		# fig = plt.figure()
		plt.plot(lamdas['Mg2SiO4'], intflux_sum, label=r'T={0}'.format(amor_temp))
	
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.legend()
	plt.title(r'{0} Overall Spectrum r=0.1 $\mu$m'.format(disk))
	plt.savefig("Overall_spectrum_r0.1.png")
	plt.show()
		
	# ~ printfffff
		# ~ # Plotting the overall spectrum considering multiple radii limits
		# ~ R_list = np.round(R_arr[np.linspace(0, len(R_arr)-1, 7).astype(int)], 3)
		# ~ Rmin_list = R_list[:-1]
		# ~ Rmax_list = R_list[1:]
		# ~ fig = plt.figure()
		
		# ~ for i in range(len(Rmax_list)):			 			
			# ~ intflux_sum_mr = calculate_spectra(F_map_sum, R_arr, Rmin_list[i], Rmax_list[i], dist_pc)			
			# ~ plt.plot(lamdas['Mg2SiO4'], intflux_sum_mr, label=r"($R_{{min}}$,$R_{{max}}$) = ({0},{1}) AU".format(Rmin_list[i].value, Rmax_list[i].value))
				
		# ~ plt.xlabel(r'$\lambda$ ($\mu$m)')
		# ~ plt.ylabel('Flux (Jy)')
		# ~ plt.title(r'{0} Overall spectrum for multiple radii $T_{{am}}$ = {1}K'.format(disk, amor_temp.value))
		# ~ plt.legend()	
		# ~ plt.savefig(folder + "Overall_spectrum_multiple_radii_limits_T{0}.png".format(amor_temp.value))
		# ~ plt.show()
		
		# ~ # Plotting the correlated flux density for multiple baselines against wavelengths
		# ~ fig = plt.figure()
		
		# ~ for Bl in B_small:
			
			# ~ corr_flux_absB = hankel_transform(F_map_sum, R_arr, lamdas['Mg2SiO4'], wl, Bl, dist_pc, wl_array = True) 
			# ~ plt.plot(lamdas['Mg2SiO4'], corr_flux_absB, label = 'B = {0} m'.format(Bl.value))
		
		# ~ plt.xlabel(r'$\lambda$ ($\mu$m)')
		# ~ plt.ylabel('Correlated flux (Jy)')             
		# ~ plt.title(r'{0} Correlated flux for multiple baselines, rv = 0.1, $T_{{am}}$ = {1}K'.format(disk, amor_temp.value))
		# ~ plt.legend()
		# ~ plt.savefig(folder + 'Correlated_flux_multB_rv0.1_T{0}.png'.format(amor_temp.value))
		# ~ plt.show()
		
		# ~ # Plotting the correlated flux density for a single wavelength against baselines
		# ~ inter_flux = np.zeros(len(B)) * u.Jy
		
		# ~ for wl in wl_list:
			# ~ for bl in range(len(B)):			
				# ~ inter_flux[bl] = hankel_transform(F_map_sum, R_arr, lamdas['Mg2SiO4'], wl, B[bl], dist_pc, wl_array = False) 
			# ~ plt.plot(B, inter_flux, label=r"{0} $\mu$m".format(wl.value))	
			# ~ inter_flux = np.zeros(len(B)) * u.Jy	
		
		# ~ plt.xlabel(r'Baseline B (m)')
		# ~ plt.ylabel('Correlated flux (Jy)')
		# ~ plt.title(r'{0} Correlated flux for multiple wavelengths, rv = 0.1, $T_{{am}}$ = {1}K'.format(disk, amor_temp.value))
		# ~ plt.legend()
		# ~ plt.savefig(folder + "Correlated_flux_multWL_rv0.1_T{0}.png".format(amor_temp.value))
		# ~ plt.show()
	
if __name__ == "__main__":
	main()
