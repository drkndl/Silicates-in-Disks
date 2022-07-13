import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.constants import astropyconst20 as const
from fancy_name import latex_name


def Qcurve_plotter(lamda, kappa, mineral, rv, fmax, folder):
	
	"""
	Plots the wavelength (micron) vs opacity (cm^2/g) curve for a solid
	
	Parameters:
	
	lamda        : 1D array of wavelengths (shape: (lsize,)) in microns (float)
	kappa        : 1D array of opacities (shape: (lsize,)) in cm^2/g (float)
	mineral      : Name of the solid whose Qcurve is plotted, used for plot name (string)
	rv           : Grain radius of the solid in microns, used for plot name (float)
	fmax         : Maximum emptiness fraction of the solid grains according to the DHS theory (float)
	"""
	
	# Defining the names of the mineral for plot formatting. If the materials are amorphous (i.e. olivine or pyroxene), latex_name will not work and therefore their names in the plots are specifically defined
	if mineral == 'MgOlivine':
		fancy = 'MgOlivine'
	elif mineral == 'MgPyroxene':
		fancy = 'MgPyroxene'
	else:
		fancy = latex_name(mineral)		
		
	plt.plot(lamda, kappa)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel(r'$\kappa_{abs}$ ($cm^2/g$)')
	plt.title(r"Q-curve for {0}, r = {1}, $f_{{max}}$ = {2}".format(fancy, rv, fmax))
	plt.savefig(folder + "Qcurve_{0}_r{1}_f{2}.png".format(mineral, rv, fmax), bbox_inches = 'tight')
	plt.show()
	
	
	
def plot_Bv(lamda, I, solid, folder):
	
	"""
	Plots the spectral density of the solid against the wavelengths 
	
	Parameters:
	
	lamda        : 1D array of wavelengths (shape: (lsize,)) in microns (float)
	I            : 2D array (shape: (lsize, NPOINT)) of spectral radiance of the solid in erg/(s Hz sr cm^2) i.e. CGS units (float)
	solid 		 : Name of the solid for which the spectral density is plotted (string)
	folder       : Path where the output plots are saved (string)
	"""
	
	plt.plot(lamda, I[:, 0])
	plt.title("{0} Spectral radiance vs Wavelength".format(latex_name(solid)))
	plt.ylabel("Spectral radiance ({0})".format(I.unit))
	plt.xlabel(r"$\lambda$ ($\mu$m)")
	plt.savefig(folder + "Bv_vs_wl_{0}_rv0.1.png".format(solid))
	plt.show()
	
	
	
def plot_tau(tau, solid, folder):
	
	"""
	Plots the optical depth values for each solid
	
	Parameters:
	
	tau 		 : 2D array of the optical depths (shape: (NPOINT, lsize)) of the solid (float). The quantity is unitless.
	solid 		 : Name of the solid for which the spectral density is plotted (string)
	folder       : Path where the output plots are saved (string)
	"""
	
	plt.imshow(tau)
	plt.colorbar()
	plt.title("Optical depth - {0}".format(latex_name(solid)))
	plt.savefig(folder + "OptDepth_{0}.png".format(solid))
	plt.show()
	
	
	
def plot_fluxmap(solid_name, rv, fmax, F_map, lamda, R_arr, folder):
	
	"""
	Plots the flux map for a solid as a colormap
	
	Parameters:
	
	solid_name 	  : Name of the solid for which the spectral density is plotted (string)
	rv            : Grain radius of the solid in microns, used for plot name (float)
	fmax          : Maximum emptiness fraction of the solid grains according to the DHS theory (float)
	F_map         : 2D array (shape: (NPOINT, lsize)) of the flux map for the solid in erg/(s Hz sr cm^2) i.e. CGS units (float)
	lamda         : 1D array of wavelengths (shape: (lsize,)) of the solid in microns (float)
	R_arr         : 1D array of radii in AU obtained from the temperature array in GGchem output based on the power law model (float)
	folder        : Path where the output plots are saved (string)
	"""
	
	fig, ax = plt.subplots(1,1)
	img = ax.imshow(F_map, cmap='plasma', interpolation='none')
	
	# Axes formatting    
	x_axis_locations = np.linspace(0, len(lamda)-1, 8).astype(int)
	ax.set_xticks(x_axis_locations)
	x_axis_labels = np.round(lamda[x_axis_locations], 1)
	ax.set_xticklabels(x_axis_labels.value)
	ax.set_xlabel(r'$\lambda$ ($\mu$m)')
	
	y_axis_locations = np.linspace(0, len(R_arr)-1, 8).astype(int)
	ax.set_yticks(y_axis_locations)
	y_axis_labels = np.round(R_arr[y_axis_locations], 3)
	ax.set_yticklabels(y_axis_labels.value)
	ax.set_ylabel('R (AU)')
	
	ax.set_title(r"Flux Map for {0}, r = {1}, $f_{{max}}$ = {2}".format(latex_name(solid_name), rv, fmax))
	fig.colorbar(img, label=r'{0}'.format(F_map.unit))
	plt.savefig(folder + "Fluxmap_{0}_r{1}_fmax{2}.png".format(solid_name, rv, fmax), bbox_inches = 'tight')
	plt.show()
	
	
	
def plot_spectra(lamda, summ, solid_name, rv, fmax, Rmin, Rmax, folder):
	
	"""
	Plots the spectrum for a solid based on the integrated flux calculated in calculate_spectra()
	
	Parameters:
	
	lamda         : 1D array of wavelengths (shape: (lsize,)) of the solid in microns (float)
	summ   		  : 1D array (shape: (lsize,)) of the integrated flux in Jy (float)
	solid_name 	  : Name of the solid for which the spectral density is plotted (string)
	rv            : Grain radius of the solid in microns, used for plot name (float)
	fmax          : Maximum emptiness fraction of the solid grains according to the DHS theory (float)
	Rmin 		  : Minimum radius limit over which the flux map is integrated (float)
	Rmax 		  : Maximum radius limit over which the flux map is integrated (float)
	folder        : Path where the output plots are saved (string)
	"""
	
	fig = plt.figure()
	plt.plot(lamda, summ)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title(r'Spectrum {0} r={1} $\mu$m $f_{{max}}$={2} R={3}-{4} AU'.format(latex_name(solid_name), rv, fmax, Rmin.value, Rmax.value))
	plt.savefig(folder + "Spectrum_{0}_r{1}_f{2}_R{3}-{4}.png".format(solid_name, rv, fmax, Rmin.value, Rmax.value))
	plt.show()
