import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.constants import astropyconst20 as const
from astropy.modeling.models import BlackBody
from fancy_name import latex_name
from molmass import Formula
from diskprop import midplaneT_profile
from scipy.interpolate import UnivariateSpline
from scipy.special import j0
plt.rcParams['axes.titlesize'] = 10
plt.rcParams["figure.figsize"] = (8, 6)


# Some constants in CGS
Na = const.N_A.cgs                    		# Avogadro's number in /mol
h = const.h.cgs                   			# Planck's constant in cm^2 g s-1
c = const.c.cgs              				# Speed of light in cm/s              
k = const.k_B.cgs                   		# Boltzmann constant in cm^2 g s^-2 K^-1
AU = const.au.cgs                			# 1 astronomical unit in cm
dist_pc = 100 * u.pc  			            # Assuming a distance to the Sun-like star in parsec
      
 
      
def molecular_weight(solids):
    
    """
    Calculates the molecular weight of the given solids in g/mol

    Parameters:

    solids       : A 1D list of the names (strings) of the top 5 most abundant solids at all radii 

    Returns a dictionary of the molecular weights (in g/mol) for each of the given solids where the keys are the solid names and the values are the molecular weights (float)
    """
    
    molwt = {key: None for key in solids}

    for solid in molwt.keys():
        f = Formula(solid)
        molwt[solid] = f.mass * u.g / u.mol / Na
    
    return molwt



def surface_density(solids, molwt, top_abunds, nHtot, H):
	
	"""
	Calculates the surface density of the given solids. Note that the calculated "surface" densities are actually (g/cm^3), but we assume column height to be 1 cm, so the surface density can be g/cm^2
	
	Parameters:
	
	solids        : A 1D list of the names (strings) of the top 5 most abundant solids at all radii
	molwt         : A dictionary of the molecular weights (in g/mol) for each of the given solids where the keys are the solid names and the values are the molecular weights (float)
	top_abunds    : A dictionary of the abundances (in log10(nsolid/nH)) for each of the given solids where the keys are the solid names and the values are the array (shape: (NPOINT,)) of abundances (float)
	nHtot         : A 1D array (shape: (NPOINT,)) of total Hydrogen nuclei particle density (/cm^3) at each radius point (float)
	
	Returns a dictionary of surface density values in g/cm^2 for each solid at every radius where the keys are the solid names and the values are the array (shape: (NPOINT,)) of surface densities (float)
	"""
	
	# Transposing 1D row array (shape: (1, NPOINT)) into 2D column array (shape: (NPOINT, 1)) for matrix multiplication
	# ~ nHtot = nHtot[np.newaxis]
	# ~ nHtot = nHtot.T
	
	H = 1 * u.cm
	surf_dens = {key: None for key in solids}
	
	for solid in surf_dens.keys():
		n_solid = nHtot * 10**top_abunds[solid]
		surf_dens[solid] = molwt[solid] * n_solid
		# surf_dens[solid] = H * np.sqrt(2*np.pi) * surf_dens[solid] * np.exp(0.5)    # Assuming a column height of H AU
		# surf_dens[solid] = H * u.cm * np.sqrt(2*np.pi) * surf_dens[solid] * np.exp(0.5)
		surf_dens[solid] = H.to(u.cm) * surf_dens[solid]
		
	return surf_dens


def r_to_rad(R_arr):
	
	"""
	Converting radius to arcseconds to radians
	"""
	
	rad_arr = R_arr/dist_pc * (1 * u.pc)/(1 * u.AU)/ 3600.0 * np.pi / 180.0 * u.rad
	return rad_arr
	
	
	
def slice_lQ(lamda, Q, lmin, lmax, lsize):
	
	"""
	Slices the arrays to the given limits lmin, lmax and of the given size lsize
	"""
	
	# Taking only the indices where the wavelength is between lmin to lmax microns
	indices = np.where(np.logical_and(lamda >= lmin, lamda <= lmax))
	lamda = lamda[indices]
	Q = Q[indices]
	
	# Ensuring that there are an equal number of wavelength i.e. lsize points across all solids
	if len(lamda) > lsize:
		
		# If there are more than lsize points of lamda, roughly equally spaced lsize points are sampled from the array and lamda is reassigned accordingly
		idx = np.round(np.linspace(0, len(lamda) - 1, lsize)).astype(int) 
		lamda = lamda[idx]
		Q = Q[idx]
	
	elif len(lamda) < lsize:
		
		# If there are less than lsize points of lamda, lsize points are created through interpolation and lamda is reassigned accordingly
		old_indices = np.arange(0,len(lamda))
		new_indices = np.linspace(0, len(lamda)-1, lsize)
		
		spl1 = UnivariateSpline(old_indices, lamda, k=3, s=0)
		lamda = spl1(new_indices)
		
		spl2 = UnivariateSpline(old_indices, Q, k=3, s=0)
		Q = spl2(new_indices)
		
	return lamda, Q
	
		
		
def get_l_and_k(opfile, dens, gs, lmin, lmax, lsize):
	
	"""
	Obtains the wavelength and Q_abs values for a condensate from the given opfile, and calculates the corresponding kappa values
	
	Parameters:
	
	opfile       : Path to the file that contains wavelength and Qabs values for a solid (string)
	dens         : Material density of the solid in g/cm^3 (float)
	gs           : Grain radius in cm (float)
	lmin         : Lower limit of wavelength in microns (float)
	lmax         : Upper limit of wavelength in microns (float)
	lsize        : Number of wavelength (and kappa) points
	
	Returns:
	
	mineral      : The name of the solid as extracted from the opfile (string)
	rv           : Grain radius used in the opfile in microns (float)
	fmax         : Maximum emptiness fraction value used for the solid opacities according to the DHS theory (float)
	lamda        : 1D array of wavelengths for the solid from lmin to lmax (shape: (lsize,)) (float)
	kappa        : 1D array of opacities of the solid in cm^2/g (shape: (lsize,)) (float)
	
	Note that the returned variables are then stored as values for the corresponding solids in the dictionaries in main() 
	"""
	
	# Obtaining the opacity file metadata such as the mineral name, grain size value and fmax 
	filename = opfile.split('/')[1]
	values = filename.split('_')
	mineral = values[1]
	
	for i in range(len(values)):
	
		if values[i].startswith('rv'):
			rv = values[i][2:5] 
		elif values[i].startswith('fmax'):
			fmax = values[i][4:7]
	
	# Since the opfile for H2O is created using optool which has a different output file, its wavelength and opacity are extracted differently 
	if mineral == 'H2O':
		
		Q_curve = pd.read_csv(opfile, delimiter='\s+', skiprows = 1, names = ['wavelength', 'Kabs', 'Ksca', 'g_asymm'])
		lamda = u.Quantity(Q_curve['wavelength'].to_numpy(), u.micron)
		kappa = Q_curve['Kabs'].to_numpy() * u.cm**2 / u.g
		lamda, kappa = slice_lQ(lamda, kappa, lmin, lmax, lsize)
		
		return mineral, rv, fmax, lamda, kappa
		
	Q_curve = pd.read_csv(opfile, delimiter='\s+', skiprows=1, names=['wavelength','Q'])
	
	lamda = Q_curve['wavelength'].to_numpy() * u.micron
	Q = Q_curve['Q'].to_numpy()
	
	lamda, Q = slice_lQ(lamda, Q, lmin, lmax, lsize)
	
	# Obtaining k_abs from Q_abs
	kappa = 3 * Q / (4 * gs * dens)
	
	return mineral, rv, fmax, lamda, kappa
	
	
	
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
	
	plt.plot(lamda, kappa)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel(r'$\kappa_{abs}$ ($cm^2/g$)')
	plt.title(r"Q-curve for {0}, r = {1}, $f_{{max}}$ = {2}".format(latex_name(mineral), rv, fmax))
	plt.savefig(folder + "Qcurve_{0}_r{1}_f{2}.png".format(mineral, rv, fmax), bbox_inches = 'tight')
	plt.show()



def Plancks(T0, R_arr, R_in, lamda, q, solid, folder):
	
	"""
	Calculates the spectral radiance (spectral emissive power per unit area, per unit solid angle, per unit frequency for particular radiation frequencies) of the solid using Planck's law for a range of radii (converted from temperatures based on the disk model) and a range of wavelengths in erg/(s Hz sr cm^2) 
	
	Parameters:
	
	T0           : Dust sublimation temperature in K (float)
	R_arr        : 1D array of radii in AU obtained from the temperature array in GGchem output based on the power law model (float)
	R_in         : Inner radius of the planet forming disk inside which dust grains cannot condense out from the gas in AU (float)
	lamda        : 1D array of wavelengths (shape: (lsize,)) of the solid in microns (float)
	
	Returns:
	
	I            : 2D array (shape: (lsize, NPOINT)) of spectral radiance of the solid in erg/(s Hz sr cm^2) i.e. CGS units 
	"""
		
	# Transposing 1D row array (shape: (lsize,)) into 2D column array (shape: (lsize, 1)) for matrix multiplication
	lamda = lamda[np.newaxis]
	lamda = lamda.T
	
	lamda_cm = lamda.to(u.cm)
	
	# Planck's function in terms of frequency
	v = c/lamda_cm
	T = midplaneT_profile(R_in, T0, R_arr, q)
	bb = BlackBody(temperature=T, scale=1.0)
	I = bb(v)
	
	plt.plot(lamda, I)
	plt.title("{0} Spectral radiance vs Wavelength".format(solid))
	plt.ylabel("Spectral radiance ({0})".format(I.unit))
	plt.xlabel(r"$\lambda$ ($\mu$m)")
	plt.savefig(folder + "Bv_vs_wl_{0}_rv0.1.png".format(solid))
	plt.show()
	
	return I



def tau_calc(sigma, kappa, solid, folder):
    
    """
    Calculates the optical depth tau for a solid (unitless)
    
    Parameters:
    
    sigma         : 1D array of the surface density (shape: (NPOINT,)) of the solid in g/cm^2 (float)
    kappa         : 1D array of the opacities (shape: (lsize,)) of the solid in cm^2/g (float)
    
    Returns:
    
    tau           : 2D array of the optical depths (shape: NPOINT, lsize)) of the solid (float). The quantity is unitless.
    """

    # Transposing 1D row array (shape: (NPOINT,)) into 2D column array (shape: (NPOINT, 1)) for matrix multiplication
    sigma = sigma[np.newaxis]
    sigma = sigma.T

    # tau = sigma * kappa * 10**-8     # Assuming the dust scale height is 8 times lesser than the gas scale height. Why? Because I'm desperate.
    tau = sigma * kappa
    
    # Plotting the optical depths
    plt.imshow(tau)
    plt.colorbar()
    plt.title("Optical depth - {0}".format(solid))
    plt.savefig(folder + "OptDepth_{0}.png".format(solid))
    plt.show()
    
    return tau
    
    

def flux_map(tau, I):

    """
    Calculates the flux map for a solid (r, lambda)
    """
	
    I = I.T                                             # Transposing I from (lambda, r) to (r, lambda) for element-wise matrix multiplication
    F_map = (1 - np.exp(-tau)) * I                      # Calculating the flux map
    
    return F_map
    


def plot_fluxmap(solid_name, rv, fmax, F_map, lamda, R_arr, folder):
	
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
	fig.colorbar(img)
	plt.savefig(folder + "Fluxmap_{0}_r{1}_fmax{2}.png".format(solid_name, rv, fmax), bbox_inches = 'tight')
	plt.show()

    
    
def calculate_spectra(tau, F_map, I, R_arr, Rmin, Rmax, last_id):
	
	"""
	Calculates the integrated flux for the range of wavelength by integrating over radius
	
	Summ shape (lamda, 1)
	"""
	
	# Finding the indices of Rmin and Rmax by taking the first instance of where they are in the rounded R_arr 
	R_rounded = np.round(R_arr, 3)	
	Rmin_id = np.where(R_rounded == Rmin)[0][0]
	# Rmax_id = np.where(R_rounded == Rmax)[0][0]
	Rmax_id = last_id 
	
	# Transposing the matrices for element-wise matrix multiplication
	rad_arr = r_to_rad(R_arr) 						# Converting radius (AU) to arcseconds to radians
	rad_arr = rad_arr[np.newaxis]
	rad_arr = rad_arr.T 
	I = I.T
	
	summ = np.trapz(2 * np.pi * rad_arr * tau * I, x = rad_arr, axis = 0)
	
	# delr = (rad_arr[Rmin_id+1: Rmax_id+1, :] - rad_arr[Rmin_id: Rmax_id, :])	
	# f1 = rad_arr[Rmin_id: Rmax_id, :] * tau[Rmin_id: Rmax_id, :] * I[Rmin_id: Rmax_id, :] * 2 * np.pi * delr * 0.5
	# f2 = rad_arr[Rmin_id+1: Rmax_id+1, :] * tau[Rmin_id+1: Rmax_id+1, :] * I[Rmin_id+1: Rmax_id+1, :] * 2 * np.pi * delr * 0.5
	
	# Using same formula as Hankel transform one
	# f1 = rad_arr[Rmin_id: Rmax_id, :] * F_map[Rmin_id: Rmax_id, :] * 2 * np.pi * delr * 0.5
	# f2 = rad_arr[Rmin_id+1: Rmax_id+1, :] * F_map[Rmin_id+1: Rmax_id+1, :]  * 2 * np.pi * delr * 0.5
	# f1, f2 shape: (NPOINT - 1, lsize)
	
	# temp = (f1 + f2)
	# summ = temp.sum(axis=0)       # summ shape: (lsize, 1)
	summ = summ.to(u.Jy, equivalencies = u.dimensionless_angles())

	# Inefficient for-loop method
	# ~ summ = 0
	# ~ for r1 in range(Rmin_id, Rmax_id-1):
		# ~ for r2 in range(r1+1, Rmax_id):
	
			# ~ # Numerical integration using the trapezoidal rule
			# ~ delr = R_arr[r2] - R_arr[r1]
			# ~ fr1 = f(tau[r1, :], I[:, r1], R_arr[r1], delr)
			# ~ fr2 = f(tau[r2, :], I[:, r2], R_arr[r2], delr)
			# ~ summ += (fr1 + fr2)
	
	return summ
	
	
	
def plot_spectra(lamda, summ, solid_name, rv, fmax, Rmin, Rmax, folder):
	
	fig = plt.figure()
	plt.plot(lamda, summ)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux (Jy)')
	plt.title(r'Spectrum {0} r={1} $\mu$m $f_{{max}}$={2} R={3}-{4} AU'.format(latex_name(solid_name), rv, fmax, Rmin.value, Rmax.value))
	plt.savefig(folder + "Spectrum_{0}_r{1}_f{2}_R{3}-{4}.png".format(solid_name, rv, fmax, Rmin.value, Rmax.value))
	plt.show()



def hankel_transform(F_map, R_arr, lamda, wl, B, wl_array):
	
	"""
	Performs Fourier transform in cylindrical coordinates to get the interferometric flux map (i.e. correlated flux density)
	I (lsize, NPOINT)
	"""
	
	rad_arr = r_to_rad(R_arr) 						# Converting radius (AU) to arcseconds to radians	
	rad_arr = rad_arr[np.newaxis]
	rad_arr = rad_arr.T 							# Shape: (NPOINT, 1)
	rad_arr = rad_arr.to('', equivalencies=u.dimensionless_angles()) 
	
	# Considering F_map only at the given wavelength because I am out of ideas
	# lamda_rounded = np.round(lamda, 1)	
	# wl_id = np.where(lamda_rounded == wl)[0][0]
	# F_map = F_map[:, wl_id]
	
	# One more idea, summing up the wavelength axis
	# F_map = np.sum(F_map, axis=1)
	
	bessel = j0(2.0 * np.pi * rad_arr * B / (lamda.to(u.m)))           							# Shape: (NPOINT, lsize)
	inter_flux = 2.0 * np.pi * np.trapz(rad_arr * bessel * F_map, x = rad_arr, axis = 0)   		# Shape: (lsize,)
	inter_flux_abs = np.abs(inter_flux) * u.rad**2
	inter_flux_abs = inter_flux_abs.to(u.Jy, equivalencies = u.dimensionless_angles())
	
	if wl_array: 
		
		return inter_flux_abs
		
	else:
		
		# Finding the index where the array of wavelengths matches the given wavelength
		lamda_rounded = np.round(lamda, 1)	
		wl_id = np.where(lamda_rounded == wl)[0][0]
		
		inter_flux_abs = inter_flux_abs[wl_id]                                             # Single float value
		
		return inter_flux_abs
