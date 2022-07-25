# Calculates opacities, flux maps, spectra, correlated fluxes and all other related properties

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
      
 
      
def molecular_weight(solids):
	
	"""
	Calculates the molecular weight of the given solids in g
	
	Parameters:
	
	solids       : A 1D list of the names (strings) of the top 5 most abundant solids at all radii 
	
	Returns:
	
	molwt 		 : A dictionary of the molecular weights (in g) for each of the given solids where the keys are the solid names and the values are the molecular weights (float)
	"""
	
	molwt = {key: None for key in solids}
	
	for solid in molwt.keys():
		
		if solid == 'Olivine':
			f = Formula('MgFeSiO4')
			molwt[solid] = f.mass * u.g / u.mol / Na
			
		elif solid == 'Pyroxene':
			molwt[solid] = 109.8495 * u.g / u.mol / Na
		
		else:
			f = Formula(solid)
			molwt[solid] = f.mass * u.g / u.mol / Na
	
	return molwt



def surface_density(solids, molwt, top_abunds, nHtot, H, add_gap, R_arr, rgap, wgap, sgap):
	
	"""
	Calculates the surface density of the given solids. Note that the calculated "surface" densities are actually (g/cm^3), but we assume column height to be 1 cm, so the surface density can be g/cm^2. The surface densities can also be calculated by assuming a column height of H AU, but this does not give accurate results at the moment
	
	Parameters:
	
	solids        : A 1D list of the names (strings) of the top 5 most abundant solids at all radii
	molwt         : A dictionary of the molecular weights (in g) for each of the given solids where the keys are the solid names and the values are the molecular weights (float)
	top_abunds    : A dictionary of the abundances (in log10(nsolid/nH)) for each of the given solids where the keys are the solid names and the values are the array (shape: (NPOINT,)) of abundances (float)
	nHtot         : A 1D array (shape: (NPOINT,)) of total Hydrogen nuclei particle density (/cm^3) at each radius point (float)
	H 			  : 1 1D array (shape: (NPOINT,)) of pressure scale heights at each radius point in AU (float)
	
	Returns:
	
	surf_dens 	  : A dictionary of surface density values in g/cm^2 for each solid at every radius where the keys are the solid names and the values are the array (shape: (NPOINT,)) of surface densities (float)
	"""
	
	surf_dens = {key: None for key in solids}
	
	for solid in surf_dens.keys():
		
		if solid == 'Olivine':
			n_solid = nHtot * 10**top_abunds['Mg2SiO4']
			surf_dens[solid] = molwt[solid] * n_solid
			# surf_dens[solid] = H * u.cm * np.sqrt(2*np.pi) * surf_dens[solid] * np.exp(0.5)   # Assuming a column height of H AU and using the rho-Sigma formula for protodisks
			surf_dens[solid] = H.to(u.cm) * surf_dens[solid]
			
		elif solid == 'Pyroxene':
			n_solid = nHtot * 10**top_abunds['MgSiO3']
			surf_dens[solid] = molwt[solid] * n_solid
			# surf_dens[solid] = H * u.cm * np.sqrt(2*np.pi) * surf_dens[solid] * np.exp(0.5)   # Assuming a column height of H AU and using the rho-Sigma formula for protodisks
			surf_dens[solid] = H.to(u.cm) * surf_dens[solid]
		
		else:
			n_solid = nHtot * 10**top_abunds[solid]
			surf_dens[solid] = molwt[solid] * n_solid
			# surf_dens[solid] = H * u.cm * np.sqrt(2*np.pi) * surf_dens[solid] * np.exp(0.5)   # Assuming a column height of H AU and using the rho-Sigma formula for protodisks
			surf_dens[solid] = H.to(u.cm) * surf_dens[solid]
		
		# Adding a gap in the disk	
	if add_gap:
		
		# Find indices corresponding to gap radius and width
		print((rgap - wgap/2.0), (rgap + wgap/2.0))
		print(np.round(rgap - wgap/2.0, 1), np.round(rgap + wgap/2.0, 1))
		rgap_ind = np.where(np.round(R_arr, 1) == rgap)[0][0]
		wgap_ind1 = np.where(np.round(R_arr, 1) == np.round(rgap - wgap/2.0, 1))[0][0]
		wgap_ind2 = np.where(np.round(R_arr, 1) == np.round(rgap + wgap/2.0, 1))[0][0]
		
		for solid in surf_dens.keys():
			
			for i in range(wgap_ind1, wgap_ind2 + 1):
				
				surf_dens[solid][i] = sgap * surf_dens[solid][i]
				# ~ if i <= rgap_ind:
					
					# ~ m = (sgap * surf_dens[solid][rgap_ind] - surf_dens[solid][wgap_ind1])/(R_arr[rgap_ind] - R_arr[wgap_ind1])
					# ~ surf_dens[solid][i] = surf_dens[solid][wgap_ind1] + m * (R_arr[i] - R_arr[wgap_ind1])
				
				# ~ else:
					
					# ~ m = (surf_dens[solid][rgap_ind] - surf_dens[solid][wgap_ind2])/(R_arr[rgap_ind] - R_arr[wgap_ind2])
					# ~ surf_dens[solid][i] = surf_dens[solid][rgap_ind] + m * (R_arr[i] - R_arr[rgap_ind])
		
	return surf_dens


def r_to_rad(R_arr, dist_pc):
	
	"""
	Converts radius to arcseconds to radians. This is required since, for spectra, we integrate over the solid angle instead of the radius
	(R_arr / dist_pc) * (pc / AU) converts radii to arseconds. This value is divided by 3600 to obtain the angle in degrees. This is then multiplied by π/180 to represent the angle in radians.
	
	Parameters:
	
	R_arr        : 1D array of radii in AU obtained from the temperature array in GGchem output based on the power law model (float)
	
	Returns:
	
	rad_arr      : 1D array of angles in radians (float)
	"""
	
	rad_arr = R_arr/dist_pc * (1 * u.pc)/(1 * u.AU) / 3600.0 * np.pi / 180.0 * u.rad
	return rad_arr
	
	
	
def slice_lQ(lamda, Q, lmin, lmax, lsize):
	
	"""
	Slices the wavelength array to the given limits lmin, lmax and of the given size lsize and also obtains the coresponding Q/Kappa values. This is done because the given opacity curves have different types of wavelength sampling, so to ensure uniformity while performing array operations, all the opacities are to be taken between the same wavelengths (lmin, lmax) and be of the same size: lsize.
	
	Parameters:
	
	lamda        : 1D array of wavelengths for the solid in microns (size varies with the solid)
	Q  			 : 1D array of the absorption efficiencies of the solid (size varies with the solid)
	lmin         : Lower limit of wavelength in microns (float)
	lmax         : Upper limit of wavelength in microns (float)
	lsize        : Number of wavelength (and Q/kappa) points (int)
	
	Returns:
	
	lamda        : 1D array of wavelengths (shape: (lsize,)) of the solid from lmin to lmax in microns (float) 
	Q 			 : 1D array of absorption efficiencies (or kappas) (shape: (lsize,)) of the solid (if kappa, unit = (cm^2/g), otherwise, unitless) (float)
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
		old_indices = np.arange(0, len(lamda))
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
	kappa = 3 * Q / (4 * gs.to(u.cm) * dens)
	
	return mineral, rv, fmax, lamda, kappa



def Plancks(lamda, Tg):
	
	"""
	Calculates the spectral radiance (spectral emissive power per unit area, per unit solid angle, per unit frequency for particular radiation frequencies) of the solid using Planck's law for a range of radii (converted from temperatures based on the disk model) and a range of wavelengths in erg/(s Hz sr cm^2) 
	
	Parameters:
	
	lamda        : 1D array of wavelengths (shape: (lsize,)) of the solid in microns (float)
	Tg  		 : 1D array of gas temperatures (shape: (NPOINT,)) in K (float)
	
	Returns:
	
	I            : 2D array (shape: (lsize, NPOINT)) of spectral radiance of the solid in erg/(s Hz sr cm^2) i.e. CGS units (float)
	"""
		
	# Transposing 1D row array (shape: (lsize,)) into 2D column array (shape: (lsize, 1)) for matrix multiplication
	lamda = lamda[np.newaxis]
	lamda = lamda.T
	
	lamda_cm = lamda.to(u.cm)
	
	# Planck's function in terms of frequency
	v = c/lamda_cm
	bb = BlackBody(temperature = Tg, scale = 1.0)
	I = bb(v)
	
	return I



def tau_calc(sigma, kappa):
    
    """
    Calculates the optical depth tau for a solid (unitless)
    
    Parameters:
    
    sigma         : 1D array of the surface density (shape: (NPOINT,)) of the solid in g/cm^2 (float)
    kappa         : 1D array of the opacities (shape: (lsize,)) of the solid in cm^2/g (float)
    
    Returns:
    
    tau           : 2D array of the optical depths (shape: (NPOINT, lsize)) of the solid (float). The quantity is unitless.
    """

    # Transposing 1D row array (shape: (NPOINT,)) into 2D column array (shape: (NPOINT, 1)) for matrix multiplication
    sigma = sigma[np.newaxis]
    sigma = sigma.T

    tau = sigma * kappa
    
    return tau
    
    

def tau_calc_amorphous(sigma, sigma_am, kappa, kappa_am, Tg, amor_temp):
	
	"""
	Calculates the optical depths tau for a solid (unitless) by considering the amorphous opacities for temperatures below 500K
	
	Parameters:
	
	sigma         : 1D array of the surface density (shape: (NPOINT,)) of the solid in g/cm^2 (float)
	sigma_am 	  : 1D array of the amorphous surface density (shape: (NPOINT,)) of the solid in g/cm^2 (float)
	kappa         : 1D array of the opacities (shape: (lsize,)) of the solid in cm^2/g (float)
	kappa_am      : 1D array of the amorphous opacities (shape: (lsize,)) of the solid in cm^2/g (float)
	Tg  		  : 1D array of gas temperatures (shape: (NPOINT,)) in K (float)
	amor_temp 	  : Temperature in K at and above which the grains are considered to be crystalline and below which the grains are amorphous (float)
	
	Returns:
	
	tau           : 2D array of the optical depths (shape: (NPOINT, lsize)) of the solid (float). The quantity is unitless.
	"""
	
	tau = np.zeros((len(sigma), len(kappa)))
	l_ind = np.where(Tg < (amor_temp))[0][0]     # Taking the first index where the temperature goes below 500K. Due to the power law, we can assume every index after this also has T < 500K
	
	# Transposing 1D row array (shape: (NPOINT,)) into 2D column array (shape: (NPOINT, 1)) for matrix multiplication
	sigma = sigma[np.newaxis]
	sigma = sigma.T
	sigma_am = sigma_am[np.newaxis]
	sigma_am = sigma_am.T
	
	# Considering the crystalline kappas for temperatures at amor_temp and above, and amorphous kappas for temperatures below amor_temp
	tau[:l_ind, :] = sigma[:l_ind, :] * kappa
	tau[l_ind:, :] = sigma_am[l_ind:, :] * kappa_am
	
	return tau
    
    

def flux_map(tau, I):

    """
    Calculates the flux map for a solid in erg/(s Hz sr cm^2)
    
    Parameters:
    
    tau          : 2D array of the optical depths (shape: (NPOINT, lsize)) of the solid (float). The quantity is unitless.
    I 			 : 2D array (shape: (lsize, NPOINT)) of spectral radiance of the solid in erg/(s Hz sr cm^2) i.e. CGS units (float)
    
    Returns:
    
	F_map 		 : 2D array (shape: (NPOINT, lsize)) of the flux map for the solid in erg/(s Hz sr cm^2) i.e. CGS units (float)
    """
	
    I = I.T                                             # Transposing I from (lambda, r) to (r, lambda) for element-wise matrix multiplication
    F_map = (1 - np.exp(-tau)) * I                      # Calculating the flux map
    
    return F_map

    
    
def calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc):
	
	"""
	Calculates the integrated flux for the range of wavelength by integrating over the radius limits given (Rmin, Rmax)
	
	Parameters:
	
	F_map 		 : 2D array (shape: (NPOINT, lsize)) of the flux map for a single solid in erg/(s Hz sr cm^2) i.e. CGS units (float)
	R_arr        : 1D array of radii in AU obtained from the temperature array in GGchem output based on the power law model (float)
	Rmin 		 : Minimum radius limit over which the flux map is integrated (float)
	Rmax 		 : Maximum radius limit over which the flux map is integrated (float)
	
	Returns: 
	
	summ   		 : 1D array (shape: (lsize,)) of the integrated flux in Jy (float)
	"""
	
	# Finding the indices of Rmin and Rmax by taking the first instance of where they are in the rounded R_arr 
	R_rounded = np.round(R_arr, 3)	
	Rmin_id = np.where(R_rounded == Rmin)[0][0]
	Rmax_id = np.where(R_rounded == Rmax)[0][0] 
	
	# Obtaining only that part of the arrays that lie between Rmin and Rmax (both included)
	rad_arr = r_to_rad(R_arr[Rmin_id: Rmax_id+1], dist_pc) 						# Converting radius (AU) to arcseconds to radians
	F_map_req = F_map[Rmin_id: Rmax_id+1, :]
	
	# Transposing the matrices for element-wise matrix multiplication
	rad_arr = rad_arr[np.newaxis]
	rad_arr = rad_arr.T 
	
	# Calculating the integrated flux using the Trapezoidal numerical integration
	summ = np.trapz(2 * np.pi * rad_arr * F_map_req, x = rad_arr, axis = 0)
	summ = summ.to(u.Jy, equivalencies = u.dimensionless_angles())

	return summ * 10**4



def hankel_transform(F_map, R_arr, lamda, wl, B, dist_pc, wl_array):
	
	"""
	Calculates the absolute correlated flux density (in Jy) using the combined flux map of all the solids. This involves finding the Fourier transform in cylindrical coordinates using the Bessel function of the first kind zeroth order (j0). It is assumed that the correlated fluxes (shape: (lsize,)) are to be plotted against the wavelengths. If, however, it is to be plotted against the baselines, then the correlated flux density (a single float value) corresponding to a specific wavelength (wl) is extracted and returned from the array (shape: (lsize,)). This can be specified with wl_array.
	
	Parameters: 
	
	F_map 		  : 2D array (shape: (NPOINT, lsize)) of the combined flux map for all solids in erg/(s Hz sr cm^2) i.e. CGS units (float)
	R_arr         : 1D array of radii in AU obtained from the temperature array in GGchem output based on the power law model (float)
	lamda         : 1D array of wavelengths (shape: (lsize,)) of the solid in microns (float)
	wl 			  : A single wavelength (in micron) where the correlated flux is required if the fluxes are plot against baseline (float). It is used only if wl_array = False
	B 			  : The interferometric baseline value in m (float)
	wl_array 	  : If True, the returned correlated fluxes is an array of shape (lsize,) to be plotted against the wavelengths. If False, the correlated flux for a single wavelength value is returned to be plot against the baselines (Boolean)
	
	Returns:
	
	inter_flux_abs: The absolute correlated flux density in Jy. If wl_array = True, a float array of shape (lsize,) is returned, else a single float value is returned
	"""
	
	rad_arr = r_to_rad(R_arr, dist_pc) 						 											# Converting radius (AU) to arcseconds to radians	
	rad_arr = rad_arr[np.newaxis]
	rad_arr = rad_arr.T 																		# Shape: (NPOINT, 1)
	rad_arr = rad_arr.to('', equivalencies=u.dimensionless_angles()) 
	
	# Finding the Bessel function of the first kind, zeroth order
	bessel = j0(2.0 * np.pi * rad_arr * B / (lamda.to(u.m)))           							# Shape: (NPOINT, lsize)
	
	# Calculating the absolute interferometric flux
	inter_flux = 2.0 * np.pi * np.trapz(rad_arr * bessel * F_map, x = rad_arr, axis = 0)   		# Shape: (lsize,)
	inter_flux_abs = np.abs(inter_flux) * u.rad**2
	inter_flux_abs = inter_flux_abs.to(u.Jy, equivalencies = u.dimensionless_angles())
	
	if wl_array: 
		
		return inter_flux_abs * 10**4
		
	else:
		
		# Finding the index where the array of wavelengths matches the given wavelength wl
		lamda_rounded = np.round(lamda, 1)	
		wl_id = np.where(lamda_rounded == wl)[0][0]
		
		inter_flux_abs = inter_flux_abs[wl_id]                                             		# Single float value
		
		return inter_flux_abs * 10**4
