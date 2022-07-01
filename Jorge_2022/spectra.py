import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from molmass import Formula
from pyvalem.formula import Formula
from jorge_diskprop import inner_radius, r_from_T
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from scipy.interpolate import UnivariateSpline
from scipy.special import j0
plt.rcParams['axes.titlesize'] = 10



# Some constants in CGS
Na = 6.022E23                    # Avogadro's number in /mol
h = 6.6261e-27                   # Planck's constant in cm^2 g s-1
c = 2.99792458e10                # Speed of light in cm/s                
k = 1.3807e-16                   # Boltzmann constant in cm^2 g s^-2 K^-1
bar = 1.E+6                      # 1 bar in dyn/cm^2
AU = 1.496E13                    # 1 astronomical unit in cm
Jy = 1.E-23                      # Converting flux from CGS units to Jy



def latex_name(solid):
	
	"""
	Returns a string of latex friendly name of a compound for better plot formatting
	
	Parameters:
	
	solids         : The formula of the condensate (string)
	
	Example:
	
	>>> latex_name('CaMgSi2O6')
	>>> $CaMgSi_2O_6$
	"""
	
	f = Formula(solid)
	fancy = "$" + f.latex + "$"
	raw_solid = r"{}".format(fancy)
	
	return raw_solid
      
 
      
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
        molwt[solid] = f.mass / Na
    
    return molwt



def surface_density(solids, molwt, top_abunds, nHtot):
    
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

    surf_dens = {key: None for key in solids}

    for solid in surf_dens.keys():
        n_solid = nHtot * 10**top_abunds[solid]
        surf_dens[solid] = molwt[solid] * n_solid
    
    return surf_dens



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
	
	# Since the opfile for H2O is created using optool which has a different output file, its wavelength and opacity is extracted differently 
	if mineral == 'H2O':
		
		Q_curve = pd.read_csv(opfile, delimiter = '\s+', skiprows = 1, names = ['wavelength', 'Kabs', 'Ksca', 'g_asymm'])
		lamda = Q_curve['wavelength'].to_numpy()
		kappa = Q_curve['Kabs'].to_numpy()
		lamda, kappa = slice_lQ(lamda, kappa, lmin, lmax, lsize)
		
		return mineral, rv, fmax, lamda, kappa
		
	Q_curve = pd.read_csv(opfile, delimiter='\s+', skiprows=1, names=['wavelength','Q'])
	
	lamda = Q_curve['wavelength'].to_numpy()
	Q = Q_curve['Q'].to_numpy()
	
	lamda, Q = slice_lQ(lamda, Q, lmin, lmax, lsize)
	
	# Obtaining k_abs from Q_abs
	kappa = 3 * Q / (4 * gs * dens)
	
	return mineral, rv, fmax, lamda, kappa
	
	
	
def Qcurve_plotter(lamda, kappa, mineral, rv, fmax):
	
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
	plt.savefig("Qcurves_trimmed/Qcurve_{0}_r{1}_f{2}.png".format(mineral, rv, fmax), bbox_inches = 'tight')
	plt.show()



def Plancks(T0, R_arr, R_in, lamda):
	
	"""
	Calculates the spectral radiance (spectral emissive power per unit area, per unit solid angle, per unit frequency for particular radiation frequencies) of the solid using Planck's law for a range of radii (converted from temperatures based on the disk model) and a range of wavelengths in erg/(s cm^3 sr)
	
	Parameters:
	
	T0           : Dust sublimation temperature in K (float)
	R_arr        : 1D array of radii in AU obtained from the temperature array in GGchem output based on the power law model (float)
	R_in         : Inner radius of the planet forming disk inside which dust grains cannot condense out from the gas in AU (float)
	lamda        : 1D array of wavelengths (shape: (lsize,)) of the solid in microns (float)
	
	Returns:
	
	I            : 2D array (shape: (lsize, NPOINT)) of spectral radiance of the solid in erg/(s cm^3 sr) i.e. CGS units 
	"""
		
	# Transposing 1D row array (shape: (lsize,)) into 2D column array (shape: (lsize, 1)) for matrix multiplication
	lamda = lamda[np.newaxis]
	lamda = lamda.T
	
	lamda_cm = lamda * 10**-4
	
	# Planck's function in terms of frequency
	v = c/lamda_cm
	T = T0 * (R_arr/R_in)**(-3/4)    
	I = 2.0*h*v**3.0 / (c**2 * (np.exp(h*v/(k*T)) - 1.0))
	
	return I



def tau_calc(sigma, kappa):
    
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

    tau = sigma * kappa
    
    return tau
    
    

def flux_map(tau, I):

    """
    Calculates the flux map for a solid (r, lambda)
    """
	
    I = I.T                                             # Transposing I from (lambda, r) to (r, lambda) for element-wise matrix multiplication
    F_map = (1 - np.exp(-tau)) * I                      # Calculating the flux map
    
    return F_map
    


def plot_fluxmap(solid_name, rv, fmax, F_map, lamda, R_arr):

    fig, ax = plt.subplots(1,1)
    img = ax.imshow(F_map, cmap='plasma', interpolation='none')

    # Axes formatting    
    x_axis_locations = np.linspace(0, len(lamda)-1, 8).astype(int)
    ax.set_xticks(x_axis_locations)
    x_axis_labels = np.round(lamda[x_axis_locations], 1)
    ax.set_xticklabels(x_axis_labels)
    ax.set_xlabel(r'$\lambda$ ($\mu$m)')

    y_axis_locations = np.linspace(0, len(R_arr)-1, 8).astype(int)
    ax.set_yticks(y_axis_locations)
    y_axis_labels = np.round(R_arr[y_axis_locations], 3)
    ax.set_yticklabels(y_axis_labels)
    ax.set_ylabel('R (AU)')

    ax.set_title(r"Flux Map for {0}, r = {1}, $f_{{max}}$ = {2}".format(latex_name(solid_name), rv, fmax))
    fig.colorbar(img)
    plt.savefig("Flux_Maps_trimmed/{0}_fluxmap_r{1}_fmax{2}.png".format(solid_name, rv, fmax), bbox_inches = 'tight')
    plt.show()

    
    
def calculate_spectra(tau, I, R_arr, Rmin, Rmax):
	
	"""
	Plots the integrated flux vs wavelength
	
	Summ shape (lamda, 1)
	"""
	
	# Finding the indices of Rmin and Rmax by taking the first instance of where they are in the rounded R_arr 
	R_rounded = np.round(R_arr, 2)	
	Rmin_id = np.where(R_rounded == Rmin)[0][0]
	Rmax_id = np.where(R_rounded == Rmax)[0][0]
	
	# Transposing the matrices for element-wise matrix multiplication
	R_arr = R_arr * AU
	R_arr = R_arr[np.newaxis]
	R_arr = R_arr.T 
	I = I.T
	
	delr = (R_arr[Rmin_id+1: Rmax_id+1, :] - R_arr[Rmin_id: Rmax_id, :])	
	f1 = R_arr[Rmin_id: Rmax_id, :] * tau[Rmin_id: Rmax_id, :] * I[Rmin_id: Rmax_id, :] * 2 * np.pi * delr * 0.5
	f2 = R_arr[Rmin_id+1: Rmax_id+1, :] * tau[Rmin_id+1: Rmax_id+1, :] * I[Rmin_id+1: Rmax_id+1, :] * 2 * np.pi * delr * 0.5
	# f1, f2 shape: (NPOINT - 1, lsize)
	
	temp = (f1 + f2)
	summ = temp.sum(axis=0)       # summ shape: (lsize, 1)

	# Inefficient for loop method
	# ~ summ = 0
	# ~ for r1 in range(Rmin_id, Rmax_id-1):
		# ~ for r2 in range(r1+1, Rmax_id):
	
			# ~ # Numerical integration using the trapezoidal rule
			# ~ delr = R_arr[r2] - R_arr[r1]
			# ~ fr1 = f(tau[r1, :], I[:, r1], R_arr[r1], delr)
			# ~ fr2 = f(tau[r2, :], I[:, r2], R_arr[r2], delr)
			# ~ summ += (fr1 + fr2)
	
	return summ / Jy
	
	
	
def plot_spectra(lamda, summ, solid_name, rv, fmax, Rmin, Rmax):
	
	fig = plt.figure()
	plt.plot(lamda, summ)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux [Jy]')
	plt.title(r'Spectrum {0} r={1} $\mu$m $f_{{max}}$={2} R={3}-{4} AU'.format(latex_name(solid_name), rv, fmax, Rmin, Rmax))
	plt.savefig("Spectra/Spectrum_{0}_r{1}_f{2}_R{3}-{4}.png".format(solid_name, rv, fmax, Rmin, Rmax))
	plt.show()



def hankel_transform(F_map, R_arr, wl, B):
	
	"""
	Performs Fourier transform in cylindrical coordinates to get the interferometric flux
	I (lsize, NPOINT)
	"""
	
	R_arr = R_arr * AU * 10**-2                 # Converting radii from AU to m
	# R_arr = R_arr[np.newaxis]
	# R_arr = R_arr.T
	
	# Considering F_map only at the given wavelength because I am out of ideas
	# lamda_rounded = np.round(lamda, 1)	
	# wl_id = np.where(lamda_rounded == wl)[0][0]
	# F_map = F_map[:, wl_id]
	
	# One more idea, summing up the wavelength axis
	F_map = np.sum(F_map, axis=1)
	
	bessel = j0(2.0 * np.pi * R_arr * B / (wl * 10**-6))
	# print(bessel)
	inter_flux = 2.0 * np.pi * np.trapz(R_arr * bessel * F_map, x = R_arr)
	# print(inter_flux)
	# print()
	
	return inter_flux
	


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
	T0 = 1500.0                             # Dust sublimation temperature (K)
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
	Rmax = 1.28 							# Maximum radius for spectrum plotting (AU) ENSURE IT IS ONLY 2 DECIMAL PLACES LONG
	gs = 0.1E-4                             # Grain radius (cm)
	wl = 5.5                                # Observing wavelength (microns)
	B = np.arange(20.0, 130.0, 2.0)          # 1D array of baselines (m)
	
	# All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
	minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])
	
	# Finding the most abundant condensates
	abundances, solid_names, abunds_dict = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
	top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names)
	top5_solids, topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds, abunds_dict)
	
	# Removing the solids without opfiles from top5_solids
	not_there = ['SiO', 'Mg3Si4O12H2', 'Fe3Si2O9H4', 'Ni', 'NaAlSi3O8', 'NaMg3AlSi3O12H2', 'CaAl2Si2O8']
	top5_solids = np.setdiff1d(top5_solids, not_there)
	
	# Calculating the surface density
	molwt = molecular_weight(top5_solids)
	surf_dens = surface_density(top5_solids, molwt, topabunds_radii, nHtot)
	
	# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
	opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64, 'Qcurve_inputs/Q_H2O_rv0.1_fmax0.7.dat' : 0.92}
	
	# Initializing dictionaries for the wavelengths, opacities, grain sizes and emptiness fractions used in the opfiles for each solid 	
	lamdas = {key: None for key in top5_solids}
	kappas = {key: None for key in top5_solids}
	rvs = {key: None for key in top5_solids}
	fmaxs = {key: None for key in top5_solids}
	
	for opfile, density in opfile_dens.items():
		mineral, rv, fmax, lamda, kappa = get_l_and_k(opfile, density, gs, lmin, lmax, lsize)
		# Qcurve_plotter(lamda, kappa, mineral, rv, fmax)
		lamdas[mineral] = lamda
		kappas[mineral] = kappa
		rvs[mineral] = rv
		fmaxs[mineral] = fmax		
	
	# Plotting the flux map and calculating the integrated flux for each solid
	I = {key: None for key in top5_solids}
	tau = {key: None for key in top5_solids}
	F_map_sum = np.zeros((NPOINT, lsize))
	intflux_sum = np.zeros(lsize)
	
	for solid in top5_solids:
			
			I[solid] = Plancks(T0, R_arr, R_in, lamdas[solid]) 
			tau[solid] = tau_calc(surf_dens[solid], kappas[solid])
			
			F_map = flux_map(tau[solid], I[solid])
			# plot_fluxmap(solid, rvs[solid], fmaxs[solid], F_map, lamdas[solid], R_arr)
			F_map_sum += F_map
			
			intflux = calculate_spectra(tau[solid], I[solid], R_arr, Rmin, Rmax)
			# plot_spectra(lamdas[solid], intflux, solid, rvs[solid], fmaxs[solid], Rmin, Rmax)
			intflux_sum += intflux
	
	# Plotting the overall flux map
	fig, ax = plt.subplots(1,1)
	img = ax.imshow(F_map_sum, cmap='plasma', interpolation='none')
	
	x_axis_locations = np.linspace(0, lsize-1, 8).astype(int)
	ax.set_xticks(x_axis_locations)
	x_axis_labels = np.round(lamda[x_axis_locations], 1)
	ax.set_xticklabels(x_axis_labels)
	ax.set_xlabel(r'$\lambda$ ($\mu$m)')
	
	y_axis_locations = np.linspace(0, len(R_arr)-1, 8).astype(int)
	ax.set_yticks(y_axis_locations)
	y_axis_labels = np.round(R_arr[y_axis_locations], 3)
	ax.set_yticklabels(y_axis_labels)
	ax.set_ylabel('R (AU)')
	
	ax.set_title(r"Overall Flux Map for r=0.1 microns")
	fig.colorbar(img)
	plt.savefig("Flux_Maps_trimmed/overall_fluxmap.png", bbox_inches = 'tight')
	plt.show()

	# Plotting the overall spectrum
	fig = plt.figure()
	plt.plot(lamdas['Mg2SiO4'], intflux_sum)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux [Jy]')
	plt.title(r'Overall Spectrum r=0.1 $\mu$m R={0}-{1} AU'.format(Rmin, Rmax))
	plt.savefig("Spectra/Overall_spectrum_r0.1_R{0}-{1}.png".format(Rmin, Rmax))
	plt.show()
	
	# Plotting the overall spectrum considering multiple radii together
	Rmin_list = [0.5, 0.75]
	Rmax_list = [1.0, 1.28]
	intflux_sum_mr = np.zeros(lsize)
	fig = plt.figure()
	
	for i in range(len(Rmax_list)):		
		for solid in top5_solids:		 			
			intflux_sum_mr += calculate_spectra(tau[solid], I[solid], R_arr, Rmin_list[i], Rmax_list[i])
			
		plt.plot(lamdas['Mg2SiO4'], intflux_sum_mr, label=r"($R_{{min}}$, $R_{{max}}$) = ({0},{1}) AU".format(Rmin_list[i], Rmax_list[i]))
		intflux_sum_mr = np.zeros(lsize)
			
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux [Jy]')
	plt.title(r'Overall spectrum for multiple radii - Zoomed')
	plt.legend()	
	plt.savefig("Spectra/Overall_spectrum_multiple_radii_limits_zoomed.png")
	plt.show()
	
	# Plotting the interferometric flux
	inter_flux = np.zeros(len(B))
	
	for bl in range(len(B)):
			
			inter_flux[bl] = hankel_transform(F_map_sum, R_arr, wl, B[bl])
		
	plt.plot(B, inter_flux)
	plt.xlabel(r'Baseline B (m)')
	plt.ylabel('Flux')
	plt.title(r'Interferometric flux $\lambda$ = {0} $\mu$m'.format(wl))
	plt.savefig("Spectra/Interferometric_flux_wl{0}.png".format(wl))
	plt.show()


if __name__ == "__main__":
    main()
