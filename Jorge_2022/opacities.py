import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from molmass import Formula
from pyvalem.formula import Formula
from jorge_diskprop import inner_radius, r_from_T
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii
from scipy.interpolate import UnivariateSpline
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
	Creates latex friendly names of compounds for better plot formatting
	
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

    solids       : A 1D list of the names (strings) of the top 5 most abundant solids at every radius 

    Returns a dictionary of the molecular weights (in g) for each of the given solids where the keys are the solid names and the values are the molecular weights
    """
    
    molwt = {key: None for key in solids}

    for solid in molwt.keys():
        f = Formula(solid)
        molwt[solid] = f.mass / Na
    
    return molwt



def surface_density(solids, molwt, top_abunds, nHtot):
    
    """
    Calculates the surface density of the given solids. Note that the calculated "surface" densities are actually (g/cm^3), but we assume column height to be 1 cm, so the surface density can be g/cm^2
    Shape of surf_dens = (NPOINT, 1)
    """

    # Transposing 1D row array (shape: (1, NPOINT)) into 2D column array (shape: (NPOINT, 1)) for matrix multiplication
    # ~ nHtot = nHtot[np.newaxis]
    # ~ nHtot = nHtot.T

    surf_dens = {key: None for key in solids}

    for solid in surf_dens.keys():
        n_solid = nHtot * 10**top_abunds[solid]
        surf_dens[solid] = molwt[solid] * n_solid
    
    return surf_dens



def Qcurve_plotter(opfile, dens, gs, lmin, lmax, lsize):
	
	"""
	Plots the Qcurve
	"""
	
	# Obtaining the mineral name, grain size value and fmax for plot labelling
	filename = opfile.split('/')[1]
	values = filename.split('_')
	mineral = values[1]
	
	for i in range(len(values)):
	
		if values[i].startswith('rv'):
			rv = values[i][2:5]
		elif values[i].startswith('fmax'):
			fmax = values[i][4:7]
	
	Q_curve = pd.read_csv(opfile, delimiter='\s+', skiprows=1, names=['wavelength','Q'])
	
	lamda = Q_curve['wavelength'].to_numpy()
	Q = Q_curve['Q'].to_numpy()
	
	# Taking only the indices where the wavelength is between 4 to 16 microns
	indices = np.where(np.logical_and(lamda >= lmin, lamda <= lmax))
	lamda = lamda[indices]
	Q = Q[indices]
	
	# Ensuring that there are an equal number of wavelength i.e. 600 points across all solids
	if len(lamda) > lsize:
		
		idx = np.round(np.linspace(0, len(lamda) - 1, lsize)).astype(int) 
		lamda = lamda[idx]
		Q = Q[idx]
	
	elif len(lamda) < lsize:
		
		old_indices = np.arange(0,len(lamda))
		new_indices = np.linspace(0, len(lamda)-1, lsize)
		
		spl1 = UnivariateSpline(old_indices, lamda, k=3, s=0)
		lamda = spl1(new_indices)
		
		spl2 = UnivariateSpline(old_indices, Q, k=3, s=0)
		Q = spl2(new_indices) 
	
	kappa = 3 * Q / (4 * gs * dens)
	
	plt.plot(lamda, kappa)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel(r'$\kappa_{abs}$ ($cm^2/g$)')
	plt.title(r"Q-curve for {0}, r = {1}, $f_{{max}}$ = {2}".format(latex_name(mineral), rv, fmax))
	plt.savefig("Qcurves_trimmed/Qcurve_{0}_r{1}_f{2}.png".format(mineral, rv, fmax), bbox_inches = 'tight')
	plt.show()
	
	return mineral, rv, fmax, lamda, kappa



def Plancks(T0, R_arr, R_in, lamda):
    
    """
    Calculates Planck function for a range of radius (converted from the temperature based on the disk model) and a range of wavelength in erg/(s cm^3 sr)
    """

    # Transposing 1D row array (shape: (1, NUMPOINTS)) into 2D column array (shape: (NUMPOINTS, 1)) for matrix multiplication
    # Note that NUMPOINTS, the number of wavelength points considered (which depends on the opacity file) is not the same as NPOINT, the number of points in the GGchem simuation
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
    Calculates opacity tau
    shape (r, lamda)
    """

    # Transposing 1D row array (shape: (1, NPOINT)) into 2D column array (shape: (NPOINT, 1)) for matrix multiplication
    sigma = sigma[np.newaxis]
    sigma = sigma.T

    tau = sigma * kappa
    
    return tau
    
    

def flux_map(solid_name, rv, fmax, tau, I, lamda, R_arr):

    """
    Plotting the flux map (r, lambda)
    """

    I = I.T                                             # Transposing I from (lamda, r) to (r, lamda) for element-wise matrix multiplication
    F_map = (1 - np.exp(-tau)) * I                      # Calculating the flux map

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
    
    return F_map



def f(tau, I, r, delr):
    
    """
    f in numerical integration to calculate the integrated flux for the s
    """

    return 0.5 * delr * AU * tau * I * 2*np.pi * r * AU

    
    
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
	T0 = 1500                               # Sublimation temperature (K)
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
	
	# All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
	minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])
	
	# Finding the most abundant condensates
	abundances, solid_names, abunds_dict = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
	top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names)
	top5_solids, topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds, abunds_dict)
	
	# Calculating the surface density
	molwt = molecular_weight(top5_solids)
	surf_dens = surface_density(top5_solids, molwt, topabunds_radii, nHtot)
	
	# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
	opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64}
	
	# Plotting the Qcurves	
	lamdas = {key: None for key in top5_solids}
	kappas = {key: None for key in top5_solids}
	rvs = {key: None for key in top5_solids}
	fmaxs = {key: None for key in top5_solids}
	
	for opfile, density in opfile_dens.items():
		mineral, rv, fmax, lamda, kappa = Qcurve_plotter(opfile, density, gs, lmin, lmax, lsize)
		lamdas[mineral] = lamda
		kappas[mineral] = kappa
		rvs[mineral] = rv
		fmaxs[mineral] = fmax		
	# Since some opacity files are missing at the moment, some lamda and kappa values are None
	
	# Plotting the flux map and calculating the integrated flux for each solid
	I = {key: None for key in top5_solids}
	tau = {key: None for key in top5_solids}
	F_map = np.zeros((NPOINT, lsize))
	intflux_sum = np.zeros(lsize)
	
	for solid in top5_solids:
		
		try: 
			
			I[solid] = Plancks(T0, R_arr, R_in, lamdas[solid]) 
			tau[solid] = tau_calc(surf_dens[solid], kappas[solid])
			F_map += flux_map(solid, rvs[solid], fmaxs[solid], tau[solid], I[solid], lamdas[solid], R_arr)
			intflux = calculate_spectra(tau[solid], I[solid], R_arr, Rmin, Rmax)
			plot_spectra(lamdas[solid], intflux, solid, rvs[solid], fmaxs[solid], Rmin, Rmax)
			intflux_sum += intflux

		except:
			
			TypeError
	
	
	# Plotting the overall flux map
	fig, ax = plt.subplots(1,1)
	img = ax.imshow(F_map, cmap='plasma', interpolation='none')
	
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
	Rmin_list = [0.03, 0.5, 0.75]
	Rmax_list = [0.035, 1.0, 1.28]
	# colors = ['blue', 'black', 'red', 'darkorange', 'gold', 'darkorchid', 'aqua', 'cadetblue', 'cornflowerblue', 'chartreuse', 'limegreen', 'darkgreen']
	colors = ['blue', 'darkorchid', 'aqua']
	intflux_sum_mr = np.zeros(lsize)
	fig = plt.figure()
	
	for i in range(len(Rmax_list)):		
		for solid in top5_solids:		
			try: 			
				intflux_sum_mr += calculate_spectra(tau[solid], I[solid], R_arr, Rmin_list[i], Rmax_list[i])
	
			except:			
				TypeError
			
		plt.plot(lamdas['Mg2SiO4'], intflux_sum_mr, color = colors[i], label=r"($R_{{min}}$, $R_{{max}}$) = ({0},{1}) AU".format(Rmin_list[i], Rmax_list[i]))
		intflux_sum_mr = np.zeros(lsize)
			
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux [Jy]')
	plt.title(r'Overall spectrum for multiple radii - Zoomed x2')
	plt.legend()	
	plt.savefig("Spectra/Overall_spectrum_multiple_radii_limits_zoomed2.png")
	plt.show()


if __name__ == "__main__":
    main()
