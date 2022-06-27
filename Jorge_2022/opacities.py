import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from molmass import Formula
from jorge_diskprop import inner_radius, r_from_T
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii



# Some constants in CGS
Na = 6.022E23                    # Avogadro's number in /mol
h = 6.6261e-27                   # Planck's constant in cm^2 g s-1
c = 2.99792458e10                # Speed of light in cm/s                
k = 1.3807e-16                   # Boltzmann constant in cm^2 g s^-2 K^-1
bar   = 1.E+6                    # 1 bar in dyn/cm^2



def molecular_weight(solids):
    
    """
    Calculates the molecular weight of the given solids in g

    Parameters:

    solids       : A 2D array of the names of the top 5 most abundant solids at every radius (Shape (r, 5))

    Returns a 2D array of shape (r, 5) of the molecular weights (in g) of the given solids
    """

    shape = np.shape(solids) 
    molwt = np.zeros(shape)
    
    for r in range(len(solids)):
        for ele in range(len(solids[r])):
            f = Formula(solids[r][ele])
            molwt[r][ele] = f.mass

    molwt = molwt / Na
    
    return molwt



def surface_density(molwt, top_abunds, nHtot):
    
    """
    Calculates the surface density of the given solids. Note that the calculated "surface" densities are actually (g/cm^3), but we assume column height to be 1 cm, so the surface density can be g/cm^2
    Shape of surf_dens = (NPOINT, top)
    """

    # Transposing 1D column array (shape: (NPOINT, 1)) into 2D row array (shape: (1, NPOINT)) for matrix multiplication
    nHtot = nHtot[np.newaxis]
    nHtot = nHtot.T

    n_solid = nHtot * 10**top_abunds
    surf_dens = molwt * n_solid
    
    return n_solid, surf_dens



def Qcurve_plotter(opfile, dens, gs):
	
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
	
	# Taking only the first half of the data (don't need it till 200 microns)
	lamda = lamda[: len(lamda) // 2]
	Q = Q[: len(Q) // 2]
	
	# Still need to trim the wavelength down to about 500 points, but I need it till 20 microns. So I'll take every alternate value
	# lamda = lamda[::5]
	# Q = Q[::5]
	
	kappa = 3 * Q / (4 * gs * dens)
	
	plt.plot(lamda, kappa)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel(r'$\kappa_{abs}$ ($cm^2/g$)')
	plt.title(r"Q-curve for {0}, r = {1}, $fmax$ = {2}".format(mineral, rv, fmax))
	plt.savefig("Qcurves/Qcurve_{0}_r{1}_f{2}.png".format(mineral, rv, fmax), bbox_inches = 'tight')
	plt.show()
	
	return lamda, kappa



def Plancks(T0, R_arr, R_in, lamda):
    
    """
    Calculates Planck function for a range of radius (converted from the temperature based on the disk model) and a range of wavelength in erg/(s cm^3 sr)
    """
    
    # Transposing 1D column array (shape: (NPOINT, 1)) into 2D row array (shape: (1, NPOINT)) for matrix multiplication
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
    
    

def flux_map(tau, I, lamda, R_arr):

    """
    Plotting the flux map (r, lambda)
    """

    I = I.T                                             # Transposing I from (lamda, r) to (r, lamda) for element-wise matrix multiplication
    F_map = (1 - np.exp(-tau)) * I                      # Calculating the flux map

    fig, ax = plt.subplots(1,1)
    img = ax.imshow(F_map, cmap='jet', interpolation='none')

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

    ax.set_title(r"Flux Map for Forsterite, r = 0.1, $f_{max}$ = 1.0")
    fig.colorbar(img)
    plt.savefig("Forst_flux_map.png", bbox_inches = 'tight')
    plt.show()
    
    return F_map



def f(tau, I, r):
    
    """
    f in numerical integration to calculate the integrated flux for the s
    """

    return tau * I * 2*np.pi * r

    
    
def plot_spectra(tau, I, R_arr, lamda, Rmin, Rmax):
	
	"""
	Plots the integrated flux vs wavelength
	
	Summ shape (lamda, 1)
	"""
	
	# Finding the indices of Rmin and Rmax by taking the first instance of where they are in the rounded R_arr 
	R_rounded = np.round(R_arr, 2)	
	Rmin_id = np.where(R_rounded == Rmin)[0][0]
	Rmax_id = np.where(R_rounded == Rmax)[0][0]
	print(R_arr)
	print()
	print(R_rounded)
	print(Rmin_id, Rmax_id)
	print(R_rounded[Rmin_id], R_rounded[Rmax_id])
	print(R_arr[Rmin_id], R_arr[Rmax_id])
	
	summ = 0
	for r1 in range(Rmin_id, Rmax_id-1):
		for r2 in range(r1+1, Rmax_id):
	
			# Numerical integration using the trapezoidal rule
			delr = R_arr[r2] - R_arr[r1]
			fr1 = f(tau[r1, :], I[:, r1], R_arr[r1])
			fr2 = f(tau[r2, :], I[:, r2], R_arr[r2])
			summ += delr * 0.5 * (fr1 + fr2)
	
	fig = plt.figure()
	plt.plot(lamda, summ)
	plt.xlabel(r'$\lambda$ ($\mu$m)')
	plt.ylabel('Flux of some kind')
	plt.title(r'Spectrum $Mg_2SiO_4$ r=0.1 $\mu$m f=1.0 R={0}-{1} AU'.format(Rmin, Rmax))
	plt.savefig("Spectrum_Forst_r0.1_f1.0_R{0}-{1}.png".format(Rmin, Rmax))
	plt.show()
	
	return summ
    


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
	
	gs = 0.1E-4                             # Grain radius (cm)
	
	# All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
	minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])
	
	# Finding the most abundant condensates
	abundances, solid_names = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
	top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names)
	topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds)
	
	# Calculating the surface density
	molwt = molecular_weight(top_solids)
	n_solid, surf_dens = surface_density(molwt, top_abunds, nHtot)
	print("TOP ABUNDS: ", top_abunds)
	print()
	print("MOLECULAR WEIGHT: ", molwt)
	print()
	print("SURF_DENS: ", surf_dens)
	print()
	print("top abund, molwt, surf_dens shapes", top_abunds.shape, molwt.shape, surf_dens.shape)
	
	# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
	opfile_dens = {'Qcurve_inputs/Q_Diopside_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_Enstatite_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_Forsterite_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Serpentine_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_Spinel_rv0.1_fmax0.7.dat' : 3.64}
	
	# Plotting the Qcurves
	lamdas = []
	kappas = []
	
	for opfile, density in opfile_dens.items():
		lamda, kappa = Qcurve_plotter(opfile, density, gs)
	
	# ~ # Plotting the flux map
	# ~ I = Plancks(T0, R_arr, R_in, lamda)    
	# ~ tau = tau_calc(surf_dens[:, 1], kappa)
	# ~ F_map = flux_map(tau, I, lamda, R_arr)
	
	# ~ # Finding the integrated flux
	# ~ int_flux = plot_spectra(tau, I, R_arr, lamda, Rmin=0.03, Rmax=1.28)
	
	# ~ # TODO: Plotting fluxes of multiple radii together
	

if __name__ == "__main__":
    main()
