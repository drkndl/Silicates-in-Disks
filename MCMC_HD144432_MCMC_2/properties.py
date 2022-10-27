# HD144432
from astropy import units as u
from astropy.constants import astropyconst20 as const
import numpy as np 

T0 = 1500.0 * u.K                          				# Dust sublimation temperature (K)
Qr = 3                                  				# Ratio of absorption efficiencies 
L_star = 10**1.01 * const.L_sun.cgs         			# Stellar luminosity
T_star = 7345.138 * u.K                         		# Effective temperature of the star (K)
Sigma0 = 4500 * u.g / u.cm**2          					# Surface density with MMSN (g/cm^2)
M_star = 1.8 * 1.99E33 * u.g         					# Solar mass (g)
q = -0.373 												# Disk temperature gradient exponent
e = -0.7 												# Disk surface density gradient exponent
R_sun = 0.00465047      								# Sun's radius (AU)
M_sun = 1.99E33         								# Solar mass (g)
dist_pc = 145 * u.pc                           			# Star distance in parsec
H = 1.0 * u.cm 								    		# Scale height (cm)
add_gap = True                                         # True if adding a gap to the disk
rgap = 4.667263 * u.AU 										# The location of the gap in the disk (from the star) (AU)
wgap = 0.2016413 * u.AU 										# Total width of the gap (AU)
sgap = 10**-3  											# The amount by which the surface density is to be dampened in the gap
add_ring = True                                         # True if adding a gap to the disk
# rring = 2.543099 * u.AU 										# The location of the gap in the disk (from the star) (AU)
rring = 2.6 * u.AU
wring = 1.86739724 * u.AU 										# Total width of the gap (AU)
sring = 10**2  											# The amount by which the surface density is to be dampened in the gap

disk = 'MCMC_HD144432_MCMC_2' 										# The name of the disk
folder = disk + '/'  								    # Path where output files are saved
file = '{0}Static_Conc.dat'.format(folder)    			# Simulation output file
	
top = 5                                 	  			# Top X condensates whose abundance is the highest	
lmin = 0.0 * u.micron 						  			# Lower limit of wavelength (microns)
lmax = 20.0 * u.micron						  			# Upper limit of wavelength (microns)
lsize = 450 								  			# Number of wavelength (and kappa) points 
mass_fracs = {'Olivine': {'0.1': 0.185126159, '2.0': 0.814873},      # Mass fractions of the various silicates and grain sizes according to van Boekel (2005)
			'Pyroxene': {'0.1': 0.16969765, '2.0': 0.83030235},  
			# 'Mg2SiO4': {'0.1': 0.73, '2.0': 0.27},  
			'Mg2SiO4': {'0.1': 0.0113147494, '2.0': 0.9886852506},
			# 'MgSiO3': {'0.1': 0.4375, '2.0': 0.5625},
			'MgSiO3': {'0.1': 0.890703174, '2.0': 0.109296826},
			'Fe2SiO4': {'0.1': 0.0, '2.0': 0.0},
			'Fe3O4': {'0.1': 0.0, '2.0': 0.0},
			'Fe': {'0.1': 0.996336858, '2.0': 0.003663142},
			'FeS': {'0.1': 0.994501896, '2.0': 0.005498104},
			'Mg3Si2O9H4': {'0.1': 0.0, '2.0': 0.0},
			'MgAl2O4': {'0.1': 0.0, '2.0': 0.0},
			'CaMgSi2O6': {'0.1': 0.0, '2.0': 0.0}}
wl_list = [1.0, 2.0, 3.2, 5.5, 10.0, 12.0] * u.micron	# 1D list of wavelengths to plot correlated flux against baselines (microns)
B = np.arange(0.0, 130.0, 2.0) * u.m          			# 1D array of baselines (m)
B_small = np.linspace(0.0, 130.0, 5) * u.m    			# 1D array of a few baselines to plot correlated flux against wavelengths (m)
amor_temp = 1494.10890 * u.K 								# Temperature below which the grains are assumed to be amorphous
amor_temp_list = [400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0] * u.K
