# HD163296
from astropy import units as u
from astropy.constants import astropyconst20 as const
import numpy as np 

T0 = 1500.0 * u.K                          				# Dust sublimation temperature (K)
Qr = 3                                  				# Ratio of absorption efficiencies 
L_star = 10**1.38 * const.L_sun.cgs         			# Stellar luminosity
T_star = 8729.71 * u.K                         		# Effective temperature of the star (K)
Sigma0 = 1700 * u.g / u.cm**2          					# Surface density with MMSN (g/cm^2)
M_star = 2.0 * 1.99E33 * u.g         					# Solar mass (g)
q = -0.4 												# Disk temperature gradient exponent
e = -0.75 												# Disk surface density gradient exponent
R_sun = 0.00465047      								# Sun's radius (AU)
M_sun = 1.99E33         								# Solar mass (g)
dist_pc = 101 * u.pc                           			# Star distance in parsec (Calahan et al. 2021 MAPS)
H = 1.0 * u.cm 								    		# Scale height (cm)

disk = 'HD163296_q04_p075' 										# The name of the disk
folder = disk + '/'  								    # Path where output files are saved
file = '{0}{1}_Static_Conc.dat'.format(folder, disk)    # Simulation output file
	
top = 5                                 	  			# Top X condensates whose abundance is the highest	
lmin = 0.0 * u.micron 						  			# Lower limit of wavelength (microns) 
lmax = 20.0 * u.micron						  			# Upper limit of wavelength (microns)
lsize = 450 								  			# Number of wavelength (and kappa) points 
gs = 0.1 * u.micron                           			# Grain radius (cm)
wl_list = [1.0, 2.0, 3.2, 5.5, 10.0, 12.0] * u.micron	# 1D list of wavelengths to plot correlated flux against baselines (microns)
B = np.arange(0.0, 130.0, 2.0) * u.m          			# 1D array of baselines (m)
B_small = np.linspace(0.0, 130.0, 5) * u.m    			# 1D array of a few baselines to plot correlated flux against wavelengths (m)
amor_temp = 1000.0 * u.K 								# Temperature below which the grains are assumed to be amorphous
amor_temp_list = [400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1250.0, 1500.0] * u.K
