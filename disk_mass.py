import numpy as np 
import matplotlib.pyplot as plt 
from astropy import units as u
from astropy.constants import astropyconst20 as const
from HD144432_q04_p065_S4000_Tamf1420_massfrac.properties import *

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

def dust_to_gas(dat, keyword):
	
	"""
	Obtains the dust-to-gas ratio of the disk from the GGchem simulation output file
	"""
	
	ind = np.where(keyword == "dust/gas")[0]
	dg_ratio = dat[:,ind]
	
	return dg_ratio
	
def gas_mass(nH):
	
	"""
	Returns the gas mass in the disk by multiplying Hydrogen mass to the Hydrogen particle nuclei density obtained by the GGchem simulation
	"""
	
	gasmass = mp * nH
	return gasmass
	
	
def main():
	
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
	
	dg_ratio = dust_to_gas(dat, keyword)
	print(dg_ratio)
	
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
	
if __name__ == '__main__':
	main()
