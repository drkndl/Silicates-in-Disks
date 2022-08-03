# This program calculates the properties such as pressure, density and number density of a planet forming disk that follows the model explained in Jorge et al. (2022)

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.constants import astropyconst20 as const
from HD163296_q04_p07_S4500_Tamf1350.properties import *

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
sb = const.sigma_sb.cgs 					# Stefan-Boltzmann constant in CGS units


def star_radius(L, T):
	
	"""
	Finding the stellar radius from the luminosity formula
	"""
	
	R = np.sqrt(L / (4 * np.pi * sb * T**4))
	
	return R
	
	
def inner_radius(Qr, T0, R_star, T_star):

    """
    Calculates the inner disk radius R_in (in AU) from the power law planet forming disk model

    Parameters:
    
    Qr:           Ratio of absorption efficiencies (float)
    T0:           Dust sublimation temperature in K (float)
    R_star:       Radius of the star in AU (float)
    T_star:       Effective temperature of the star in K (float)

    Returns:
    
    R_in: 		  The inner radius of the planet forming disk inside which dust grains cannot condense out from the gas (R_in) in AU (float)
    """

    R_in = 0.5 * np.sqrt(Qr) * (T_star/T0)**2 * R_star
    return R_in


def midplaneT_profile(R_in, T0, r_arr, q):

    """
    Calculates an array of midplane temperatures T (in K) in a disk to plots its radial dependence based on the given array of radii

    Parameters:
    
    R_in:        The inner radius of the protodisk inside of which the gas cannot condense into dust grains in AU (float)
    T0:          Dust sublimation temperature in K (float)
    r_arr:       A 1D array of radii in AU (shape: (100,)) (float)
    q: 			 Temperature gradient exponent (float)

    Returns:
    
    T_arr: 		 A 1D array of temperature in K calculated using the power law relationship (shape: (100,)) (float)
    """

    T_arr = T0 * (r_arr/R_in)**(q)
    return T_arr


def r_from_T(R_in, T_arr, T0, q):

    """
    Essentially the inverse of midplaneT_profile. Calculates an array of radii of the disk R_arr (in AU) for a given array of temperature based on the power-law relationship

    Parameters:
    
    R_in:        The inner radius of the protodisk inside of which the gas cannot condense into dust grains in AU (float)
    T_arr:       1D array of disk temperatures in K (shape: (NPOINT,)) (float) where NPOINT is the number of points in the GGchem simulation
    T0:          Dust sublimation temperature in K (float)

    Returns:
    
    R_arr: 		 A 1D array of radii in AU (shape: (NPOINT,)) (float)
    """

    R_arr = R_in * (T_arr/T0)**(1/q)
    return R_arr


def surface_density(Sigma0, r, p):

    """
    Calculates the surface density Sigma (in g/cm^2)

    Parameters:
    
    Sigma0:       Surface density with MMSN in g/cm^2 (float)
    r:            1D array of radii in AU (shape: (100,))
    p: 			  Surface density gradient exponent (float)

    Returns:
    
    Sigma: 		  1D array of surface density in g/cm^2 (shape: (100,)) (float)
    """

    Sigma = Sigma0 * (r / (1 * u.AU))**(p)
    return Sigma


def scale_height(M_star, r, T):

    """
    Calculates the pressure scale height H (in cm) by calculating the isothermal speed of sound and the Keplerian angular frequency
    
    Parameters:
    
    M_star: 		Mass of the star in g (float)
    r: 				1D array of radii in AU (shape varies) (float)
    T: 				1D array of temperature at each radial point in K calculated using the power law relationship (shape varies) (float)
    
    Returns:
    
    H: 				1D array of pressure scale height in cm calculated at each radial point (float)
    """

    omega = np.sqrt( G*M_star / (r.to(u.cm))**3 )   	# Keplerian angular frequency (/s)
    cs = np.sqrt( k*T / (mp*mu) )            			# Isothermal sound speed (cm/s)
    cs = cs.to(u.cm/u.s)
    H = cs / omega
 
    return H


def density(Sigma, H):
    
    """
    Calculates the mass density (g/cm^3) and number density (/cm^3) of the gas
    
    Parameters:
    
    Sigma: 		  1D array of surface density in g/cm^2 (shape: (100,)) (float)
    H: 			  1D array of pressure scale height in cm calculated at each radial point (float)
    
    Returns:
    
    rho: 		  1D array of mass density in g/cm^3 of the gas (float)
    nH: 		  1D array of the number density in /cm^3 of the gas (float)
    """

    rho = 0.5 * np.pi * Sigma / H
    nH = 0.7 * rho / mp                # 0.7 factor assumes the disk is made of 70% Hydrogen (from Planet_forming_disk_model.ipynb by Merijn)
    return rho, nH


def pressure(rho, T):
    
    """
    Calculates the gas pressure P (in bar)
    
    Parameters:
    
    rho: 		  1D array of mass density in g/cm^3 of the gas (float)
    T: 			  1D array of temperature at each radial point in K calculated using the power law relationship (shape varies) (float)
    
    Returns:
    
    P: 			  1D array of gas pressure in bar (float)
    """
    
    P = rho * Rc * T / (mu * Na * mp)
    return P.to(u.bar)


def main():
	
	R_star = star_radius(L_star, T_star).to(u.AU)   		# Star's radius (AU) 
	r_arr = np.linspace(0.05*u.AU, 2.5*u.AU, 100)
	R_in = inner_radius(Qr, T0, R_star, T_star)
	R_out = r_from_T(R_in, 100.0 * u.K, T0, q)
	# r_arr = np.linspace(R_in, R_out, 200)
	T_arr = midplaneT_profile(R_in, T0, r_arr, q)	
	
	Sigma = surface_density(Sigma0, r_arr, e)
	H = scale_height(M_star, r_arr, T_arr)
	
	rho, nH = density(Sigma, H)
	print("Number density nH: ", nH)

	P = pressure(rho, T_arr)
	# print("Pressure: ", P)
	
	# Plotting the temperature vs radius profile of the disk
	R_label = np.round(R_star/R_sun, 1)
	M_label = np.round(M_star/M_sun, 1)
	plt.plot(r_arr, T_arr)
	plt.xlabel("Radius R [AU]")
	plt.ylabel("Temperature T [K]")
	plt.title(r"$T_{{mid}}$ vs R, $R_{{star}}$ = {0}$R_\odot$, $T_{{star}}$ = {1} K, $M_{{star}}$ = {2}$M_\odot$, $\Sigma_0$ = {3} $g/cm^2$".format(R_label.value, T_star.value, M_label.value, Sigma0.value), fontsize=10)
	plt.savefig(folder + "Tmid_vs_R.png")
	plt.show()
	
	# Plotting the radial profile of pressure and density 
	plt.semilogy(r_arr, rho, label = r"Density $\rho$ [$gm/cm^3$]")
	plt.semilogy(r_arr, P, label = "Pressure [bar]")
	plt.xlabel("Radius R [AU]")
	plt.ylabel("Properties")
	plt.title(r"Radial dependence of $\rho$, P")
	plt.legend()
	plt.savefig(folder + "Pandrho_vs_R.png")
	plt.show()
	
	# Write the disk property values required for GGchem to a file
	with open(folder + 'disk_props.dat', 'w') as f:
		f.write('R_in' + '\t' + str(R_in) + '\n')
		f.write('\n')
		f.write('Prop \t Max \t Min \n')
		f.write('P' + '\t' + str(P.max()) + '\t' + str(P.min()) + '\n')
		f.write('nH' + '\t' + str(format(nH.max(),'.3E')) + '\t' + str(format(nH.min(),'.3E')) + '\n')
	

if __name__ == "__main__":
    main()
