import numpy as np 
import matplotlib.pyplot as plt 
from astropy import units as u
from astropy.constants import astropyconst20 as const
from diskprop2 import midplaneT_profile, star_radius, inner_radius, r_from_T, surface_density, scale_height, density
from Temp.properties import *


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
	dg_ratio = 10**dat[:,ind].T[0, :]
	
	return dg_ratio
	
	
def dust_density(gasdens, dg_ratio):
	
	"""
	Returns the gas mass in the disk by multiplying Hydrogen mass to the Hydrogen particle nuclei density obtained by the GGchem simulation
	"""
	
	dustdens = dg_ratio * gasdens
	
	return dustdens
	
	
def dust_vol(dat, keyword):
	
	"""
	Calculates the dust volume per H nucleus from the GGchem output file
	"""
	
	ind = np.where(keyword == "dustVol/H")[0]
	dustVol = 10**dat[:,ind].T[0, :] * u.cm**3
	
	return dustVol
	
	
def masses(dustdens, dustvol, dg_ratio):
	
	dustmass = dustdens * dustvol
	gasmass = dustmass / dg_ratio
	
	return dustmass, gasmass
	
	
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
	
	# Extracting data from the output file
	Tg    = dat[:,0] * u.K                        # T [K]
	nHtot = dat[:,1] / u.cm**3                    # n<H> [cm^-3]          
	press = dat[:,2]                        	  # p [dyn/cm^2]
	Tmin  = np.min(Tg)                            # Minimum gas temperature
	Tmax  = np.max(Tg)                      	  # Maximum gas temperature	
	
	R_star = star_radius(L_star, T_star).to(u.AU)   		# Star's radius (AU)
	R_in = inner_radius(Qr, T0, R_star, T_star)   # Inner-most radius beyond which the dust is sublimated (AU)
	# R_arr = r_from_T(R_in, Tg, T0, q)           # 1D array of radii obtained from the power law disk model (AU)
	# print(np.max(R_arr))
	Rc = 1 * u.AU
	# ~ R_out = 10 * Rc
	R_out = r_from_T(R_in, 100.0 * u.K, T0, q)
	R_arr = np.linspace(R_in, R_out, 200)
	T_arr = midplaneT_profile(R_in, T0, R_arr, q)
	
	# plt.style.use('dark_background')
	
	plt.plot(R_arr, Tg, label = "GGchem", color='blue')
	plt.plot(R_arr, T_arr, label = "Calculated", color='red')
	plt.xlabel("Radius R [AU]")
	plt.ylabel("Temperature T [K]")
	plt.title("Temp Profile: GGchem output vs Calculated")
	plt.legend()
	plt.savefig(folder + "T_profile_GGchem_vs_calc.png")
	plt.show()
	
	q_GGchem = np.log(Tg / T0) / np.log(R_arr / R_in)
	# plt.plot(np.log(R_arr/R_in), np.log(Tg/T0))
	plt.plot(R_arr, np.log(Tg/T0) / np.log(R_arr/R_in))
	plt.xlabel(r"Radius R [AU]")
	plt.ylabel(r"$\frac{{log(T_g / T_0)}}{{log(R / R_{in})}}$")
	plt.title("q used by GGchem")
	plt.tight_layout()
	plt.savefig(folder + "q_GGchem.png")
	plt.show()
	print("q used by GGchem: ", q_GGchem)
	
	# plt.style.use('classic')
	
	T_with_ggchemq = midplaneT_profile(R_in, T0, R_arr, q_GGchem)
	plt.plot(R_arr, Tg, label = "GGchem", linewidth = 5)
	plt.plot(R_arr, T_with_ggchemq, label = "Calculated using GGchem q values")
	plt.xlabel("Radius R [AU]")
	plt.ylabel("Temperature T [K]")
	plt.title("Temp Profile: GGchem output vs Calculated with GGchem q")
	plt.legend()
	plt.savefig(folder + "T_profile_GGchem_vs_calcGGchemq.png")
	plt.show()
	
	scaleH = scale_height(M_star, R_arr, Tg)
	
	nHana = 0.5 * 0.7 * np.pi * Sigma0 * (R_arr / Rc)**e / (scaleH * mp)
	# nHana = np.pi * Sigma0 * (R_arr / (1 * u.AU))**e / (scaleH * mp)
	plt.semilogy(R_arr, nHana, label = "Calculated", color='red')
	plt.semilogy(R_arr, nHtot, label = "GGchem", color='blue')
	plt.xlabel("Radius R [AU]")
	plt.ylabel(r"nH [$/cm^3$]")
	plt.title("nH Profile: GGchem output vs Calculated")
	plt.legend()
	plt.savefig(folder + "nH_profile_GGchem_vs_calc.png")
	plt.show()
	
if __name__ == '__main__':
	main()
