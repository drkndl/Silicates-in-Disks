import numpy as np
from diskprop import inner_radius, r_from_T, scale_height, surface_density
from astropy import units as u
from astropy.constants import astropyconst20 as const



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


def rho_g(Sg, Hp):
	rhog = Sg / (np.sqrt(2 * np.pi) * Hp) * np.exp(-1/2)
	return rhog
	

def omega(r, M_star):
	omega = np.sqrt( G*M_star / (r.to(u.cm))**3 )
	return omega
	
	
def isothermal_speed(T):
	cs = np.sqrt( k*T / (mp*mu) )
	return cs
	
	
def ubar(cs):
	u_bar = cs * np.sqrt(8 / np.pi)
	return u_bar
	
	
def stopping_time(rho_s, a, rhog, u_bar):
	u_bar = u_bar.to(u.cm / u.s)
	tau_s = rho_s * a / (rhog * u_bar)
	return tau_s
	
	
def Stokes(tau_s, omega_k):
	St = tau_s * omega_k
	return St
	
	
def dust_scale_height(St, Alpha, hg):
	
	"""
	Calculates the dust scale height from the gas scale height
	If this doesn't work I will pass away
	
	Edit: I will pass away now. Why did I think this would work?
	"""
	
	hd = (1 + St/Alpha)**(-1/2) * hg
	return hd
	

def main():
	
	file   = 'HotStar/HS_Static_Conc.dat'      # Simulation output file
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
	T0 = 1500.0 * u.K                          	# Dust sublimation temperature (K)
	Qr = 1                                  	# Ratio of absorption efficiencies (assumed to be black body)
	R_star = 2 * const.R_sun.to(u.AU)             # Star's radius (AU)
	T_star = 8000 * u.K                         # Effective temperature of the star (K)
	Sigma0 = 2 * 1700 * u.g / u.cm**2          		# Surface density with MMSN (g/cm^2)
	M_star = 8 * 1.99E33 * u.g         					# Solar mass (g)
	q = -0.75
	e = -1.5
	
	R_in = inner_radius(Qr, T0, R_star, T_star)   # Inner-most radius beyond which the dust is sublimated (AU)
	R_arr = r_from_T(R_in, Tg, T0, q)                # 1D array of radii obtained from the power law disk model (AU)
	H_gas = scale_height(M_star, R_arr, Tg)
	Sigma_g = surface_density(Sigma0, R_arr, e)
	
	hg = 1.15112697e+12 * u.cm
	sg = 314.09088132 * u.g / u.cm**2
	Alpha = 0.005
	R = 4.89345643 * u.AU
	rho_s = 3 * u.g / u.cm**3
	a = (0.1 * u.micron).to(u.cm)
	T = 100 * u.K
	
	rhog = rho_g(sg, hg)
	omega_k = omega(R, M_star)
	cs = isothermal_speed(T)
	u_bar = ubar(cs)
	tau_s = stopping_time(rho_s, a, rhog, u_bar)
	St = Stokes(tau_s, omega_k)
	hd = dust_scale_height(St, Alpha, hg)
	print(hd)
	
	
if __name__ == "__main__":
	main()
