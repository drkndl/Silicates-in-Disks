# import numpy as np
# import matplotlib.pyplot as plt
# import pyvalem.formula as fo
# from molmass import Formula
import astropy.units as u
from astropy.constants import astropyconst20 as const 
# from compare_grain_sizes import get_paper_spectra
# from scipy.interpolate import UnivariateSpline
# from HD163296_q04_p07_S3700_Tamf1150_mfchange.properties import *
# from no_thoughts_just_plots import add_textbox
# import matplotlib.font_manager
# fpaths = matplotlib.font_manager.findSystemFonts()

# disk = 'MCMC_HD144432' 										# The name of the disk
# folder = disk + '/'  								    # Path where output files are saved
# diskname = disk.split('_')[1]
# datfile = folder + 'van_Boekel_' + diskname + '.dat'
# wl, flux = get_paper_spectra(datfile)
# lsize = 450

# if len(wl) < lsize:
		
# 	# If there are less than lsize points of lamda, lsize points are created through interpolation and lamda is reassigned accordingly
# 	old_indices = np.arange(0, len(wl))
# 	new_indices = np.linspace(0, len(wl)-1, lsize)
	
# 	spl1 = UnivariateSpline(old_indices, wl, k=3, s=0)
# 	wl = spl1(new_indices)
	
# 	spl2 = UnivariateSpline(old_indices, flux, k=3, s=0)
# 	flux = spl2(new_indices)

# print(len(wl), len(flux))	
# plt.plot(wl, flux)
# plt.show()

# def sin(x, **kwargs):
	
# 	fig, ax = plt.subplots(1, 1)
# 	y = np.sin(x)
# 	plt.plot(x, y)
# 	textstr = add_textbox(q, e, Qr, Sigma0.value, amor_temp.value, add_gap, rgap.value, wgap.value, sgap)	
# 	plt.text(0.15, 0.8, textstr, transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 10, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))
# 	plt.tight_layout()
# 	plt.show()

# x = np.linspace(0, 1, 100)
# kwargs = {'q': q, 'e': e, 'Qr': Qr, 'Sigma0': Sigma0, 'amor_temp': amor_temp, 'add_gap': add_gap, 'rgap': rgap, 'wgap': wgap, 'sgap': sgap}
# sin(x, **kwargs)

# refjref

# for i in fpaths:
#     f = matplotlib.font_manager.get_font(i)
#     print(f.family_name)

# opfile = 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat'
# filename = opfile.split('/')[1]
# values = filename.split('_')
# mineral = values[1]

# for i in range(len(values)):

# 	if values[i].startswith('rv'):
# 		rv = values[i][2:5] 
# 	elif values[i].startswith('fmax'):
# 		fmax = values[i][4:7]
		
# gs = float(rv) * u.micron
# print(gs)

# mass_fracs = {'Olivine': {'0.1': 0.527, '2': 0},  'Pyroxene': {'0.1': 0, '2': 0.423},  'Mg2SiO4': {'0.1': 0.019, '2': 0.007},  'MgSiO3': {'0.1': 0.007, '2': 0.009}}
# olisilicates = ['Olivine', 'Mg2SiO4']
# pyrosilicates = ['Pyroxene', 'MgSiO3']

# for solid in pyrosilicates:
# 	other = np.setdiff1d(pyrosilicates, solid)[0]
# 	print(solid, other)
# 	solid_mass = (mass_fracs[solid]['0.1'] + mass_fracs[solid]['2']) / (mass_fracs[solid]['0.1'] + mass_fracs[solid]['2'] + mass_fracs[other]['0.1'] + mass_fracs[other]['2'])
# 	print(solid_mass)

# printfff
# # Checking Gauss for 
# r = np.linspace(0.1, 375, 3000)
# sigma = 3000 * r**(-0.75)
# print(np.round(r, 1))

# rgap = 0.7
# wgap = 0.7
# sgap = 10**-2
# # sgap = ((wgap/2.0) * np.sqrt(2 * np.pi))**-1

# rgap_ind = np.where(np.round(r, 1) == rgap)[0][0]
# print(rgap_ind)

# print(rgap - wgap/2.0)
# wgap_ind1 = np.where(np.round(r, 1) == np.round(rgap - wgap/2.0, 1))[0][0]
# wgap_ind2 = np.where(np.round(r, 1) == np.round(rgap + wgap/2.0, 1))[0][0]
# print(wgap_ind1, wgap_ind2)

# for i in range(wgap_ind1, wgap_ind2 + 1):
	
# 	if i <= rgap_ind:
		
# 		m = (sgap * sigma[rgap_ind] - sigma[wgap_ind1])/(r[rgap_ind] - r[wgap_ind1])
# 		sigma[i] = sigma[wgap_ind1] + m * (r[i] - r[wgap_ind1])
	
# 	else:
		
# 		m = (sigma[rgap_ind] - sigma[wgap_ind2])/(r[rgap_ind] - r[wgap_ind2])
# 		sigma[i] = sigma[rgap_ind] + m * (r[i] - r[rgap_ind])
		 
	
	
# 	# sigma[i] = 1 - (sgap * np.exp(-(r[i] - rgap)**2 / (2 * (wgap/2.0)**2)))
# 	# sigma[i] = 1 - (np.exp(-(r[i] - rgap)**2 / (2 * (wgap/2.0)**2)))	
# 	# sigma[i] = sigma[i] * sgap

# plt.plot(r, sigma)
# plt.show()

# printff

# disk = 'HD144432' 										# The name of the disk
# Folder = disk + '/'  								    # Path where output files are saved
# file = '{0}{1}_Static_Conc.dat'.format(Folder, disk)    # Simulation output file
# print(file)

# print(p)
# print(q)

# Na = const.N_A.cgs                    		# Avogadro's number in /mol
# solid = 'Mg0.7Fe0.3SiO3'
# f = Formula(solid)
# molwt = f.mass * u.g / u.mol / Na
# print(molwt)

# f = fo.Formula(solid)
# fancy = "$" + f.latex + "$"
# raw_solid = r"{}".format(fancy)
# print(raw_solid)

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

P = 10E-1 * u.bar
T = [1500, 100] * u.K
rho = P * (mu * Na * mp) / (Rc * T)
nH = rho / mp
print(nH)
