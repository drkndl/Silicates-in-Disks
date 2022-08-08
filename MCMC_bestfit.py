# Tutorial by Imad Pasha (https://prappleizer.github.io/Tutorials/MCMC/MCMC_Tutorial.html)


import numpy as np
import emcee
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.constants import astropyconst20 as const
from scipy.interpolate import UnivariateSpline
from spectra import molecular_weight, surface_density, Plancks, tau_calc, tau_calc_amorphous, flux_map, calculate_spectra, get_l_and_k
from all_solids import get_all_solids
from molmass import Formula
from diskprop import inner_radius, r_from_T, scale_height, star_radius
from T_plot import T_plot
from radial_plot import R_plot
from compare_grain_sizes import get_paper_spectra
from top5_minerals import final_abundances, most_abundant, topabunds_by_radii


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


################################################### Defining all the non-parameters and other arguments for the MCMC ###########################################################################################

T0 = 1500.0 * u.K                          				# Dust sublimation temperature (K)
Qr = 3                                  				# Ratio of absorption efficiencies 
q = -0.373
e = -0.7
Sigma0 = 4500.0 * u.g/u.cm**2
L_star = 10**1.01 * const.L_sun.cgs         			# Stellar luminosity
T_star = 7345.138 * u.K                         		# Effective temperature of the star (K)
M_star = 1.8 * 1.99E33 * u.g         					# Solar mass (g)
R_sun = 0.00465047      								# Sun's radius (AU)
M_sun = 1.99E33         								# Solar mass (g)
dist_pc = 145 * u.pc                           			# Star distance in parsec
H = 1.0 * u.cm 								    		# Scale height (cm)
# rc = 1 * u.AU
add_gap = True                                         # True if adding a gap to the disk
add_ring = False
sgap = 10**-3  											# The amount by which the surface density is to be dampened in the gap
sring = 10**-3  											# The amount by which the surface density is to be dampened in the gap

disk = 'MCMC_HD144432' 										# The name of the disk
folder = disk + '/'  								    # Path where output files are saved
file = '{0}Static_Conc.dat'.format(folder)    			# Simulation output file
	
top = 5                                 	  			# Top X condensates whose abundance is the highest	
lmin = 0.0 * u.micron 						  			# Lower limit of wavelength (microns)
lmax = 20.0 * u.micron						  			# Upper limit of wavelength (microns)
lsize = 450 								  			# Number of wavelength (and kappa) points 
R_star = star_radius(L_star, T_star).to(u.AU)   		# Star's radius (AU)
	
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
# press = dat[:,2]                        	  # p [dyn/cm^2]
Tmin  = np.min(Tg)                            # Minimum gas temperature
Tmax  = np.max(Tg)                      	  # Maximum gas temperature	

R_in = inner_radius(Qr, T0, R_star, T_star)   # Inner-most radius beyond which the dust is sublimated (AU)
R_arr = r_from_T(R_in, Tg, T0, q)             # 1D array of radii obtained from the power law disk model (AU)

# H = scale_height(M_star, R_arr, Tg)

minerals = get_all_solids(keyword, dat, NELEM, NMOLE, NDUST)

# Finding the most abundant condensates
abundances, solid_names, abunds_dict = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST) 
top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names) 
top5_solids, topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds, abunds_dict)

# Removing the solids without opfiles from top5_solids
not_there = ['SiO', 'Mg3Si4O12H2', 'Fe3Si2O9H4', 'Ni', 'NaAlSi3O8', 'NaMg3AlSi3O12H2', 'CaAl2Si2O8', 'H2O', 'Ca2MgSi2O7', 'NH3', 'Al2O3', 'Ca2Al2SiO7', 'Ca3Al2Si3O12', 'CaAl2Si2O8', 'ZrO2', 'Ti3O5', 'W', 'VO', 'CaTiO3', 'NaAlSiO4']
top5_solids = np.setdiff1d(top5_solids, not_there)
top5_solids = np.concatenate((top5_solids, ['Olivine', 'Pyroxene']))

# Calculating molecular weights
molwt = molecular_weight(top5_solids)

# Creating a dictionary of Qcurve input files and the corresponding material densities in g/cm^3
opfile_dens = {'Qcurve_inputs/Q_CaMgSi2O6_rv0.1_fmaxxxx.dat' : 3.278, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv0.1.dat' : 3.2, 'Qcurve_inputs/Q_MgSiO3_Jaeger_DHS_fmax1.0_rv2.0.dat' : 3.2, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv0.1.dat' : 3.27, 'Qcurve_inputs/Q_Mg2SiO4_Sogawa_DHS_fmax1.0_rv2.0.dat' : 3.27, 'Qcurve_inputs/qval_Fe3O4_rv0.1_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe3O4_rv2.0_fmax0.7.dat' : 5.17, 'Qcurve_inputs/qval_Fe2SiO4_rv0.1_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe2SiO4_rv2.0_fmax1.0.dat' : 4.392, 'Qcurve_inputs/qval_Fe_met_rv0.1_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_Fe_met_rv2.0_fmax0.7.dat' : 7.874, 'Qcurve_inputs/qval_FeS_rv0.1_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_FeS_rv2.0_fmax0.7.dat' : 4.84, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv0.1_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_Mg3Si2O9H4_rv2.0_fmax0.7.dat' : 2.6, 'Qcurve_inputs/qval_MgAl2O4_rv0.1_fmax0.7.dat' : 3.64, 'Qcurve_inputs/qval_MgAl2O4_rv2.0_fmax0.7.dat' : 3.64, 'Qcurve_inputs/Q_Olivine_rv0.1_fmax0.7.dat': 3.71, 'Qcurve_inputs/Q_Olivine_rv2.0_fmax0.7.dat' : 3.71, 'Qcurve_inputs/Q_Pyroxene_rv0.1_fmax0.7.dat': 3.01, 'Qcurve_inputs/Q_Pyroxene_fmax0.7_rv2.0.dat': 3.01}

# Adding units to the material densities using astropy
for key, value in opfile_dens.items():
	opfile_dens[key] = value * u.g / u.cm**3

# Initializing dictionaries for the wavelengths, opacities, grain sizes and emptiness fractions used in the opfiles for each solid 	
gs_ranges = ['0.1', '2.0']
lamdas = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
kappas = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
fmaxs = {key: None for key in top5_solids}

for opfile, density in opfile_dens.items():
	mineral, rv, fmax, lamda, kappa = get_l_and_k(opfile, density, lmin, lmax, lsize)		
	# Qcurve_plotter(lamda, kappa, mineral, rv, fmax, folder)
	lamdas[mineral][rv] = lamda
	kappas[mineral][rv] = kappa
	fmaxs[mineral] = fmax
	
diskname = disk.split('_')[1]
datfile = folder + 'van_Boekel_' + diskname + '.dat'
wl, flux = get_paper_spectra(datfile)

if len(wl) != lsize:
		
	# If there are less than lsize points of lamda, lsize points are created through interpolation and lamda is reassigned accordingly
	old_indices = np.arange(0, len(wl))
	new_indices = np.linspace(0, len(wl)-1, lsize)
	
	spl1 = UnivariateSpline(old_indices, wl, k=3, s=0)
	wl = spl1(new_indices)
	
	spl2 = UnivariateSpline(old_indices, flux, k=3, s=0)
	flux = spl2(new_indices)

flux_err = 0.05 * np.mean(flux)
# data = (q, e, Sigma0, Tg, R_in, R_arr, H, nHtot, molwt, top5_solids, add_gap, add_ring, sgap, sring, lamdas, kappas, topabunds_radii, dist_pc, flux, flux_err)
data = (lamdas['Mg2SiO4']['0.1'], flux, flux_err)
nwalkers = 24
niter = 3000
initial = np.array([0.5, 0.2, 0.5, 0.2, 1350.0, 1.0, 0.0, 0.25, 0.35, 0.5, 0.5])
ndim = len(initial)
p0 = [np.array(initial) + 1e-2 * np.random.randn(ndim) for i in range(nwalkers)]


###############################################################################################################################################################################################################


def model(theta, q=q, p=e, Sigma0=Sigma0, T=Tg, Rin=R_in, R_arr=R_arr, H=H, nHtot=nHtot, mol=molwt, solids=top5_solids, add_gap=add_gap, add_ring=add_ring, sgap=sgap, sring=sring, lamdas=lamdas, kappas=kappas, topabunds=topabunds_radii, dist_pc=dist_pc):
	
	# Model parameters
	# rgap, wgap, rring, wring, amor_temp, mf_ol_01, mf_py_01, mf_fors_01, mf_ens_01, mf_Fe_01, mf_FeS_01 = theta
	rgap_ul, wgap_ul, rring_ul, wring_ul, amor_temp_ul, mf_ol_01, mf_py_01, mf_fors_01, mf_ens_01, mf_Fe_01, mf_FeS_01 = theta
	print(rgap_ul, wgap_ul, amor_temp_ul, mf_ol_01, mf_py_01, mf_fors_01, mf_ens_01, mf_Fe_01, mf_FeS_01)
	
	# Adding units to the model
	rgap = rgap_ul * u.AU
	wgap = wgap_ul * u.AU
	rring = rring_ul * u.AU
	wring = wring_ul * u.AU
	amor_temp = amor_temp_ul * u.K
	
	# Defining mass fractions in a dictionary
	mass_fracs = {'Olivine': {'0.1': mf_ol_01, '2.0': 1.0-mf_ol_01},      
			'Pyroxene': {'0.1': mf_py_01, '2.0': 1.0-mf_py_01},    
			'Mg2SiO4': {'0.1': mf_fors_01, '2.0': 1.0-mf_fors_01},
			'MgSiO3': {'0.1': mf_ens_01, '2.0': 1.0-mf_ens_01},
			'Fe2SiO4': {'0.1': 0.0, '2.0': 0.0},
			'Fe3O4': {'0.1': 0.0, '2.0': 0.0},
			'Fe': {'0.1': mf_Fe_01, '2.0': 1.0-mf_Fe_01},
			'FeS': {'0.1': mf_FeS_01, '2.0': 1.0-mf_FeS_01},
			'Mg3Si2O9H4': {'0.1': 0.0, '2.0': 0.0},
			'MgAl2O4': {'0.1': 0.0, '2.0': 0.0},
			'CaMgSi2O6': {'0.1': 0.0, '2.0': 0.0}}
			
	# Calculate some disk parameters
	Rmin = np.round(np.min(R_arr), 1) 						# Minimum radius for spectrum plotting (AU)
	Rmax = np.round(np.max(R_arr), 1)						# Maximum radius for spectrum plotting (AU)
	
	# Calculate the surface densities of the individual solids
	if add_gap:
		surf_dens, rgap_ind, wgap_ind1, wgap_ind2 = surface_density(top5_solids, mol, topabunds, nHtot, H, add_gap, R_arr, rgap, wgap, sgap, add_ring, rring, wring, sring)
	elif add_ring:
		surf_dens, rring_ind, wring_ind1, wring_ind2 = surface_density(top5_solids, mol, topabunds, nHtot, H, add_gap, R_arr, rgap, wgap, sgap, add_ring, rring, wring, sring)
	else:
		surf_dens = surface_density(top5_solids, mol, topabunds, nHtot, H, add_gap, R_arr, rgap, wgap, sgap, add_ring, rring, wring, sring)
	
	# Calculate spectra
	gs_ranges = ['0.1', '2.0']
	I = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	tau = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	F_map_sum = np.zeros((NPOINT, lsize)) * u.erg / (u.s * u.Hz * u.sr * u.cm**2)
	intflux = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in top5_solids}
	intflux_sum = np.zeros(lsize) * u.Jy
	
	for solid in top5_solids:
		for size in gs_ranges:
			
			if lamdas[solid][size] == None:
				# print("Skipping: ", solid, size)
				continue
		
			if solid == 'Olivine':
				continue
			
			elif solid == 'Pyroxene':
				continue
						
			elif solid == "Mg2SiO4":
					
				I[solid][size] = Plancks(lamdas[solid][size], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid][size] = tau_calc_amorphous(surf_dens[solid], surf_dens['Olivine'], kappas[solid][size], kappas['Olivine'][size], Tg, amor_temp, mass_fracs[solid][size])
				# plot_tau(tau[solid][size], solid, size, folder)
				
				F_map = flux_map(tau[solid][size], I[solid][size])
				# plot_fluxmap(solid, size, fmaxs[solid], F_map, lamdas[solid][size], R_arr, folder, **kwargs) 
				F_map_sum += F_map
				
				intflux[solid][size] = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid][size], intflux[solid][size], solid, size, fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux[solid][size]
				
			elif solid == "MgSiO3":
					
				I[solid][size] = Plancks(lamdas[solid][size], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid][size] = tau_calc_amorphous(surf_dens[solid], surf_dens['Pyroxene'], kappas[solid][size], kappas['Pyroxene'][size], Tg, amor_temp, mass_fracs[solid][size])
				# plot_tau(tau[solid][size], solid, size, folder)
				
				F_map = flux_map(tau[solid][size], I[solid][size])
				# plot_fluxmap(solid, size, fmaxs[solid], F_map, lamdas[solid][size], R_arr, folder, **kwargs) 
				F_map_sum += F_map
				
				intflux[solid][size] = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid][size], intflux[solid][size], solid, size, fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux[solid][size]
				
			else:
				
				I[solid][size] = Plancks(lamdas[solid][size], Tg) 
				# plot_Bv(lamdas[solid], I[solid], solid, folder)
				
				tau[solid][size] = tau_calc(surf_dens[solid], kappas[solid][size], mass_fracs[solid][size])
				# plot_tau(tau[solid][size], solid, size, folder)
				
				F_map = flux_map(tau[solid][size], I[solid][size])
				# plot_fluxmap(solid, size, fmaxs[solid], F_map, lamdas[solid][size], R_arr, folder, **kwargs)
				F_map_sum += F_map
				
				intflux[solid][size] = calculate_spectra(F_map, R_arr, Rmin, Rmax, dist_pc)  
				# plot_spectra(lamdas[solid][size], intflux[solid][size], solid, size, fmaxs[solid], Rmin, Rmax, folder)
				intflux_sum += intflux[solid][size]
				
	return intflux_sum.value
	
	
def lnlike(theta, x, y, yerr):
    return -0.5 * np.sum(((y - model(theta, x))/yerr) ** 2)
    
    
def lnprior(theta):
	
	rgap, wgap, rring, wring, Tamf, mf_ol_01, mf_py_01, mf_fors_01, mf_ens_01, mf_Fe_01, mf_FeS_01 = theta
	if 0.4 <= rgap <= 5.0 and 0.2 <= wgap <= 2.0 and 0.4 <= rring <= 5.0 and 0.2 <= wring <= 2.0 and 500.0 <= Tamf <= 1500.0 and 0.0 <= mf_ol_01 <= 1.0 and 0.0 <= mf_py_01 <= 1.0 and 0.0 <= mf_fors_01 <= 1.0 and 0.0 <= mf_ens_01 <= 1.0 and 0.0 <= mf_Fe_01 <= 1.0 and 0.0 <= mf_FeS_01 <= 1.0:
		return 0.0
	return -np.inf

	
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
    
    
def main(p0,nwalkers,niter,ndim,lnprob,data):
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)

    print("Running burn-in...")
    p0, _, _ = sampler.run_mcmc(p0, 300, progress=True)
    sampler.reset()

    print("Running production...")
    pos, prob, state = sampler.run_mcmc(p0, niter, progress=True)

    return sampler, pos, prob, state


sampler, pos, prob, state = main(p0,nwalkers,niter,ndim,lnprob,data)
samples = sampler.flatchain

theta_max  = samples[np.argmax(sampler.flatlnprobability)]
print(theta_max)
best_fit_model = model(theta_max)

diskname = disk.split('_')[0]
datfile = folder + 'van_Boekel_' + diskname + '.dat'
wl, flux = get_paper_spectra(datfile)

if add_gap:
	
	textstr = '\n'.join((
				r'$q=%.3f$' % (q, ),
				r'$p=%.2f$' % (e, ),
				r'$Q_r=%d$' % (Qr, ),
				r'$\Sigma_0=%.1f g/cm^2$' % (Sigma0.value, ),
				r'$T_{amf}=%.1f K$' % (theta_max[4], ),
				r'$r_{gap}=%.1f AU$' % (theta_max[0], ),
				r'$w_{gap}=%.1f AU$' % (theta_max[1], ),
				r'$s_{gap}=%.4f$' % (sgap, ),
				r'$mf_{Ol,0.1}=%.3f$' % (theta_max[5], ),
				r'$mf_{Py,0.1}=%.3f$' % (theta_max[6], ),
				r'$mf_{Fors,0.1}=%.3f$' % (theta_max[7], ),
				r'$mf_{Ens,0.1}=%.3f$' % (theta_max[8], ),
				r'$mf_{Fe,0.1}=%.3f$' % (theta_max[9], ),
				r'$mf_{FeS,0.1}=%.3f$' % (theta_max[10], )))
				
elif add_ring:
	
	textstr = '\n'.join((
				r'$q=%.3f$' % (q, ),
				r'$p=%.2f$' % (e, ),
				r'$Q_r=%d$' % (Qr, ),
				r'$\Sigma_0=%.1f g/cm^2$' % (Sigma0.value, ),
				r'$T_{amf}=%.1f K$' % (theta_max[4], ),
				r'$r_{ring}=%.1f AU$' % (theta_max[2], ),
				r'$w_{ring}=%.1f AU$' % (theta_max[3], ),
				r'$s_{ring}=%.4f$' % (sring, ),
				r'$mf_{Ol,0.1}=%.3f$' % (theta_max[5], ),
				r'$mf_{Py,0.1}=%.3f$' % (theta_max[6], ),
				r'$mf_{Fors,0.1}=%.3f$' % (theta_max[7], ),
				r'$mf_{Ens,0.1}=%.3f$' % (theta_max[8], ),
				r'$mf_{Fe,0.1}=%.3f$' % (theta_max[9], ),
				r'$mf_{FeS,0.1}=%.3f$' % (theta_max[10], )))
				
else:
	
	textstr = '\n'.join((
				r'$q=%.3f$' % (q, ),
				r'$p=%.2f$' % (e, ),
				r'$Q_r=%d$' % (Qr, ),
				r'$\Sigma_0=%.1f g/cm^2$' % (Sigma0.value, ),
				r'$T_{amf}=%.1f K$' % (theta_max[4], ),
				r'No gap or ring',
				r'$mf_{Ol,0.1}=%.3f$' % (theta_max[5], ),
				r'$mf_{Py,0.1}=%.3f$' % (theta_max[6], ),
				r'$mf_{Fors,0.1}=%.3f$' % (theta_max[7], ),
				r'$mf_{Ens,0.1}=%.3f$' % (theta_max[8], ),
				r'$mf_{Fe,0.1}=%.3f$' % (theta_max[9], ),
				r'$mf_{FeS,0.1}=%.3f$' % (theta_max[10], )))

plt.plot(lamdas['Mg2SiO4']['0.1'], best_fit_model, color="black", label="Best Fit")
plt.plot(wl, flux, color="grey", label = "Data")
plt.text(0.15, 0.85, textstr, transform=axs[0].transAxes, horizontalalignment='center', verticalalignment='center', fontsize = 10, bbox = dict(boxstyle='round', facecolor = 'white', alpha = 0.5))
plt.xlabel(r"Wavelength $\lambda$ ($\mu$m)")
plt.ylabel(r"Flux (Jy)")
plt.title("Spectrum with MCMC Best Fit")
plt.legend()
plt.savefig(folder + "bets_fit_spectrum.png")
plt.show()
