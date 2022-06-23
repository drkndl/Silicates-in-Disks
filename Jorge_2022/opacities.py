import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from molmass import Formula
from jorge_diskprop import inner_radius, r_from_T
from top5_minerals import final_abundances, most_abundant


def molecular_weight(solids):
    
    """
    Calculates the molecular weight of the given solids in g
    """

    shape = np.shape(solids)    
    molwt = np.zeros(shape)
    
    for r in range(len(solids)):
        for ele in range(len(solids[r])):
            f = Formula(solids[r][ele])
            molwt[r][ele] = f.mass
    
    return molwt


def nH_in_Rbins(nHtot, R_arr, R_bins):
    
    """
    Finds the hydrogen particle density nH in the specific radial bins
    """

    R_arr = np.round(np.array(R_arr), 4)
    indices = []                                                # Finding the indices where the bin radii correspond to the radius array
    
    for r in range(len(R_arr)):
        if R_arr[r] in R_bins:
            indices.append(r)

    nH_rbins = nHtot[indices]

    return nH_rbins


def surface_density(molwt, top_abunds, nH_rbins):
    
    """
    Calculates the surface density of the given solids
    !!!!!!!!!!! NOTE: CURRENTLY THE DENSITY IS G/CM^3 !!!!!!!!!!!!!!!!!!!
    """

    # Transposing 1D column array (shape: (Nbins, 1)) into 2D row array (shape: (1, NBins)) for matrix multiplication
    nH_rbins = nH_rbins[np.newaxis]
    nH_rbins = nH_rbins.T

    n_solid = nH_rbins * 10**top_abunds
    surf_dens = molwt * n_solid

    # !!!!!!!!!!!!!!!!!! NOTE: CURRENTLY ADDING UP THE INDIVIDUAL SOLID DENSITIES AT EACH RADIAL BIN !!!!!!!!!!!!!!!!!!!!!!!!!!!
    surf_dens_sum = np.sum(surf_dens, axis=1)
    
    return n_solid, surf_dens_sum


def Plancks(T, lamda, h, c, k):
    
    """
    Calculates Planck function for an effective temperature and a range of wavelength
    I units J/(s m^3 sr)
    """

    # Converting lamda to m
    lamda = lamda * 10**6
    
    I = 2*h*c**2 / (lamda**5 * np.exp( h*c/(lamda*k*T) ) - 1)
    return I


def Qcurve_plotter(file):

    """
    Plots the Qcurve
    List lengths 2001
    """

    Q_curve = pd.read_csv(file, delimiter='\s+', skiprows=1, names=['wavelength','kappa'])

    lamda = Q_curve['wavelength'].to_numpy()
    kappa = Q_curve['kappa'].to_numpy()

    plt.plot(lamda, kappa)
    plt.show()

    return lamda, kappa


def tau_calc(sigma, kappa):
    
    """
    Calculates opacity tau
    shape (r, lamda)
    """

    # Transposing 1D row array (shape: (1, NBins)) into 2D column array (shape: (NBins, 1)) for matrix multiplication
    sigma = sigma[np.newaxis]
    sigma = sigma.T

    tau = sigma * kappa
    
    return tau
    
    

def flux_map(tau, I, lamda, kappa, R_bins, surf_dens):

    """
    Plotting the flux map
    """

    # Transposing 1D row array (shape: (1, 2001)) into 2D column array (shape: (2001, 1)) for matrix multiplication
##    I = I[np.newaxis]
##    I = I.T

    F_map = (1 - np.exp(-tau)) * I
    print(np.shape(F_map))

    temp = np.array([[1, 2], [4,3]])
    plt.imshow(temp)
    plt.colorbar()
    plt.show()
    
    return 
    


def main():

    file   = 'Sun/sun_Static_Conc.dat'          # Simulation output file
    data   = open(file)
    dummy  = data.readline()                    # Ignoring first line
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

    # Some constants
    h = 6.62607004e-34               # Planck's constant in Js
    c = 299792458.0                  # Speed of light in m/s
    k = 1.380649e-23                 # Boltzmann constant in J/K
    bar   = 1.E+6                    # 1 bar in dyn/cm^2

    # Extracting data from the output file
    Tg    = dat[:,0]                 # T [K]
    nHtot = dat[:,1]                 # n<H> [cm^-3]
    lognH = np.log10(nHtot)          
    press = dat[:,2]                 # p [dyn/cm^2]
    Tmin  = np.min(Tg)               # Minimum gas temperature
    Tmax  = np.max(Tg)               # Maximum gas temperature

    # Converting temperatures to corresponding radii
    T0 = 1500               # Sublimation temperature (K)
    Qr = 1                  # Ratio of absorption efficiencies (assumed to be black body)
    R_sun = 0.00465047      # Sun's radius (AU)
    T_sun = 5780            # Effective temperature of the sun (K)

    R_in = inner_radius(Qr, T0, R_sun, T_sun)
    R_arr = r_from_T(R_in, Tg, T0)

    NBins = 13              # Number of radial bins at which the most abundant species are extracted
    top = 5                 # Top X condensates whose abundance is the highest

    # All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
    minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])

    # Finding the most abundant condensates
    abundances, solid_names = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
    R_bins, top_abunds, top_solids = most_abundant(top, NBins, abundances, R_arr, solid_names)

    # Calculating the surface density
    molwt = molecular_weight(top_solids)
    nH_rbins = nH_in_Rbins(nHtot, R_arr, R_bins)
    n_solid, surf_dens = surface_density(molwt, top_abunds, nH_rbins)

    # Plotting the Qcurves
    opfile = 'qval_Fe2SiO4_rv2.0fmax_1.0.dat'
    lamda, kappa = Qcurve_plotter(opfile)

    # Plotting the flux map
    I = Plancks(T_sun, lamda, h, c, k)
    tau = tau_calc(surf_dens, kappa)
    flux_map(tau, I, lamda, kappa, R_bins, surf_dens)
    

if __name__ == "__main__":
    main()