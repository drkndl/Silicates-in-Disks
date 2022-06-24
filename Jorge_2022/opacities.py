import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from molmass import Formula
from jorge_diskprop import inner_radius, r_from_T
from top5_minerals import final_abundances, most_abundant


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
    
    return molwt


def surface_density(molwt, top_abunds, nHtot):
    
    """
    Calculates the surface density of the given solids
    !!!!!!!!!!! NOTE: CURRENTLY THE DENSITY IS G/CM^3 !!!!!!!!!!!!!!!!!!!
    """

    # Transposing 1D column array (shape: (NPOINT, 1)) into 2D row array (shape: (1, NPOINT)) for matrix multiplication
    nHtot = nHtot[np.newaxis]
    nHtot = nHtot.T

    n_solid = nHtot * 10**top_abunds
    surf_dens = molwt * n_solid
    
    return n_solid, surf_dens


def Qcurve_plotter(file):

    """
    Plots the Qcurve
    List lengths 2001
    """

    Q_curve = pd.read_csv(file, delimiter='\s+', skiprows=1, names=['wavelength','kappa'])

    lamda = Q_curve['wavelength'].to_numpy()
    kappa = Q_curve['kappa'].to_numpy()

    # Taking only the first half of the data (don't need it till 200 microns)
    lamda = lamda[: len(lamda)//2]
    kappa = kappa[: len(kappa)//2]

    # Still need to trim the wavelength down to about 500 points, but I need it till 20 microns. So I'll take every alternate value
    lamda = lamda[::5]
    kappa = kappa[::5]
    
    print("KAPPA LAMDA", np.shape(kappa), np.shape(lamda))
    
    plt.plot(lamda, kappa)
    plt.xlabel(r'$\lambda$ ($\mu$m)')
    plt.ylabel(r'$\kappa_{abs}$ ($cm^2/g$)')
    plt.title(r"Q-curve for Fe, r = 0.1, $f_{max}$ = 0.7")
    plt.savefig("Qcurve_Fe_r0.1_f0.7.png")
    plt.show()

    return lamda, kappa


def Plancks(T0, R_arr, R_in, lamda, h, c, k):
    
    """
    Calculates Planck function for a range of radius (converted from the temperature based on the disk model) and a range of wavelength in J/(s m^3 sr)
    """
    
    # Transposing 1D column array (shape: (NPOINT, 1)) into 2D row array (shape: (1, NPOINT)) for matrix multiplication
    lamda = lamda[np.newaxis]
    lamda = lamda.T

    # !!!!!!!!!!!! SHOULD I DO THIS? Converting lamda to cm !!!!!!!!!!!!!!!!!!!!!!!
    lamda_cm = lamda * 10**3

    denominator = k * T0 * (R_arr/R_in)**(-3/4) * lamda_cm
    
    I = 2*h*c**2 / (lamda_cm**5 * np.exp( h*c/denominator) - 1)
    print("Planck function", np.shape(I))
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
    
    

def flux_map(tau, I, lamda, kappa, R_arr, surf_dens):

    """
    Plotting the flux map
    """

    # Transposing I from (lamda, r) to (r, lamda) for element-wise matrix multiplication
    I = I.T
    
    print("Tau", np.shape(tau))

    F_map = (1 - np.exp(-tau)) * I
    print("FMAP ", np.shape(F_map))

    fig, ax = plt.subplots(1,1)
    img = ax.imshow(F_map, origin='lower', interpolation='none')

    # Axes formatting !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SHOULD I TAKE LOGSPACE INSTEAD? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    x_axis_labels = np.arange(0, 15, 2)
    x_axis_locations = np.linspace(0, len(lamda), 8)
    ax.set_xticks(x_axis_locations)
    ax.set_xticklabels(x_axis_labels)
    ax.set_xlabel(r'$\lambda$ ($\mu$m)')

    y_axis_labels = np.round(np.linspace(0.03, 1.12, 8), 3)
    y_axis_locations = np.linspace(0, len(R_arr), 8)
    ax.set_yticks(y_axis_locations)
    ax.set_yticklabels(y_axis_labels)
    ax.set_ylabel('R (AU)')

    ax.set_title(r"Flux Map for Fe, r = 0.1, $f_{max}$ = 0.7")
    fig.colorbar(img)
    plt.savefig("Fe_flux_map.png", bbox_inches = 'tight')
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

    # Some constants in CGS 
    h = 6.6261e-27                   # Planck's constant in cm^2 g s-1
    c = 2.99792458e10                # Speed of light in cm/s                
    k = 1.3807e-16                   # Boltzmann constant in cm^2 g s^-2 K^-1
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

    top = 5                 # Top X condensates whose abundance is the highest

    # All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
    minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])

    # Finding the most abundant condensates
    abundances, solid_names = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
    top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names)

    # Calculating the surface density
    molwt = molecular_weight(top_solids)
    n_solid, surf_dens = surface_density(molwt, top_abunds, nHtot)

    # Plotting the Qcurves
    opfile = 'qval_fe_met_rv0.1fmax_0.7.dat'
    lamda, kappa = Qcurve_plotter(opfile)

    print('lamda before: ', lamda)

    # Plotting the flux map
    I = Plancks(T0, R_arr, R_in, lamda, h, c, k)
    print('lamda after: ', lamda)
    
    tau = tau_calc(surf_dens[:, 0], kappa)
    flux_map(tau, I, lamda, kappa, R_arr, surf_dens)
    

if __name__ == "__main__":
    main()
