# This program calculates the properties such as pressure, density and number density of a planet forming disk that follows the model explained in Jorge et al. (2022)

import numpy as np
import matplotlib.pyplot as plt


def inner_radius(Qr, T0, R_star, T_star):

    """
    Calculates the inner disk radius R_in (in AU) from the power law planet forming disk model

    Parameters:-
    Qr:           Ratio of absorption efficiencies (assumed to be black body)
    T0:           Dust sublimation temperature in K (float)
    R_star:       Radius of the star in AU (float)
    T_star:       Effective temperature of the star in K (float)

    Returns the inner radius of the planet forming disk inside which dust grains cannot condense out from the gas (R_in) in AU
    """

    R_in = 0.5 * np.sqrt(Qr) * (T_star/T0)**2 * R_star
    return R_in


def midplaneT_profile(R_in, T0, r_arr, q):

    """
    Calculates an array of midplane temperatures T (in K) in a disk to plots its radial dependence based on the given array of radii

    Parameters:-
    R_in:        The inner radius of the protodisk inside of which the gas cannot condense into dust grains in AU (float)
    T0:          Dust sublimation temperature in K (float)
    r_arr:       A 1D array of radii

    Returns a 1D array of temperature in K calculated using the power law relationship
    """

    T_arr = T0 * (r_arr/R_in)**(q)
    return T_arr


def r_from_T(R_in, T_arr, T0, q):

    """
    Essentially the inverse of midplaneT_profile. Calculates an array of radii of the disk R_arr (in AU) for a given array of temperature based on the power-law relationship

    Parameters:-
    R_in:        The inner radius of the protodisk inside of which the gas cannot condense into dust grains in AU (float)
    T_arr:       1D array of disk temperatures in K
    T0:          Dust sublimation temperature in K (float)

    Returns a 1D array of radii in AU 
    """

    R_arr = R_in * (T_arr/T0)**(1/q)
    return R_arr


def surface_density(Sigma0, r, e):

    """
    Calculates the surface density Sigma (in g/cm^2)

    Parameters:-
    Sigma0:       Surface density with MMSN in g/cm^2 (float)
    r:            1D array of radii

    Returns the surface density in g/cm^2 (float)
    """

    Sigma = Sigma0 * (r)**(e)
    return Sigma


def scale_height(G, M_star, r, kb, T, mp, mu):

    """
    Calculates the scale height H (in cm)
    """

    AU = 1.496E13                             # Astronomical unit (cm)
    omega = np.sqrt( G*M_star / (r*AU)**3 )   # Keplerian angular frequency (/s)
    cs = np.sqrt( kb*T / (mp*mu) )            # Isothermal sound speed (cm/s)
    H = cs / omega
    return H


def density(Sigma, H, mp):
    
    """
    Calculates the mass density (g/cm^3) and number density (/cm^3) of the gas
    """

    rho = 0.5 * np.pi * Sigma / H
    nH = rho / mp                # From Planet_forming_disk_model.ipynb
    return rho, nH


def pressure(rho, Rc, T, mu, N0, mp):
    
    """
    Calculates the gas pressure P (in bar)
    """

    bar = 10**-6                # Bar (dyne/cm^2)
    P = rho*Rc*T / (mu*N0*mp)
    return P*bar


def main():

    T0 = 1500               # Sublimation temperature (K)
    Tmax = 1500             # Maximum temperature in GGChem input (K)
    Tmin = 100              # Minimum temperature in GGChem input (K)
    Qr = 1                  # Ratio of absorption efficiencies (assumed to be black body)
    q = -0.6
    e = -0.1
    R_star = 2*0.00465047   # Star's radius (AU)
    R_sun = 0.00465047      # Sun's radius (AU)
    T_star = 8000           # Effective temperature of the star (K)
    Sigma0 = 2*1700         # Surface density with MMSN (g/cm^2)
    G = 6.6725E-8           # Gravitational constant (cm3 g^-1 s^-2)
    M_star = 8*1.99E33      # Star mass (g)
    M_sun = 1.99E33         # Solar mass (g)
    kb = 1.3806E-16 	    # Boltzmann constant (erg k^-1)
    mp = 1.6733E-24         # Mass of proton (g)
    mu = 2.3                # Mean molecular weight
    Rc = 8.314E7            # Ideal gas constant (erg K^-1 mol^-1)
    N0 = 6.022E23           # Avogadro's constant (/mol)
    Folder = "Temp/"
    
    r_arr = np.linspace(0.05, 2.5, 100)      # AU

    R_in = inner_radius(Qr, T0, R_star, T_star)
    print(R_in)
    T_arr = midplaneT_profile(R_in, T0, r_arr, q)

    Sigma = surface_density(Sigma0, r_arr, e)
    H = scale_height(G, M_star, r_arr, kb, T_arr, mp, mu)
    
    rho, nH = density(Sigma, H, mp)
    # print("Number density nH: ", nH)
    
    P = pressure(rho, Rc, T_arr, mu, N0, mp)
    # print("Pressure: ", P)

    # Plotting the temperature vs radius profile of the disk
    R_label = np.round(R_star/R_sun, 1)
    M_label = np.round(M_star/M_sun, 1)
    plt.plot(r_arr, T_arr)
    plt.xlabel("Radius R [AU]")
    plt.ylabel("Midplane Temperature T [K]")
    plt.title(r"$T_{{mid}}$ vs R, $R_{{star}}$ = {0}$R_\odot$, $T_{{star}}$ = {1} K, $M_{{star}}$ = {2}$M_\odot$, $\Sigma_0$ = {3} $g/cm^2$".format(R_label, T_star, M_label, Sigma0), fontsize=10)
    plt.savefig(Folder + "Tmid_vs_R.png")
    plt.show()

    # Plotting the radial profile of pressure and density 
    plt.semilogy(r_arr, rho, label = r"Density $\rho$ [$gm/cm^3$]")
    plt.semilogy(r_arr, P, label = "Pressure [bar]")
    plt.xlabel("Radius R [AU]")
    plt.ylabel("Properties")
    plt.title(r"Radial dependence of $\rho$, P")
    plt.legend()
    plt.savefig(Folder + "Pandrho_vs_R.png")
    plt.show()

    # Write the disk property values required for GGchem to a file
    with open(Folder + 'disk_props.dat', 'w') as f:
        f.write('Prop \t Max \t Min \n')
        f.write('P' + '\t' + str(P.max()) + '\t' + str(P.min()) + '\n')
        f.write('nH' + '\t' + str(nH.max()) + '\t' + str(nH.min()) + '\n')


if __name__ == "__main__":
    main()
