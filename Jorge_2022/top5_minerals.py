import numpy as np
import matplotlib.pyplot as plt
from jorge_diskprop import inner_radius, r_from_T


def final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST):
    
    """
    Returns the abundances (log10 n_solil/n<H>) for each condensate
    """

    abundances = []
    solid_names = []
    for i in range(4+NELEM+NMOLE, 4+NELEM+NMOLE+NDUST, 1):
        
        solid = keyword[i]
        if keyword[i] in minerals:
            
            print(' i = ',i, ' solid name = ', solid)
            solid_names.append(solid[1:])
            ind = np.where(keyword == 'n' + solid[1:])[0]           # Finds the index where nsolid data for the current solid is available
            if (np.size(ind) == 0): continue
            ind = ind[0]
            abundances.append(dat[:, ind])                          # Saves the log10 nsolid/n<H> of the current solid

    abundances = np.array(abundances).transpose()                   # Transforming row of abundances into columns
    return abundances, solid_names


def most_abundant(abundances, R_arr, min_names):

    """
    Plots the top N most abundant minerals at each radial bin
    """

    NBins = 9
    top = 3
    
    # Geting evenly spaced radii in the array to define the radial bins
    indices = np.round(np.linspace(0, len(R_arr) - 1, NBins)).astype(int)
    R_bins = R_arr[indices]
    print(R_bins)
    # Note: Weird bins for now [0.0345255  0.04781386 0.06669758 0.0923685  0.12791976 0.17844072 0.24711991 0.34471807 0.47739495 0.66113721 0.9222484  1.27720817]

    to_calc = abundances[indices]
    print(to_calc)
    print("SHAPE ", np.shape(to_calc))
    to_calc_sort = to_calc.copy()
    top_abunds = to_calc.copy()

    top_abunds.sort()
    # top_abunds = top_abunds[:, -top:]
    top_abunds = np.flip(top_abunds[:, -top:], axis=1)          # Making the sort descending order
    print("TOP ABUNDS", top_abunds)
    indices = (-to_calc_sort).argsort(axis = -1)[:, :top]    
    print(indices)

    # Finding the corresponding solid names
    top_solids = []
    for row in indices:
        for i in range(top):
            top_solids.append(min_names[row[i]])

    top_solids = np.array(top_solids)
    top_solids = np.reshape(top_solids, [NBins, top])
    print(top_solids)

    return top_abunds, top_solids    

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

    bar   = 1.E+6                    # 1 bar in dyn/cm^2 
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

    # All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
    # minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])

    # Fe based condensation sequences from Fig 4. Jorge et al. 2022:
    minerals = ['SFe', 'SFeS', 'SMnS', 'SFe2SiO4', 'SFeAl2O4', 'SFeTiO3', 'SFeAl2SiO7H2', 'SKFe3AlSi3O12H2', 'SFe3O4', 'SFe3Si2O9H4', 'SNi3S2', 'SCa3Fe2Si3O12']
    
    abundances, solid_names = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
    print(np.shape(abundances), len(solid_names))
    most_abundant(abundances, R_arr, solid_names)
    

if __name__ == "__main__":
    main()
