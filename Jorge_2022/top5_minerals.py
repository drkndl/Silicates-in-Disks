# Code to find the top N most abundant minerals at each radius packet. Does not plot right now

import numpy as np
import matplotlib.pyplot as plt
from jorge_diskprop import inner_radius, r_from_T


def final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST):
    
    """
    Returns the abundances (log10 n_solid/n<H>) for each condensate

    Parameters:-
    keyword:        1D list of the Static_Conc.dat column heading names, used to identify the various solids
    minerals:       1D list of the condensates whose abundances are extracted
    dat:            File that contains the output (headers skipped)
    NELEM:          Total number of elements used
    NMOLE:          Number of molecules created
    NDUST:          Number of condensates created

    Returns abundances, a 2D array of shape (NPOINT, no.of minerals) and the corresponding solid names as a 1D list named solid_names
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


def most_abundant(top, NPOINT, abundances, R_arr, min_names):

    """
    Finds the top N most abundant minerals at each radial bin

    Parameters:-
    top:            Top X condensates whose abundance is the highest (integer)
    NPOINT:         Number of points in the simulation, also the length of R_arr (integer)
    abundances:     2D array of abundances for all the condensates present. Shape: (NPOINT, no.of condensates)
    R_arr:          1D array of derived radii from the temperature array in the output file
    min_names:      1D list of the condensate names

    Returns the 1D array of radial bins (R_bins), 2D array (shape: (top, NBins)) of the abundances of the top X most abundant species at each radial bin (top_abunds) and the corresponding solid names (top_solids), also as a 2D array of the same shape
    """

    indices = range(len(R_arr))                                                
                   
    to_calc = abundances[indices]                               # Abundances corresponding to radial bins. Shape is (NPOINT, no.of minerals)
    to_calc_sort = to_calc.copy()
    top_abunds = to_calc.copy()

    top_abunds.sort()                                           # Finding the most abundant condensates at each radial bin
    top_abunds = np.flip(top_abunds[:, -top:], axis=1)          # Making the sort descending order and saving the top 5
    
    idx = (-to_calc_sort).argsort(axis = -1)[:, :top]    

    # Finding the corresponding solid names
    top_solids = []
    for row in idx:
        for i in range(top):
            top_solids.append(min_names[row[i]])

    top_solids = np.array(top_solids)
    top_solids = np.reshape(top_solids, [NPOINT, top])

    return top_abunds, top_solids


def topabunds_by_radii(top_solids, solid_names, top_abunds):
    
    """
    Returns a dictionary of top 5 abundance values by radius. If the solid is not present in top 5 in that radius, the abundance is aken to be -300
    """

    # Make dictionary of all minerals to save radii at which the minerals are top 5 most abundant
    topabunds_radii = {key: np.full(500, -300.0) for key in solid_names}
    
    # Radii indices where the solids are top 5 most abundant    
    for solid in topabunds_radii.keys():
        
        idx = np.where(top_solids == solid)
        radii = idx[0]
        topabunds_radii[solid][radii] = top_abunds[idx]
        
    return topabunds_radii



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

    top = 5                 # Top X condensates whose abundance is the highest

    # All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
    minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si3O12', 'SNaAlSiO4', 'SCaAl2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O12F', 'SCa2MgSi2O7', 'SCa5P3O13H', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])

    # Fe based condensation sequences from Fig 4. Jorge et al. 2022:
    # minerals = ['SFe', 'SFeS', 'SMnS', 'SFe2SiO4', 'SFeAl2O4', 'SFeTiO3', 'SFeAl2SiO7H2', 'SKFe3AlSi3O12H2', 'SFe3O4', 'SFe3Si2O9H4', 'SNi3S2', 'SCa3Fe2Si3O12']
    
    abundances, solid_names = final_abundances(keyword, minerals, dat, NELEM, NMOLE, NDUST)
    top_abunds, top_solids = most_abundant(top, NPOINT, abundances, R_arr, solid_names)
    topabunds_radii = topabunds_by_radii(top_solids, solid_names, top_abunds)

    print("\n TOP", top, "HIGHEST ABUNDANCES AT EACH RADIAL BIN : \n", top_abunds)
    print("\n CORRESPONDING TOP", top, "SPECIES AT EACH RADIAL BIN: \n", top_solids)
    
    # Write down abundances and corresponding solids element by element in a file in a human readable format 
    filename = 'Sun/sun_most_abundant_readable.dat'
    with open(filename, 'w') as f:
        
        f.write('{0:47} {1:47} {2:47} {3:47} {4:47} Radius \n'.format(str(1), str(2), str(3), str(4), str(5)))                  # Adding headers with some formatting to help readability
        for radius in range(len(top_abunds)):
            for element in range(len(top_abunds[radius])):
                f.write('{0:20} : {1:20} '.format(str(top_solids[radius][element]), str(top_abunds[radius][element])))          # Adding the top 5 most abundant solids and their corresponding abundances
                if (element+1) % top == 0:
                    f.write(str(R_arr[radius]) + '\n')                                                                          # If all 5 solids are written, then add the radius in the last column and move to the next line
                else:
                    f.write('\t ')

    # Write down abundances and corresponding solids element by element in a file with just tab spaces for later parsing 
    filename2 = 'Sun/sun_most_abundant.dat'
    with open(filename2, 'w') as f:
        
        for radius in range(len(top_abunds)):
            for element in range(len(top_abunds[radius])):
                f.write('{0}:{1}'.format(str(top_solids[radius][element]), str(top_abunds[radius][element])))          # Adding the top 5 most abundant solids and their corresponding abundances
                if (element+1) % top == 0:
                    f.write('\t' + str(R_arr[radius]) + '\n')                                                          # If all 5 solids are written, then add the radius in the last column and move to the next line
                else:
                    f.write('\t')

    print("Abundances written to {}".format(filename))
    

if __name__ == "__main__":
    main()
