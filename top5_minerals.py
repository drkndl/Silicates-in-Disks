# Code to find the top N most abundant minerals at each radius packet. Does not plot right now

import numpy as np


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
			
			# print(' i = ',i, ' solid name = ', solid)
			solid_names.append(solid[1:])            
			ind = np.where(keyword == 'n' + solid[1:])[0]           # Finds the index where nsolid data for the current solid is available
			if (np.size(ind) == 0): continue
			ind = ind[0]
			abundances.append(dat[:, ind])                          # Saves the log10 nsolid/n<H> of the current solid
	
	abunds_dict = dict(zip(solid_names, abundances))
	abundances = np.array(abundances).transpose()                   # Transforming row of abundances into columns
	
	return abundances, solid_names, abunds_dict


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


def topabunds_by_radii(top_solids, solid_names, top_abunds, abunds_dict):
	
	"""
	Returns a dictionary of top 5 abundance values by radius. If the solid is not present in top 5 in that radius, the abundance is aken to be -300
	"""
	
	top5_solids = np.unique(top_solids)
	
	# ----------------------------- THE -300 IMPLEMENTATION ----------------------------------
	
	# ~ # Make dictionary of all minerals to save radii at which the minerals are top 5 most abundant
	# ~ topabunds_radii = {key: np.full(500, -300.0) for key in top5_solids}
	
	# ~ # Radii indices where the solids are top 5 most abundant    
	# ~ for solid in topabunds_radii.keys():
		
		# ~ idx = np.where(top_solids == solid)
		# ~ radii = idx[0]
		# ~ topabunds_radii[solid][radii] = top_abunds[idx]
		
	# ----------------------------- THE NON -300 IMPLEMENTATION -------------------------------
	
	topabunds_radii = {key: None for key in top5_solids}
	
	for solid in topabunds_radii.keys():		
		topabunds_radii[solid] = abunds_dict[solid]
		    
	return top5_solids, topabunds_radii


