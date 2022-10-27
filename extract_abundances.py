import numpy as np 


def change_C_by_O(abunds, ratio):

	"""
	Function to change C/O ratio of the given abundances (abunds) to the given ratio
	"""

	C = abunds['C']
	O = abunds['O']
	total = C + O 

	new_O = total / (1.0 + ratio)
	new_C = ratio * new_O

	abunds['C'] = new_C
	abunds['O'] = new_O

	return abunds


def change_Mg_by_Si(abunds, ratio):

	"""
	Function to change Mg/Si ratio of the given abundances (abunds) to the given ratio
	"""

	Mg = abunds['Mg']
	Si = abunds['Si']
	total = Mg + Si 

	new_Si = total / (1.0 + ratio)
	new_Mg = ratio * new_Si

	abunds['Mg'] = new_Mg
	abunds['Si'] = new_Si

	return abunds



def convert_abundances(abunds):

	"""
	Converts abundances to ratio of elemental particle densities with H nuclei particle density for input to GGchem
	"""

	for key, val in abunds.items():
		dens_ratio = 10.0**(val - 12.0)
		abunds[key] = dens_ratio

	return abunds 


def convert_abundances_back(dens_ratios):

	"""
	Converts abundances to ratio of elemental particle densities with H nuclei particle density for input to GGchem
	"""

	for key, val in dens_ratios.items():
		abund = np.log10(val) + 12.0
		dens_ratios[key] = abund

	return dens_ratios



def main():

	file = 'ggchem_abundances.dat'

	# Obtaining built-in solar abundances from GGchem
	solar_abundances = np.loadtxt(file, skiprows=1, usecols=[6])
	elements = np.loadtxt(file, skiprows=1, dtype=str, usecols=[2])
	# print(solar_abundances.shape, elements.shape)    # 83 elements

	# Create a dictionary out of these elemental abundances, which will be modified according to our needs
	ggchem_abunds = dict(zip(elements, solar_abundances))

	# Using the solar abundances from Asplund et al. (2009)
	asplund_abunds = {'H': 12.00, 'He': 10.93, 'Li': 1.05, 'Be': 1.38, 'B': 2.70, 'C': 8.43, 'N': 7.83, 'O': 8.69, 'F': 5.46, 'Ne': 7.93, 'Na': 6.24, 'Mg': 7.60, 'Al': 6.45, 'Si': 7.51, 'P': 5.41, 'S': 7.12, 'Cl': 5.50, 'Ar': 6.40, 'K': 5.03, 'Ca': 6.34, 'Sc': 3.15, 'Ti': 4.95, 'V': 3.93, 'Cr': 5.64, 'Mn': 5.43, 'Fe': 7.50, 'Co': 4.99, 'Ni': 6.22, 'Cu': 4.19, 'Zn': 4.56, 'Ga': 3.04, 'Ge': 3.65, 'As': 10**-12, 'Se': 10**-12, 'Br': 10**-12,'Kr': 3.25, 'Rb': 2.52, 'Sr': 2.87, 'Y': 2.21, 'Zr': 2.58}

	# Convert abundances to ratio of element density to H nuclei particle density
	asplund_converted = convert_abundances(asplund_abunds)

	# Change C/O ratio
	# C_by_O_ratio = 2.0
	# new_abundances = change_C_by_O(asplund_converted, C_by_O_ratio)

	# Change Mg/Si ratio to corroborate curves with Jorge et al. (2022)
	dex = -0.36
	Mg_by_Si_ratio = 10**dex
	new_abundances = change_Mg_by_Si(asplund_converted, Mg_by_Si_ratio)
	
	# Convert it back in terms of log(E_X) values
	new_abunds_asplund = convert_abundances_back(new_abundances)

	# Save updated abundances to file
	with open('HD144432_MgSi_036_log/MgSi_036_log.in', 'w') as f:
		for key in new_abundances.keys():
			f.write("%s\t%s\n" %(key, new_abundances[key]))


if __name__ == "__main__":
    main()