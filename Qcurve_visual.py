# This program is used to compare how the Qcurves for different materials vary based on the grain size distributions

import numpy as np
import matplotlib.pyplot as plt 
import glob
from astropy import units as u
import pandas as pd
from fancy_name import latex_name


def gimme_lk(filename):
	"""
	Returns the wavelength and opacity (Kabs) values in microns and cm^2/g respectively for a particular solid with a particular grain size distribution
	
	Parameters:
	
	filename 	  : The path to the opacity file (generated using optool) for a solid with a particular grain size distribution (string)
	
	Returns:
	
	lamda 		  : 1D array of the wavelength values in microns (float)
	kappa 		  : 1D array of the absorption opacities in cm^2/g (float)
	"""
	
	Q_curve = pd.read_csv(filename, delimiter='\s+', skiprows = 23, names = ['wavelength', 'Kabs', 'Ksca', 'g_asymm'])
	lamda = u.Quantity(Q_curve['wavelength'].to_numpy(), u.micron)
	kappa = Q_curve['Kabs'].to_numpy() * u.cm**2 / u.g
	
	return lamda, kappa
		
	
def main():
	
	folder = "Qcurve_Comparison/"
	minerals = ['Fe2SiO4', 'Fe3O4', 'Fe', 'FeS', 'Mg2SiO4', 'MgSiO3', 'Ni', 'SiO2']
	paths = ['Qcurve_inputs_mult_GS/0.1_to_1.5/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_10/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_20/*.dat', 'Qcurve_inputs_mult_GS/0.1_to_100/*.dat']
	gs_ranges = ['0.1_to_1.5', '0.1_to_10', '0.1_to_20', '0.1_to_100']
	
	# Defining a nested dictionary where each dictionary within is named after a mineral and each such dictionary consists of grain size ranges as keys with their wavelengths/opacities as values
	ldict = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in minerals}
	kdict = {outer_k: {inner_k: None for inner_k in gs_ranges} for outer_k in minerals}
	
	for path in paths:	
		
		# Obtaining the grain size range from the given path
		size = path.split('/')[1]
		for opfile in glob.glob(path):
			
			# Obtaining the mineral from the given path
			filename = opfile.split('/')[2]
			mineral = filename.split('_')[1]
			
			ldict[mineral][size], kdict[mineral][size] = gimme_lk(opfile)
	
	# Plot the opacities for only the first 20 microns	
	for mineral in minerals:
		for size in gs_ranges:
			
			# Take the last index where the wavelength value is less than or equal to 20 microns
			index = np.where(ldict[mineral][size] <= 20 * u.micron)[0][-1]
			ldict_plot = ldict[mineral][size][:index+1]
			kdict_plot = kdict[mineral][size][:index+1]
			
			# Plot the opacity curve of the mineral for a particular size range
			plt.plot(ldict_plot, kdict_plot, label = r'{0} $\mu$m'.format(size))
		
		plt.xlabel(r'$\lambda$ ($\mu$m)')
		plt.ylabel(r'$\kappa_{abs}$ ($cm^2/g$)')
		plt.title(r'Opacity curve for multiple grain sizes - {0}'.format(latex_name(mineral)))	
		plt.legend()
		plt.savefig(folder + "Qcurve_{0}_mult_size.png".format(mineral))
		plt.show()
	
	
if __name__ == "__main__":
	main()
	
