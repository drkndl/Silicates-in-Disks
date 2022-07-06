# Returns a latex friendly name of a compound for better plot formatting

from pyvalem.formula import Formula

def latex_name(solid):
	
	"""
	Returns a string of latex friendly name of a compound for better plot formatting
	
	Parameters:
	
	solids         : The formula of the condensate (string)
	
	Example:
	
	>>> latex_name('CaMgSi2O6')
	>>> $CaMgSi_2O_6$
	"""
	
	f = Formula(solid)
	fancy = "$" + f.latex + "$"
	raw_solid = r"{}".format(fancy)
	
	return raw_solid