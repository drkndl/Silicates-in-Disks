import numpy as np
import pyvalem.formula as fo
from molmass import Formula
import astropy.units as u
from astropy.constants import astropyconst20 as const

Na = const.N_A.cgs                    		# Avogadro's number in /mol
solid = 'Mg0.7Fe0.3SiO3'
f = Formula(solid)
molwt = f.mass * u.g / u.mol / Na
print(molwt)

f = fo.Formula(solid)
fancy = "$" + f.latex + "$"
raw_solid = r"{}".format(fancy)
print(raw_solid)
