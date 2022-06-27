# I live to make my own life harder every step of the way

import numpy as np

materials = ['Fe', 'Mg2SiO4', 'FeS', 'H20', 'Mg3Si2O9H4', 'SiO', 'MgSiO3', 'Mg3Si4O12H2', 'Fe3Si2O9H4', 'Fe2SiO4', 'Fe3O4', 'Ni', 'NaAlSi3O8', 'MgAl2O4', 'NaMg3AlSi3O12H2', 'CaMgSi2O6']

with open('Sun/sun_most_abundant.dat', 'r') as f:

    for line in f:

        stuff = line.strip()
        print(stuff)
