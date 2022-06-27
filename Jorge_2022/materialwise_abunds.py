# I live to make my own life harder every step of the way

import numpy as np

materials = ['Fe', 'Mg2SiO4', 'FeS', 'H20', 'Mg3Si2O9H4', 'SiO', 'MgSiO3', 'Mg3Si4O12H2', 'Fe3Si2O9H4', 'Fe2SiO4', 'Fe3O4', 'Ni', 'NaAlSi3O8', 'MgAl2O4', 'NaMg3AlSi3O12H2', 'CaMgSi2O6']
top = 5

with open('Sun/sun_most_abundant.dat', 'r') as f:

    top1, top2, top3, top4, top5, radius = ([] for i in range(top+1))

    for line in f:

        stuff = line.strip().split('\t')
        # top1.append(float(stuff[0]))
        # top2.append(float(stuff[1]))
        # top3.append(float(stuff[2]))
        # top4.append(float(stuff[3]))
        # top5.append(float(stuff[4]))
        # radius.append(float(stuff[-1]))

# print(top1, top3, radius)

# print(len(top1), len(top2), len(top3), len(top4), len(top5), len(radius))
