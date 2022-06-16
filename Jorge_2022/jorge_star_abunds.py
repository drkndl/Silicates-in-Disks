# This program calculates the star abundances of Jorge et al. (2022) : https://arxiv.org/pdf/2202.13920.pdf

import numpy as np
from solar_abundances import get_solar_abund, paper_solars


def star_abundances(solar_abundances, dex_dict):
    
    "This function calculates the custom stellar abundances given in Jorge et al. (2022) from the solar abundances using the dex values"

    star_abunds = {}

    for k in solar_abundances.keys() & dex_dict.keys():
        print(solar_abundances[k], dex_dict[k])
        star_abunds[k] = solar_abundances[k]*10**dex_dict[k]
        
    return star_abunds


def main():

    solar_abunds = get_solar_abund("../abundances.dat")
    paper_solar_abunds = paper_solars(solar_abunds)

    star1_dex = {'H': 0, 'He': 0, 'Li': 0.277, 'C': 0.225, 'N': 0.29, 'O': 0.25, 'F': 0.277, 'Na': 0.415, 'Mg': 0.15, 'Al': 0.19, 'Si': 0.285, 'P': 0.277, 'S': 0.277, 'Cl': 0.277, 'K': 0.277, 'Ca': 0.22, 'Ti': 0.285, 'V': 0.25, 'Cr': 0.27, 'Mn': 0.33, 'Fe': 0.39, 'Ni': 0.38, 'Zr': 0.277, 'W': 0.277}
    star1_abunds = star_abundances(paper_solar_abunds, star1_dex)
    print(star1_abunds)

    # Fractional values add up to approximately 1
    # print(sum(star1_abunds.values()))

    # Writing new star abundance values into an input file
    
    with open('Star1/star1_abund.in', 'w') as f:
        
        for key, value in star1_abunds.items():
            f.write(str(key) + '\t' + str(value) + '\n')

    
if __name__== "__main__":

    main()
