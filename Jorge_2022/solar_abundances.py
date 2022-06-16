# Saves solar fractional abundances from Abundances.dat



def get_solar_abund(all_abunds):
    
    "This function saves only the solar abundances of the 92 elements from the .dat file of all abundances"

    solar_values = {}
    
    with open(all_abunds, mode='r') as file:
        header = next(file)
        
        for line in file:
            words = line.strip().split()
            # print(words, "\n", len(words))
            # solar_values.update( {words[2], words[6]} )
            solar_values[words[2]] = float(words[6])

    return solar_values



def paper_solars(solar_values):

    "This function saves only the solar abundances of the elements used in the paper Jorge et al. (2022) : https://arxiv.org/pdf/2202.13920.pdf"

    paper_elements = ['H', 'He', 'Fe', 'C', 'O', 'Mg', 'Si', 'Ca', 'Ti', 'Li', 'N', 'F', 'Na', 'Al', 'P', 'Cl', 'K', 'V', 'Cr', 'Mn', 'Ni', 'Zr', 'W', 'S']
    print(len(paper_elements))
    
    for key, value in solar_values.copy().items():
        
        if key not in paper_elements:
            del solar_values[key]

    return solar_values
        


# Note that these solar abundances are in fractional form where the sum adds up to 1

# solar_abunds = get_solar_abund("abundances.dat")
# total_solar = sum(solar_abunds.values())
# print(total_solar)

# Outputs 1.0000422738109862


            
        
