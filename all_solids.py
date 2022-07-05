# Finds and returns all the solid condensates formed in GGchem 

def get_all_solids(keyword, dat, NELEM, NMOLE, NDUST):
	
	minerals = []	
	for i in range(4+NELEM+NMOLE,4+NELEM+NMOLE+NDUST,1):
		
		if keyword[i].startswith('S') and not keyword[i].endswith('[l]'):
			
			print(' i = ',i, ' solid name = ', keyword[i])
			minerals.append(keyword[i])                        # Saves the name of the current solid
	
	return minerals

