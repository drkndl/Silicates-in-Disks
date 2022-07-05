# This code is adapted from the plotting scripts in GGchem/tools and Plot_different_materials.ipynb by Merijn Lambregts. It is used to plot the condensation data by temperature for only the required solids.


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
from fancy_name import latex_name
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams["figure.figsize"] = (14, 7)


def T_plot(minerals, dat, keyword, Tg, Tmin, Tmax, folder, NELEM, NMOLE, NDUST):
	
	# Some stylistic choices    
	colo = ['blue', 'black', 'red', 'darkorange', 'gold', 'darkorchid', 'aqua', 'cadetblue', 'cornflowerblue', 'chartreuse', 'limegreen', 'darkgreen', 'chocolate', 'darkgoldenrod', 'darkkhaki', 'pink', 'moccasin', 'darkolivegreen', 'darkmagenta', 'aquamarine', 'coral', 'burlywood', 'silver', 'darkorange', 'crimson', 'darkcyan', 'bisque', 'indigo']
	#, 'peru', 'sienna', 'orangered', 'lightskyblue', 'navy', 'paleturquoise', 'deepskyblue', 'springgreen', 'plum', 'darkslateblue', 'mediumslateblue', 'goldenrod', 'gray', 'royalblue', 'cornflowerblue', 'lightcoral', 'rosybrown', 'saddlebrown', 'lime', 'forestgreen', 'lavender', 'hotpink', 'deeppink', 'gainsboro', 'peachpuff', 'beige']
	Ncolor = len(colo)
	colo = colo*10
	styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
	print(styl)
	widt = [2]*Ncolor*10
	
	ymin  = -12.0                # Minimum exponent on y-axis
	ymax = -2.0                  # Maximum exponent on y-axis
	csize = 5
	
	filename = folder + 'abundances_vs_T.pdf'
	
	points = np.where((Tg>Tmin) & (Tg<Tmax))[0]             # Excluding the indices of the maximum and minimum gas temperature
	solids = []
	smean = []
	
	for i in range(4+NELEM+NMOLE,4+NELEM+NMOLE+NDUST,1):
		
		solid = keyword[i]
		if keyword[i] in minerals:
			
			print(' i = ',i, ' solid name = ',solid)
			solids.append(solid[1:])                        # Saves the name of the current solid without the S at the beginning
			smean.append(np.mean(dat[points,i]))            # Calculates the mean of the above solid data (excluding the maximum and minimum temperatures)
	
	pp = PdfPages(filename)
	
	# Creating the plots
	fig, ax = plt.subplots()
	indices = np.argsort(smean)
	
	colors = []
	count = 0
	
	for isolid in reversed(indices):                        # Going through the given solids in descending order of average log 10 nsolid/n<H>
		
	    solid = solids[isolid]
	    ind = np.where(keyword == 'n'+solid)[0]               # Finds the index where nsolid data for the current solid is available
	    if (np.size(ind) == 0): continue
	    ind = ind[0]
	    yy = dat[:,ind]                                       # Saves the log10 nsolid/n<H> of the current solid
	  
	    if (np.max(yy[points])>ymin):
		
		    plt.plot(Tg[points], yy[points], c = colo[count], ls = styl[count], lw = widt[count], label = latex_name(solid))
		    colors.append(colo[count])
		    count = count + 1
	   
	# Plot formatting    
	plt.title('Abundances vs T', fontsize=20)
	plt.xlabel(r'$T\ \mathrm{[K]}$', fontsize=12)
	plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$', fontsize=12)
	plt.xlim(Tmin.value, Tmax.value)
	plt.ylim(ymin, ymax)
	plt.tick_params(bottom=True, top=True, left=True, right=True)
	plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
	plt.tick_params(axis='y',which='both',direction='in',length=7,right=True)
	plt.tick_params(axis='x',which='both',direction='in',length=7,right=True)
	plt.setp(ax.spines.values(), linewidth=1.5)
	ax.xaxis.set_tick_params(width=1.5)
	ax.yaxis.set_tick_params(width=1.5)
	
	leg = plt.legend(loc='upper right', fontsize=9, fancybox=True, handlelength=0.5, prop={'size':csize}, ncol=3)       # Legend properties
	for color, text in zip(colors, leg.get_texts()):		
	    text.set_color(color)
	    text.set_size(10)
	  
	plt.tight_layout()
	plt.savefig(pp,format='pdf')
	plt.clf()
	pp.close()

