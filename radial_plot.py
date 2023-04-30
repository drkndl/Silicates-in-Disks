# This program plots the condensation sequences as a function of radius instead of temperature using the power law relationship in Jorge et al. (2022)
# NOTE: CHANGE PLOT FORMATTING BASED ON YOUR NEED

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from diskprop import midplaneT_profile, r_from_T
from fancy_name import latex_name
# plt.rcParams['axes.linewidth'] = 1.5
# plt.rcParams["figure.figsize"] = (14, 7)


def R_plot(minerals, dat, keyword, R_arr, R_in, Rmin, Rmax, T0, q, folder, disk, NELEM, NMOLE, NDUST):
	
	# Some stylistic choices
	colo = ['blue', 'black', 'red', 'darkorange', 'green', 'darkorchid', 'aqua', 'cadetblue', 'cornflowerblue', 'chartreuse', 'limegreen', 'gold', 'chocolate', 'darkgoldenrod', 'darkolivegreen', 'darkmagenta', 'crimson', 'darkcyan', 'springgreen', 'darkslateblue', 'hotpink']
	#, 'bisque', 'indigo', 'peru', 'sienna', 'orangered', 'lightskyblue', 'navy', 'paleturquoise', 'deepskyblue', 'springgreen', 'plum', 'darkslateblue', 'mediumslateblue', 'goldenrod', 'gray', 'royalblue', 'cornflowerblue', 'lightcoral', 'rosybrown', 'saddlebrown', 'lime', 'forestgreen', 'lavender', 'hotpink', 'deeppink', 'gainsboro', 'peachpuff', 'beige']
	Ncolor = len(colo)
	colo = colo*10
	styl = ['solid']*Ncolor + ['dotted']*Ncolor + ['dashed']*Ncolor + ['dashdot']*Ncolor + ['(0, (1, 10))']*Ncolor*7
	widt = [3]*Ncolor*10
	
	ymin  = -5.4                # Minimum exponent on y-axis
	ymax = -4.4                  # Maximum exponent on y-axis
	csize = 5
	# Tlimit = 200 * u.K
	# T0 = 4000 * u.K
	# Rlimit = r_from_T(R_in, Tlimit, T0, q)                                            # Radius limit for zoomed in plot (AU)
	Rlimit = 5.0 * u.AU
	
	filename = folder + 'abundances_vs_R.png'
	
	points = np.where((R_arr>Rmin) & (R_arr<Rlimit))[0]             # Excluding the indices of the maximum and minimum gas temperature
	solids = []
	smean = []
	
	for i in range(4+NELEM+NMOLE,4+NELEM+NMOLE+NDUST,1):
		
		solid = keyword[i]
		if keyword[i] in minerals:
			
			# print(' i = ',i, ' solid name = ',solid)
			solids.append(solid[1:])                        # Saves the name of the current solid without the S at the beginning
			smean.append(np.mean(dat[points,i]))            # Calculates the mean of the above solid data (excluding the maximum and minimum temperatures)
	
	# pp = PdfPages(filename)
	
	# Creating the plots
	fig, ax = plt.subplots(figsize=(12, 6))
	ax2 = ax.twiny()                                        # Adding the temperatures as an X-axis on top of the plot
	indices = np.argsort(smean)
	
	colors = []
	count = 0
	
	print("Solids plotted: \n")
	for isolid in reversed(indices):                        # Going through the given solids in descending order of average log 10 nsolid/n<H>
		
		solid = solids[isolid]
		ind = np.where(keyword == 'n'+solid)[0]               # Finds the index where nsolid data for the current solid is available
		if (np.size(ind) == 0): continue
		ind = ind[0]
		yy = dat[:,ind]                                       # Saves the log10 nsolid/n<H> of the current solid
		
		if (np.max(yy[points])>ymin):
			print(solid)
			plt.plot(R_arr[points], yy[points], c = colo[count], ls = styl[count], lw = widt[count], label = latex_name(solid))
			colors.append(colo[count])
			count = count + 1
	
		
	# Plot formatting
	# plt.title('{0} Abundances of Condensates vs R'.format(disk), fontsize=14)
	ax.set_xlabel(r'$R\ \mathrm{[AU]}$', fontsize=19)
	ax.set_ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$', fontsize=19)
	ax.set_xlim(Rmin.value, Rlimit.value)
	ax.set_ylim(ymin, ymax)
	
	# Adding the temperature axis on top
	# Rmin = 0.01 * u.AU 
	x_ticks = np.linspace(Rmin, Rlimit, 7) 
	print(x_ticks)
	T_ticks = midplaneT_profile(R_in, T0, x_ticks, q).astype(int)
	print(T_ticks)
	ax2.set_xlim(ax.get_xlim())
	ax2.set_xticks(x_ticks.value)
	ax2.set_xticklabels(T_ticks.value)
	ax2.set_xlabel(r"$T \mathrm{[K]}$", fontsize=19)

	ax.tick_params(axis='x', labelsize=17)
	ax.tick_params(axis='y', labelsize=17)
	ax2.tick_params(axis='x', labelsize=17)
	
	leg = plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=19, fancybox=True, handlelength=0.5, prop={'size':csize}, ncol=1)       # Legend properties
	for color, text in zip(colors, leg.get_texts()):
		text.set_color(color)
		text.set_size(19)
	  
	plt.tight_layout()
	# plt.savefig(pp,format='pdf')
	plt.savefig(filename)
	plt.show()
	plt.clf()
	# pp.close()

