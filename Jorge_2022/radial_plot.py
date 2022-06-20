# This program plots the condensation sequences as a function of radius instead of temperature using the power law relationship in Jorge et al. (2022)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
from jorge_diskprop import inner_radius, r_from_T, midplaneT_profile
from pyvalem.formula import Formula
plt.rcParams['axes.linewidth'] = 1.5


def main():
    
    file   = 'Sun/sun_Static_Conc.dat'          # Output file
    data   = open(file)
    dummy  = data.readline()                # Ignoring first line
    dimens = data.readline()                
    dimens = np.array(dimens.split())
    NELEM  = int(dimens[0])                 # Total number of elements used
    NMOLE  = int(dimens[1])                 # Number of molecules created
    NDUST  = int(dimens[2])                 # Number of condensates created
    NPOINT = int(dimens[3])                 # Number of points in simulation
    header = data.readline()                # Saves parameter names such as molecule names
    data.close()

    dat = np.loadtxt(file,skiprows=3)
    keyword = np.array(header.split())      # Array of parameter names such as molecule names

    bar   = 1.E+6                    # 1 bar in dyn/cm^2 
    Tg    = dat[:,0]                 # T [K]
    nHtot = dat[:,1]                 # n<H> [cm^-3]
    lognH = np.log10(nHtot)          
    press = dat[:,2]                 # p [dyn/cm^2]
    Tmin  = np.min(Tg)               # Minimum gas temperature
    Tmax  = np.max(Tg)               # Maximum gas temperature

    # Converting temperatures to corresponding radii
    T0 = 1500               # Sublimation temperature (K)
    Qr = 1                  # Ratio of absorption efficiencies (assumed to be black body)
    R_sun = 0.00465047      # Sun's radius (AU)
    T_sun = 5780            # Effective temperature of the sun (K)

    R_in = inner_radius(Qr, T0, R_sun, T_sun)
    R_arr = r_from_T(R_in, Tg, T0)
    Rmin = np.min(R_arr)    # Minimum radius
    Rmax  =np.max(R_arr)    # Maximum radius

    # Some stylistic choices
    colo = ['blue','black', 'red', 'darkorange', 'gold', 'darkorchid', 'aqua', 'cadetblue', 'cornflowerblue', 'chartreuse', 'limegreen', 'darkgreen', 'chocolate', 'darkgoldenrod', 'darkkhaki', 'pink', 'moccasin', 'darkolivegreen', 'darkmagenta', 'aquamarine', 'coral', 'burlywood', 'silver', 'beige', 'darkorange', 'crimson', 'darkcyan', 'bisque']
    Ncolor = len(colo)
    colo = colo*10
    styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
    widt = [2]*Ncolor
  
    ymin  = -12.0                # Minimum exponent on y-axis
    ymax = -4.0                  # Maximum exponent on y-axis
    csize = 5

    # All 52 condensates for Sun from Fig C.1 Jorge et al. 2022:
    # minerals = (['SZrSiO4', 'SV2O3', 'SCaTiSiO5', 'SCr2O3', 'SCaMgSi2O6', 'SMg2SiO4','SMgSiO3','SMg3Si2O9H4', 'SMgCr2O4', 'SMnTiO3', 'SNi', 'SFe', 'SZrO2', 'SFeS', 'SCa3Al2Si2O8', 'SMgAl2O4', 'SFeTiO3', 'SMnS', 'SNaAlSi3O8', 'SW', 'SCaTiO3', 'SMn3Al2Si3O12', 'SKAlSi3O8', 'SNi3S2', 'SNaCl', 'SVO', 'SFeAl2O4', 'SAlO2H', 'SFe2SiO4', 'SCa5P3O13H', 'SCa2MgSiO7', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2', 'SLi2SiO3', 'SWO3', 'SLiCl', 'SMg3Si4O12H2', 'SMnAl2SiO7H2', 'SFeAl2SiO7H2', 'SFe3O4', 'SCa3Fe2Si3O12', 'STi3O5', 'STi4O7', 'SSiO', 'SKFe3AlSi3O12H2', 'SCr', 'SMg3Si2O9H4', 'SCaAl2Si2O10H4', 'SH2O', 'SFe3Si2O9H4'])

    # Fe based condensation sequences from Fig 4. Jorge et al. 2022:
    minerals = ['SFe', 'SFeS', 'SMnS', 'SFe2SiO4', 'SFeAl2O4', 'SFeTiO3', 'SFeAl2SiO7H2', 'SKFe3AlSi3O12H2', 'SFe3O4', 'SFe3Si2O9H4', 'SNi3S2', 'SCa3Fe2Si3O12']

    # Mg based condensation sequences from Fig 3. Jorge et al. 2022:
    # minerals = ['SMg2SiO4', 'SMgSiO3', 'SSiO', 'SMgCr2O4', 'SCaMgSi2O6', 'SMgAl2O4', 'SCa2MgSi2O7', 'SMg3Si2O9H4', 'SMg3Si4O12H2', 'SKMg3AlSi3O12H2', 'SNaMg3AlSi3O12H2']
    
    filename = 'Sun/sun_Fe_condensates_radial.pdf'

    points = np.where((R_arr>Rmin) & (R_arr<Rmax))[0]             # Excluding the indices of the maximum and minimum gas temperature
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
    ax2 = ax.twiny()                                        # Adding the temperatures as an X-axis on top of the plot
    indices = np.argsort(smean)
    
    colors = []
    count = 0
    
    for isolid in reversed(indices):                        # Going through the given solids in descending order of average log 10 nsolid/n<H>
        
      solid = solids[isolid]
      ind = np.where(keyword == 'n'+solid)[0]               # Finds the index where nsolid data for the current solid is available
      if (np.size(ind) == 0): continue
      ind = ind[0]
      yy = dat[:,ind]                                       # Saves the log10 nsolid/n<H> of the current solid

      # Latex formula for plot legend
      f = Formula(solid)
      fancy = "$" + f.latex + "$"
      raw_solid = r"{}".format(fancy)
      
      if (np.max(yy[points])>ymin):
        plt.plot(R_arr[points], yy[points], c = colo[count], ls = styl[count], lw = widt[count], label = raw_solid)
        colors.append(colo[count])
        count = count + 1

        
    # Plot formatting
    Rlimit = 0.6                                            # Radius limit for zoomed in plot (AU)
    plt.title('Fe Condensates (Sun) vs Radius', fontsize=12)
    ax.set_xlabel(r'$R\ \mathrm{[AU]}$', fontsize=10)
    ax.set_ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$', fontsize=10)
    # ax.set_xlim(Rmin, Rlimit)
    ax.set_ylim(ymin, ymax)

    # Adding the temperature axis on top
    x_ticks = np.array([0.05, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    # x_ticks = np.array([0.05, 0.2, 0.4, 0.6])
    T_ticks = midplaneT_profile(R_in, T0, x_ticks).astype(int)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(T_ticks)
    ax2.set_xlabel(r"$T \mathrm{[K]}$", fontsize=10)
    
##    plt.tick_params(bottom=True, top=True, left=True, right=True)
##    plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
##    plt.tick_params(axis='y',which='both',direction='in',length=7,right=True)
##    plt.tick_params(axis='x',which='both',direction='in',length=7,right=True)
##    plt.setp(ax.spines.values(), linewidth=1.5)
##    ax.xaxis.set_tick_params(width=1.5)
##    ax.yaxis.set_tick_params(width=1.5)
    
    leg = plt.legend(loc='upper right', fontsize=9, fancybox=True, handlelength=0.5, prop={'size':csize}, ncol=3)       # Legend properties
    for color, text in zip(colors, leg.get_texts()):
      text.set_color(color)
      text.set_size(10)
      
    plt.tight_layout()
    plt.savefig(pp,format='pdf')
    plt.clf()
    pp.close()


if __name__ == "__main__":
    main()
