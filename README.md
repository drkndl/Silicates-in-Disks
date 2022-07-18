# Silicates-in-Disks

Determining the radial dependence of silicates in protoplanetary disks as a part of LEAPS 2022

### Dependencies

* `Numpy`
* `Matplotlib`
* `Pandas`
* `Molmass`
* `Pyvalem`
* `Astropy`

### Files

1. `Default` contains the input file and output plots for the default configuration of GGchem.

2. `Jorge_2022` contains some scripts to calculate Star1 abundances from the solar abundances as given in [Jorge et al. (2022)](https://arxiv.org/pdf/2202.13920.pdf). For example:

   * `jorge_star_abunds.py` calculates the initial chemical abundances of the elements for the GGchem input file corresponding to Star 1 in Jorge et al. (2022).
   
3. `Sun` contains the input file, output plots such as the Q-curves, flux maps, spectra and the correlated fluxes, and other files corresponding to the solar GGchem simulation.

4. `Jorge_2022/Star1` contains the input file, output plots and other files corresponding to the Star 1 GGchem simulation (Jorge et al. (2022)). Currently, the output files here are not correct.

5. `HotStar*` contains the input files, star and disk property files, GGchem output file and output plots such as Q-curves, flux maps, spectra and correlated fluxes for a Herbig Ae type star with varying q, p and Qr values to determine the optimal parameters to obtain a well-resolved disk. The default values were taken as $q=-0.75$, $p=-1.5$ and $Q_r=-1$.

6. `Amorphous_*` contains the output plots for a HotStar where for disk temperatures lesser than *K, the opacity values are taken from the amorphous Olivine and Pyroxene opacities.

7. `Amorphous_Temps_Comparison` contains the overall flux maps, spectra, spectra with varying radii limits and correlated flux plots for various $T_{am}$ i.e. temperature below which the solids are assumed to be amorphous. These plots are obtained by assuming a single grain size of 0.1 $\mu m$ (therefore using the opacity files from `Qcurve_inputs`) for the disk HD144432. 

8. `Gas_scale_height*` contains the output plots assuming the gas scale height formula used to calculate the surface density. This results in very high optical depths of the order $10^8$. Therefore, some plots are created by multiplying the optical depths with $10^{-6}$ and $10^{-8}$ just to experience what it's like to have the right plots. But of course, all the plots in these folders are wrong. The other folders assume a gas scale height of 1cm. 

9. `Qcurve_inputs` contains the opacity files for several crystalline condensates as well as amorphous olivine and pyroxene, assuming a grain size of 0.1 $\mu m$ and $f_{max}$ equal to 0.7 or 1.0. The files contain wavelength values $\lambda$ and absorption efficiencies Q.

10. `Qcurve_inputs_mult_gs` contains the opacity files for several crystalline condensates assuming different grain size distributions such as $0.1-1.5 \mu m$, $0.1-10 \mu m$, $0.1-20 \mu m$ and $0.1-100 \mu m$. These opacity files are obtained using `optool`, the opacity tool [software](https://github.com/cdominik/optool) ([Dominik, C. et al (2021)](https://ui.adsabs.harvard.edu/abs/2021ascl.soft04010D)). Currently, the files are only for crystalline condensates.

11. `Qcurve_Comparison` contains the Q-curves for the condensates from `Qcurve_inputs_mult_gs` comparing how the opacities vary with the grain size distribution. 

12. `HD144432` contains the star and disk property files, GGchem input and output files and the abundance plot, spectra and correlated fluxes for the disk HD144432, considering multiple grain size distributions (using files from `Qcurve_inputs_mult_GS`). Currently, it only considers crystalline grains.

The scripts are as follows:

   * `all_solids.py` is a small script to obtain all the condensates formed in a GGchem simulation automatically, without the need to input the mineral list. Note that all the condensates obtained may not be plot in the abundance curve since some of them are less than the minimum required abundance.
   
   * `compare_Qcurve.py` is used to compare how the opacity curves vary with changing grain size distributions. The outputs of this script are stored in `Qcurves_Comparison`. 
   
   * `compare_amorf_temps.py` is used to compare how the flux maps, spectra and the correlated fluxes vary with changing the temperature below which the gains are assumed to be amorphous, $T_{am}$. The outputs of this script are stored in `Amorphous_Temps_Comparison`.

   * `T_plot.py` plots the abundance curves for all the indicated condensates with (relatively) clean formatting as a function of temperature.
    
   * `diskprop.py` calculates parameters such as temperature, pressure, density, number density for a planet forming disk model as indicated in Jorge et al. (2022) and plots the radial distribution of these parameters. Change the star parameters in this file to obtain the input properties for the GGchem simulation.
   
   * `fancy_name.py` is a small script to create LaTeX friendly names of the condensate formulae.
   
   * `main.py` is where it all comes together: the output file is loaded, quantities initialized, properties calculated and overall plots created.
   
   * `no_thoughts_just_plots.py` contains plotting functions for individual Qcurves, spectral densities, optical depth maps, flux maps and spectra.
   
   * `optool_check.py` checks whether the forsterite Qcurve obtained using optool (Dominik, C. et al (2021)) matches the forsterite Qcurve from Suto (2006) with a grain size of $0.1 \mu m$ and $f_{max} = 1.0$. At the moment, it doesn't really match.  
    
   * `radial_plot.py` plots the chemical abundances as a function of radius as well as the temperature (on the top axis).
   
   * `recreate_spectra.py` is used to plot the spectra and the correlated fluxes for multiple grain size distributions for a disk from [van Boekel et al. (2005)](https://pure.uva.nl/ws/files/2167353/46459_211024y.pdf). Essentially, it is used to corroborate whether my model fits the observations.
   
   * `spectra.py` calculates the opacities, flux maps, correlated fluxes, spectra and other dependent quantities.
   
   * `top5_minerals.py` finds the top 5 most abundant condensates at the various radius points.
   
   * `surf_dens_involved.py` calculates the dust scale height (and all the relevant properties therein) to see if the flux values improve by considering the dust scale height in the surface density calculations in `spectra.py` instead of the gas scale height. Regrettably, it does not help. Not even a little. 
   
   * `temp.py` is simply a temporary 'playground' script I use to quickly check something.

### License 

This repository is distributed by the GNU General Public Licences v3.0. You can read more about it [here](https://github.com/drkndl/Silicates-in-Disks/blob/main/LICENSE).
