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

5. `HotStar*` contains the input files, star and disk property files, GGchem output file and output plots such as Q-curves, flux maps, spectra and correlated fluxes for a Herbig Ae type star with varying q and p values to determine the optimal parameters to obtain a well-resolved disk. 

4. `Amorphous_500K` contains the output plots for a HotStar where the opacity values lesser than 500K are taken from the amorphous Olivine and Pyroxene opacities.

6. `Gas_scale_height*` contains the output plots assuming the gas scale height formula used to calculate the surface density. This results in very high optical depths of the order $10^8$. Therefore, some plots are created by multiplying the optical depths with $10^{-6}$ and $10^{-8}$ just to experience what it's like to have the right plots. But of course, all the plots in these folders are wrong. The other folders assume a gas scale height of 1cm. 

7. `Qcurve_inputs` contains the opacity files for several condensates, used to extract the wavelength values and $\kappa_{abs}$.

The scripts are as follows:

   * `all_solids.py` is a small script to obtain all the condensates formed in a GGchem simulation automatically, without the need to input the mineral list. Note that all the condensates obtained may not be plot in the abundance curve since some of them are less than the minimum required abundance.

   * `clean_plots.py` plots the abundance curves for all the indicated condensates with (relatively) clean formatting as a function of temperature.
    
   * `diskprop.py` calculates parameters such as temperature, pressure, density, number density for a planet forming disk model as indicated in Jorge et al. (2022) and plots the radial distribution of these parameters. Change the star parameters in this file to obtain the input properties for the GGchem simulation.
   
   * `fancy_name.py` is a small script to create LaTeX friendly names of the condensate formulae.
   
   * `main.py` is where it all comes together: the output file is loaded, quantities initialized, properties calculated and overall plots created.
   
   * `no_thoughts_just_plots.py` contains plotting functions for individual Qcurves, spectral densities, optical depth maps, flux maps and spectra.
    
   * `radial_plot.py` plots the chemical abundances as a function of radius as well as the temperature (on the top axis).
   
   * `spectra.py` calculates the opacities, flux maps, correlated fluxes, spectra and other dependent quantities.
   
   * `top5_minerals.py` finds the top 5 most abundant condensates at the various radius points.
   
   * `surf_dens_involved.py` calculates the dust scale height (and all the relevant properties therein) to see if the flux values improve by considering the dust scale height in the surface density calculations in `spectra.py` instead of the gas scale height. Regrettably, it does not help. Not even a little. 

### License 

This repository is distributed by the GNU General Public Licences v3.0. You can read more about it [here](https://github.com/drkndl/Silicates-in-Disks/blob/main/LICENSE).
