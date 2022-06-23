# Silicates-in-Disks

Determining the radial dependence of silicates in protoplanetary disks as a part of LEAPS 2022

### Dependencies

* `Numpy`
* `Matplotlib`
* `Pandas`
* `Molmass`
* `Pyvalem`

### Files

1. `Default` contains the input file and output plots for the default configuration of GGchem

2. `Jorge_2022` contains the input file, output plots (plotted using tools/plot_logT.py in GGchem) and other post-processing scripts and files for the solar GGchem simulation in [Jorge et al. (2022)](https://arxiv.org/pdf/2202.13920.pdf). Some relevant scripts contained within are as follows:

    * `clean_plots.py` plots the abundance curves for all the indicated condensates with (relatively) clean formatting.
    * `jorge_diskprop.py` calculates the parameters such as pressure, density, number density for a planet forming disk model as indicated in Jorge et al. (2022).
    * `jorge_star_abunds.py` calculates the initial chemical abundances of the elements for the GGchem input file corresponding to Star 1 in Jorge et al. (2022).
    * `opacities.py` calculates the opacities, fluxes and plots the spectra corresponding to the various condensates.
    * `radial_plot.py` plots the chemical abundances as a function of radius instead of the temperature.
    * `top5_minerals.py` finds the top 5 most abundant condensates at the various radius points.

3. `Jorge_2022/Sun` contains the input file, output plots and other files corresponding to the solar GGchem simulation.

4. `Jorge_2022/Star1` contains the input file, output plots and other files corresponding to the Star 1 GGchem simulation (Jorge et al. (2022)). Currently, the output files here are not correct.

### License 

This repository is protected by the GNU General Public Licences v3.0. You can read more about it [here](https://github.com/drkndl/Silicates-in-Disks/blob/main/LICENSE).
