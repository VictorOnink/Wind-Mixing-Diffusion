# *Wind Mixing Diffusion* Repository

## Overview
This repository contains all the code for the vertical wind mixing experiments. Using the [parcels](http://oceanparcels.org/) (Probably A Really Computationally Efficient Lagrangian Simulator) package, a 1D model has been developed to model the vertical concentration profile of buoyant particles dependent on the surface wind stress, the particle rise velocities, and the parametrization of the vertical diffusion coefficient profile. 

### Particle rise velocities
The particle rise velocities are assigned by the user in the main.py file. The corresponding elliptical particle size for a given rise velocity can be computed based on [Poulain et al. (2019)](https://doi.org/10.1021/acs.est.8b05458)

### Vertical diffusion profile
Two approaches have been taken for the vertical diffusion profile, and in the code the diffusion parametrizations are referred as *Kukulka* and *KPP* diffusion
- **Kukulka**: This is a combination of the wind-driven surface turbulent mixing parametrization from [Kukulka et al. (2012)](https://doi.org/10.1029/2012GL051116) and the work presented by Marie Poulain-Zarcos at MICRO2020 regarding the depth dependence of the eddy viscosity.
- **KPP**: The K-Profile Parametrization (KPP) approach was first published by [Large et al. (1994)](https://doi.org/10.1029/94RG01872) and is commonly used in large-scale ocean models. For the 1D case, we consider a simpler approach that neglects convective fluxes presented by [Boufadel et al. (2020)](https://doi.org/10.1029/2019JC015727).

### Stochastic transport modelling
To model stochastic particle transport, we have considered both Markov 0 and Markov 1 approaches, where the main difference between the two is that the Markov 0 approach implements diffusion a random particle displacement on top of the displacement due to the particle rise velocity, while the Markov 1 approach implements a diffusion as a random velocity perturbation (with a memory of preceeding random perturbations). The Markov 1 approach requires using a Lagrangian integral timescale, where that is computed based on [Denman & Gargett (1983)](https://doi.org/10.4319/lo.1983.28.5.0801), while the equation is based on [Koszalka et al. (2013)](http://dx.doi.org/10.1016/j.dsr2.2012.07.035). The stochastic vertical transport equation for the Markov 0 approach is from ([Ross & Sharples, 2004](https://doi.org/10.4319/lom.2004.2.289)).

### Boundary conditions
A number of different boundary conditions have been tried for dealing with particles at the ocean/atmosphere interface. These include:
- Zero ceiling boundary condition: If a particle tries to cross the ocean surface, the particle depth is set to 0 (This approach was tested but has been removed from the code due to poor performance. However, the code is still present in the repository history).
- Reflecting boundary condition: A particle is reflected off of the ocean surface, where for a given overshoot *x* of the ocean surface, the final particle position is given as |*x*|.
- Mixed layer boundary condition: If the particle is within a distance *L* of the ocean surface, the particle particle depth is set at *R L*, where *R* is a random number between 0 and 1. The mixing depth $L$ is set based on the maximum depth at which a particle displacement for a timestep would still be able to cross the ocean surface.
- Reduce dt: If a particle is to cross the ocean surface, then the time step halved if the timestep is greater than a minimum timestep. If so, the particle position is recalculated using the reduced timestep. If the halved timestep is less than the minimum timestep, a reflecting boundary condition is applied.

### Field measurements
In order to verify the modelled distributions, field data of vertical plastic concentration profiles have been collected from a number of sources:
- [Kukulka et al. (2012)](https://doi.org/10.1029/2012GL051116): Data shared by Tobias Kukulka gave concentration depth profiles measured using trawls at a number of stations in the North Atlantic. The shared data also contained the surface wind speed and the mixed layer depth (MLD) determined from CTD using the [de Boyer Montegut et al. (2004)](https://doi.org/10.1029/2004JC002378) approach.
- [Pieper et al. (2020)](https://doi.org/10.1007/978-3-030-45909-3_21): Depth profiles measured during the PE442 cruise (August 2018) on the RV Pelagia between Terceira (the Azores) and Catania (Sicily). Unlike the other data, the concentrations were measured from Niskin bottles, and the data sampled much greater depths (up to 1150m) below the surface, but only limited data close to the surface. Also available was the wind speed at the time of the measurements and CTD profiles for almost all statons.
- [Kooi et al., 2016](https://doi.org/10.1038/srep33882): This is depth profiles measured for the first 5m of the water column, where the measurements were taken with a multi-depth manta trawl. Data is available for a range of wind conditions, and for each station CTD data is available for the computation of the MLD.
- PE448 data: This is data collected during the PE448 cruise on the RV Pelagia in the South Atlantic (January 2019) using a multi-stage net for sub-surface measurements and a manta trawl for the surface measurements. Wind speed data is available, but currently not CTD data for determining the MLD.


## Code setup
The following files are contained within the repository:
- main.py: This is the main file, and running this will first load all the field data and output standardized formats, followed by running parcels simulations for the given parameters. This is then followed by creating figures. For the data standardization and the parcels simulations, the code will check if the output files already exist, and will only run the full code for these stages if the output files are not available.
- settings.py: This file contains all the main model parameters, ranging from setting output/data/figure directories to parcels integration settings (dt, particle number) to basic oceanographic parameters. All other files call to the settings file to load these values.
- field_data.py: This contains the data standardization functions for all the field data. The output of these standardization function is a dictionary containing the measurement depth, wind data, normalized concentration (with reference max concentration along a given station depth profile) and the depth normalized by the MLD (if this data is available).
- parcels_simulation_functions.py: This contains all the functions for running the 1D particle model. To run a simulation, the main.py file calls the vertical_diffusion_run function, which requires the surface wind speed, the particle rise velocity, the diffusion type and the boundary condition (which is also what sets if a simulation is Markov 0 or Markov 1).
- analysis.py: This contains any functions necessary for processing the parcels output.
- visualization.py: This contains all the functions for creating the various figures that main.py calls.
- utils_visualization.py: This contains a number of utility functions used for creating the various figures, such as labels for the figure legends, file names for the final figures and options to add diffusion profiles or field measurements to a given plot axis.
- utils.py: This contains further utility functions, excluding any that are used exclusively for creating figures. Examples include functions to check if a file exists, determining the significant wave height or wind stress for given wind conditions and computing the vertical diffusion/diffusion gradient profiles.
