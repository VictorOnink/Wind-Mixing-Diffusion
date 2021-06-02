# *Wind Mixing Diffusion* Repository

## Overview

This repository contains all the code for the vertical wind mixing experiments described in Onink et al. (in prep). Using
the [parcels](http://oceanparcels.org/) v2.2.1 (**P**robably **A** **R**eally **C**omputationally **E**fficient **L**agrangian **S**imulator) package, a 1D
model has been developed to model the vertical concentration profile of buoyant particles in the ocean surface mixed layer depending on the surface wind
stress, the particle rise velocities, and the parametrization of vertical turbulent mixing.

## Model Components
### Particle rise velocities

The particle rise velocities are assigned by the user in the main.py file. The corresponding spherical particle size
for a given rise velocity can be computed based on [Enders et al. (2015)](https://doi.org/10.1016/j.marpolbul.2015.09.027)

### Vertical diffusion profile

Two approaches have been taken for the vertical diffusion profile, and in the code the diffusion parametrizations are
referred as *SWB* and *KPP* diffusion

- **SWB**: The **S**urface **W**ave **B**reaking parametrization is a combination of the wind-driven surface turbulent mixing parametrization
  from [Kukulka et al. (2012)](https://doi.org/10.1029/2012GL051116) and the work of [Poulain (2020)](https://catalogue-archipel.univ-toulouse.fr/primo-explore/fulldisplay?docid=33TOUL_ALMA51536870940004116&context=L&vid=33UT1_VU1&lang=fr_FR&search_scope=default_scope&adaptor=Local%20Search%20Engine&tab=default_tab&query=any,contains,marie%20poulain%20zarcos&offset=0) regarding the depth dependence of the eddy viscosity.
- **KPP**: The K-Profile Parametrization (KPP) approach was first published
  by [Large et al. (1994)](https://doi.org/10.1029/94RG01872) and is commonly used in large-scale ocean models. For the
  1D case, we consider a simpler approach that neglects convective fluxes presented
  by [Boufadel et al. (2020)](https://doi.org/10.1029/2019JC015727).

### Stochastic transport modelling

To model stochastic particle transport, we consider two stochastic transport formulations:
- Markov - 0 (M-0): In this formulation, diffusion is implemented as a random particle displacement, and total stochastic transport is determined by this random walk together with the particle rise velocity. We use the M-0 formulation described by [Grawe et al. (2012)](https://doi.org/10.1007/s10236-012-0523-y).
- Markov - 1 (M-1): In this formulation, a degree of autocorrelation set by parameter $\alpha$ is assumed in the particle velocity perturbations, and the total stochastic transport is determined by the turbulent velocity perturbation and the particle rise velocity. We follow the M-1 formulation described by [Mofakham & Ahmadi (2020)](https://doi.org/10.1016/j.ijmultiphaseflow.2019.103157)
In both cases, the stochastic transport is discretized using an Euler-Maruyama formulation [(Maruyama, 1955)](https://doi.org/10.1007/BF02846028)

### Boundary conditions
The model domain is $z \in [0, 100]$ m, where the vertical axis $z$ is negative upwards with z=0 indicating the ocean surface (although in the figures and the paper we follow the convention that $z$ is positive upward with z=0 at the ocean surface). We have tested a number of boundary conditions (BCs) for the ocean surface:
- Ceiling BC: If a particle crosses the ocean surface (z=0), the particle depth is set to depth = 0.
- Reflect BC: A particle is reflected off of the ocean surface, where for a given overshoot *\partial z* of
  the ocean surface, the final particle position is given as |*\partial z*|.
- Mixed layer boundary condition: If the particle is within a distance *L* of the ocean surface, the particle particle
  depth is set at *R L*, where *R* is a random number between 0 and 1. The mixing depth $L$ is set based on the maximum
  depth at which a particle displacement for a timestep would still be able to cross the ocean surface.
- Reduce dt: If a particle is to cross the ocean surface, then the time step halved if the timestep is greater than a
  minimum timestep. If so, the particle position is recalculated using the reduced timestep. If the halved timestep is
  less than the minimum timestep, a reflecting boundary condition is applied.
If a particle were to cross the lower domain boundary, then we apply a reflecting boundary condition. However, since the maximum depth z=100 m is significantly deeper than the Mixed Layer Depth (MLD), particles in practice never reached the bottom of the model domain.

### Field measurements

In order to verify the modelled distributions, field data of vertical plastic concentration profiles have been collected
from a number of sources:

- [Kukulka et al. (2012)](https://doi.org/10.1029/2012GL051116): Data shared by Tobias Kukulka gave concentration depth
  profiles measured using trawls at a number of stations in the North Atlantic. The shared data also contained the
  surface wind speed and the mixed layer depth (MLD) determined from CTD using
  the [de Boyer Montegut et al. (2004)](https://doi.org/10.1029/2004JC002378) approach.
- [Pieper et al. (2020)](https://doi.org/10.1007/978-3-030-45909-3_21): Depth profiles measured during the PE442
  cruise (August 2018) on the RV Pelagia between Terceira (the Azores) and Catania (Sicily). Unlike the other data, the
  concentrations were measured from Niskin bottles, and the data sampled much greater depths (up to 1150m) below the
  surface, but only limited data close to the surface. Also available was the wind speed at the time of the measurements
  and CTD profiles for almost all statons.
- [Kooi et al., 2016](https://doi.org/10.1038/srep33882): This is depth profiles measured for the first 5m of the water
  column, where the measurements were taken with a multi-depth manta trawl. Data is available for a range of wind
  conditions, and for each station CTD data is available for the computation of the MLD. The data is publically available [here](https://figshare.com/articles/dataset/Data_from_The_effect_of_particle_properties_on_the_depth_profile_of_buoyant_plastics_in_the_ocean_/3427862)
- Amaral-Zettler (unpublished data): This is data collected during the PE448 cruise on the RV Pelagia in the South Atlantic (January 2019)
  using a multi-stage net for sub-surface measurements and a manta trawl for the surface measurements. Wind speed data
  is available, as well as CTD data for determining the MLD.
- [Egger et al. (2020)](https://doi.org/10.1038/s41598-020-64465-8): This is data collected at various depths with a multi-stage neuston net in the North Pacific. Unlike the data in the paper, this data was provided by Matthias Egger without the depth correction that had been applied in their analysis.
With the exception of the [Kooi et al. (2016)](https://doi.org/10.1038/srep33882) data, these data sets are not currently publically available and we recommend contacting the corresponding authors of the studies if one is interested in acquiring this data.

## Code setup 
### Running the model
#### Setting up the run parameters
All the commands for running the model go via the ```main.py``` file. At the beginning of the file, you can set the wind speeds ```w_10```, the rise velocities ```w_rise``` and M-1 memory terms ```alpha``` for which you want to run the simulation (```alpha``` values are not required to run a M-0 simulation, but ```alpha``` must still be defined). All these run parameters are written as lists, since the code will run through all elements of the list to carry out the requested simulations.

The variable ```diffusion``` sets the diffusion profile that will be used in the simulation, and must be set either as ```'KPP'``` or ```'SWB'```. The boundary condition is set by ```boundary```, where the options are either ```'Ceiling'```, ```'Reflect'```, ```'Mixed'``` or ```'Reduce_dt'``` (the BC condition used for the results in the paper is typically ```'Ceiling'``` unless otherwise stated). If you want to run a M-1 simulation, then you add ```'_Markov'``` to the ```diffusion``` string (e.g. ```boundary='Ceiling_Markov'``` would run a M-1 simulation with the Ceiling BC). Otherwise, the default is to run a M-0 simulation.

A number of run parameters (e.g. integration and model output timesteps, physical constants, particle densities, the MLD) are defined in the ```settings.py``` file. These parameters can not be changed from the ```main.py``` file. However, prior to carrying out any simulation on your local computer, please update ```hostname``` to the hostname of your local computer.

#### Carrying out the simulations
Running the ```main.py``` file will run the simulations. First, the code will run ```field_data_processing()```, which will check if the standardized field data files exist and print some basic statistical properties of the field data. If these files are not present, it will create these standardized files. Unless you are in possession of the original field data files, please comment out ```field_data_processing()``` prior to running any simulations.

Next, the function ```ubelix_submission.ubelix_synchronization(update=False)``` uses ```rsync``` to update the files on the laptop of the user with any potential output files on the [Ubelix HPC cluster](https://hpc-unibe-ch.github.io/) from the Universit\"{a}t Bern. Unless one has their own Ubelix account, please keep ```update=False``` so that the function doesn't run. Otherwise, adapt ```ubelix_submission.py``` to work on your own computer.

Next, for all the wind, rise velocity and alpha values in ```w_10```, ```w_rise``` and ```alpha```, the function ```parcels_simulations(wind=wind, rise=rise, alpha=alpha_val)``` is called. Within ```parcels_simulations()```, if one wishes to overwrite any prior model output that exists for the same set of run parameters, set ```conduct = True``` in ```utils.remove_file(conduct=True, file_name=concentration_file)```. Next, if an output file doesn't already exist the parcels simulation for the given parameters will be carried out, and the parcels output will be converted to a non-normalized concentration profile by ```analysis.depth_concentration()```. Once the concentration profile is calculated, we delete the original parcels output to save storage, but if you wish to keep the original parcels file for further analysis set ```remove=False```. Note: this procedure assumes that the simulations are carried out on a local computer, and not the Ubelix cluster. On the Ubelix cluster, the overall steps are the same, but the simulations are submitted as individual jobs using SLURM with the function ```ubelix_submission.ubelix_submission(diffusion, boundary, wind, rise, alpha, submission='parcels')```.

Finally, the function ```plotting()``` calls all plotting functions that are set within ```plotting()```. The visualization functions are all contained within the ```visualization``` directory, and a brief explanation of each plotting function is given either within the ```plotting()``` function or otherwise in the function documentation of the respective visualization functions.


### Overview of files within the respository
- ```main.py```: This is the main file used to define and run the simulations
- ```settings.py```: This is a file used to set a number of constants, dictionaries and directory paths that are used throughout the functions in the directory.
- ```parcels_simulation_functions.py```: This file contains all the functions to run the parcels simulations.
- ```eulerian_simulation_functions.py```: This file contains functions to run Eulerian model simulations. These were ultimately not used in the paper.
- ```field_data.py```: This contains all the functions to convert the original field data files into one standardized format.
- ```ubelix_submission.py```: This contains functions used to run the parcels simulations on the Ubelix server.
- ```analysis```: This directory contains all the analysis functions, split into three files:
  - ```function_concentration.py```: function for converting the parcels model output to concentration profiles.
  - ```function_RMSE.py```: all analysis functions related to calculating RMSE values between field data and concentration profiles
  - ```function_statistics.py```: analysis functions related to calculating statistics of the field data
- ```utils```: This directory contains a number of utility functions, defined as useful functions that are used at multiple points through the model code. The functions are split into three files:
  - ```utils_filenames.py```: These functions are used to name the parcels output, the concentration profile output and the standardized field data
  - ```utils_files.py```: These functions are used for a number of tasks involving files, such as checking if a file exists and saving/loading data using the ```pickle``` module.
  - ```utils_physics.py```: These functions are used to determine physical quantities, such as diffusion profiles, drag coefficients, surface wind stress, etc.
- ```visualization```: This directory contains all the files to visualize the data and create figures. Each file contains the necessary code to create one specific figure described in each function documentation, so please refer there for descriptions of what figure each function creates. The file ```utils_visualization.py``` contains a number of utility functions for creating the figures, such as functions for creating standardized ```matplotlib.pyplot``` axes, determining linestyles and loading concentration profile data such that it can be plotted.

### Installation requirements
In order to run the model, a number of packages need to be installed. I've listed the versions that I used to run all my simulations, and in most cases its probably fine to run with newer package versions. However, this is not guaranteed to always be the case.
- parcels: v2.2.1 or higher. For the installation requirements, please consult the section _Installing Parcels_ at [oceanparcels.org](https://oceanparcels.org/#installing). I did all my simulations with v2.2.1, and generally then the code should also work for more recent versions.
- matplotlib: v3.3.4
- progressbar: v3.37.1
- numpy: v1.19.1
- seabird: v0.11.5: This was only used to analysis CTD data to standardize the field data, so if you are not running this part of the code you don't require this package
- scipy: v1.5.2
- pandas: v1.1.3
- netCDF4: 1.5.4
