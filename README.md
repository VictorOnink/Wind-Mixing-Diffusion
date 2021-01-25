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

### Field measurements
In order to verify 


## Code setup
The following files
