from datetime import timedelta
import utils

# directory for the output files of the parcels simulations
output_dir = '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/parcels_output/'
conc_dir = '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/concentration_output/'
figure_dir = '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/Figures/'

# Timesteps for integration and for printing to the output file
dt_out = timedelta(seconds=3600)
dt_int = timedelta(seconds=30)

# Runtime for the entire simulation
runtime = timedelta(seconds=168*3600)

# Number of particles in a simulation
p_number = 100000

# Some basic physical parameters
rho_w = 1027                                                # density sea water (kg/m^3)
rho_a = 1.22                                                # density air (kg/m^3)
vk = 0.4                                                    # von Karman constant
wave_age = 35                                               # assuming fully developed sea, Kukulka et al. (2012)
g = 9.81                                                    # acceleration due to gravity (m s^-2)
MLD = 150                                                   # Ocean mixing layer depth (m)
phi = 0.9                                                   # Stability function in Monin-Obukov boundary layer theory
z0 = 1                                                      # roughness scale of turbulence