from datetime import timedelta
import numpy as np

# directory for the output files of the parcels simulations
output_dir = '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/parcels_output/'
conc_dir = '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/concentration_output/'
figure_dir = '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/Figures/'
data_dir = '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/Data/'

# Timesteps for integration and for printing to the output file
dt_out = timedelta(seconds=3600)
dt_int = timedelta(seconds=30)

# Runtime for the entire simulation
runtime = timedelta(seconds=12*3600)

# Number of particles in a simulation
p_number = 100000
p_start_depth = 0.0                                         # starting depth of the particles
seed = 1

# Some basic physical parameters
rho_w = 1027                                                # density sea water (kg/m^3)
rho_a = 1.22                                                # density air (kg/m^3)
vk = 0.4                                                    # von Karman constant
wave_age = 35                                               # assuming fully developed sea, Kukulka et al. (2012)
g = 9.81                                                    # acceleration due to gravity (m s^-2)
MLD = 20.                                                   # Ocean mixing layer depth (m)
max_depth = 100.                                            # Maximum depth in our two layer model
phi = 0.9                                                   # Stability function in Monin-Obukov boundary layer theory
mu = 1e-3                                                   # dynamic viscosity
rho_p = 920                                                 # density polypropylene (kg/m^3)
latitude = 35 * np.pi / 180                                 # Latitude for the determination
bulk_diff = 3.7e-4                                          # Dianeutral diffusion below MLD (m^2/s) (Ganachaud, 2003)
w_prime = 0.01                                              # Magnitude of initial w_prime (m/s)
