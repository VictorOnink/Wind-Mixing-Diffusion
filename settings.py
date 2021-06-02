from datetime import timedelta
import numpy as np
import socket

# Checking which computer this is all running on

if socket.gethostname() == 'Victors-MBP.home':
    server = 'laptop'
else:
    server = 'ubelix'

# directory for the output files of the parcels simulations
root_direc = {'laptop': '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/',
              'ubelix': r'/storage/homefs/vo18e689/Wind-Mixing/'}
output_dir = root_direc[server] + 'parcels_output/'
eulout_dir = root_direc[server] + 'eulerian_output/'
conc_dir = root_direc[server] + 'concentration_output/'
figure_dir = root_direc[server] + 'Figures/'
data_dir = root_direc[server] + 'Data/'
code_dir = root_direc[server] + 'Codes/Wind-Mixing-Diffusion/'
bin_dir = root_direc[server] + r'Codes/Wind-Mixing-Diffusion/bin/'

# Timesteps for integration and for printing to the output file
dt_out = timedelta(seconds=3600)  # timedelta(seconds=100)
dt_int = timedelta(seconds=30)  # timedelta(seconds=1)

# Runtime for the entire simulation
runtime = timedelta(seconds=12*3600)  # timedelta(seconds=3600)

# Number of particles in a simulation
p_number = 100000
p_start_depth = 0.0                                         # starting depth of the particles
seed = 1

# The number of depth levels used in calculating the Kz profiles
depth_levels = 1000

# Some basic physical parameters
# Density of plastic polymers from Brignac et al. (2019) at https://pubs.acs.org/doi/abs/10.1021/acs.est.9b03561
rho_w = 1027                      # density sea water (kg/m^3)
rho_a = 1.22                      # density air (kg/m^3)
vk = 0.4                          # von Karman constant
wave_age = 35                     # assuming fully developed sea, Kukulka et al. (2012)
g = 9.81                          # acceleration due to gravity (m s^-2)
MLD = 20.                         # Ocean mixing layer depth (m)
max_depth = 100.                  # Maximum depth in our two layer model (m)
phi = 0.9                         # Stability function in Monin-Obukov boundary layer theory (Boufadel et al., 2019)
mu = 1e-3                         # dynamic viscosity
nu = 1.1e-6                       # kinematic viscosity of sea water (Enders et al., 2015)
rho_p_pp = 850                    # density polypropylene (kg/m^3)
rho_p_pe = 980                    # density high density polyethylene (kg/m^3)
bulk_diff = 3e-5                  # Dianeutral diffusion below MLD (m^2/s) (Waterhouse et al., 2014)

