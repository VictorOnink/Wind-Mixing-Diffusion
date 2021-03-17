from datetime import timedelta
import numpy as np

# Computer this is running on. 0 = my laptop, 1 = Ubelix
server_dict = {0: 'laptop', 1: 'ubelix'}
server = server_dict[1]

# directory for the output files of the parcels simulations
root_direc = {'laptop': '/Users/victoronink/Desktop/Bern Projects/Wind Mixing/',
              'ubelix': r'/home/ubelix/climate/vo18e689/Wind-Mixing/'}
output_dir = root_direc[server] + 'parcels_output/'
conc_dir = root_direc[server] + 'concentration_output/'
figure_dir = root_direc[server] + 'Figures/'
data_dir = root_direc[server] + 'Data/'
code_dir = root_direc[server] + 'Codes/Wind-Mixing-Diffusion/'
bin_dir = root_direc[server] + r'Codes/Wind-Mixing-Diffusion/bin/'

# Timesteps for integration and for printing to the output file
dt_out = timedelta(seconds=3600)  # timedelta(seconds=100)
dt_int = timedelta(seconds=1)  # timedelta(seconds=1)

# Runtime for the entire simulation
runtime = timedelta(seconds=6*3600)  # timedelta(seconds=3600)

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
bulk_diff = 3e-5                                            # Dianeutral diffusion below MLD (m^2/s) (Waterhouse et al., 2014)
w_prime = 0.001                                             # Magnitude of initial w_prime (m/s)
alpha = 0.95                                                    # Lagrangian integral timescale as multiple of dt_int
