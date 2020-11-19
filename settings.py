from datetime import timedelta

# directory for the output files of the parcels simulations
output_dir = 'ctoronink/Desktop/Bern Projects/Wind Mixing/parcels_output/'

# Timesteps for integration and for printing to the output file
dt_out = timedelta(days=14)
dt_int = timedelta(seconds=30)

# Runtime for the entire simulation
runtime = timedelta(days=365*1)