# A number of analysis functions for determining the e-folding scale of the vertical distribution of particles
import settings
import utils
from netCDF4 import Dataset
import numpy as np
import math


def depth_concentration(k_z, w_10, w_rise, diffusion_type, boundary):
    dataset = Dataset(utils.get_parcels_output_name(k_z, w_10, w_rise, diffusion_type, boundary=boundary))
    time = dataset.variables['time'][0, :]
    # Bins, at 0.5 meter intervals
    depth_bins = np.arange(0, settings.MLD, 0.1)
    # Saving the depth profiles
    output_dir = {}
    for t in range(len(time)):
        depth = dataset.variables['z'][:, t]
        # The concentration
        concentrations, bin_edges = np.histogram(depth, bins=depth_bins)
        output_dir[t] = concentrations
    output_dir['bin_edges'] = bin_edges
    output_dir['last_time_slice'] = t
    utils.save_obj(filename=utils.get_concentration_output_name(k_z, w_10, w_rise, diffusion_type, boundary),
                   object=output_dir)


