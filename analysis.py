# A number of analysis functions for determining the e-folding scale of the vertical distribution of particles
import settings
import utils
from netCDF4 import Dataset
import numpy as np
import math


def depth_concentration(k_z, w_10, w_rise):
    dataset = Dataset(utils.get_parcels_output_name(k_z, w_10, w_rise))
    depth = dataset.variables['z'][:, -1]
    # Bins, at 0.5 meter intervals
    depth_bins = np.arange(0, 100, 0.1)
    # The concentration
    concentrations, bin_edges = np.histogram(depth, bins=depth_bins)
    # Saving the depth profiles
    output_dir = {'concentration': concentrations, 'bin_edges': bin_edges}
    utils.save_obj(filename=utils.get_concentration_output_name(k_z, w_10, w_rise), object=output_dir)


