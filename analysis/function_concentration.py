import settings
import utils
from netCDF4 import Dataset
import numpy as np
import scipy.stats as stats


def depth_concentration(w_10, w_rise, diffusion_type, boundary, alpha, print_removed=False):
    parcels_file = utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary=boundary, mld=settings.MLD,
                                                 alpha=alpha)
    dataset = Dataset(parcels_file)
    time = dataset.variables['time'][0, :]
    # Bins, at 0.5 meter intervals
    depth_bins = np.arange(0, settings.max_depth, 0.2)
    # Saving the depth profiles
    output_dir = {}
    for t in range(len(time)):
        depth = dataset.variables['z'][:, t]
        # The concentration
        concentrations, bin_edges = np.histogram(depth, bins=depth_bins)
        output_dir[t] = concentrations
        if print_removed:
            removed_frac = np.sum(depth.mask) / depth.shape[0] * 100.
            print(removed_frac)
    output_dir['bin_edges'] = bin_edges
    output_dir['last_time_slice'] = t

    utils.save_obj(filename=utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary, alpha=alpha),
                   item=output_dir)
    utils.remove_file(conduct=True, file_name=parcels_file)