import settings
import utils
from netCDF4 import Dataset
import numpy as np


def depth_concentration(w_10, w_rise, diffusion_type, boundary, alpha, bin_size=0.2, remove_file=True):
    """
    Binning the parcels output for a given run into vertical bins of size bin_size.
    note: the concentrations are not normalized, that is done later when the concentrations are loaded
    :param w_10: 10m wind speed
    :param w_rise: rise velocity
    :param diffusion_type: diffusion type, either SWB or KPP
    :param boundary: boundary condition, and whether M-0 or M-1
    :param alpha: memory term for M-1
    :param remove_file: if True, remove the original parcels output file
    :return:
    """
    # Loading the relevant parcels file
    parcels_file = utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary=boundary, mld=settings.MLD,
                                                 alpha=alpha)
    dataset = Dataset(parcels_file)
    time = dataset.variables['time'][0, :]
    # Setting the depth bins and creating the output directory
    depth_bins = np.arange(0, settings.max_depth, bin_size)
    output_dir = {}
    # Looping through all the saved timesteps in the parcels iutput
    for t in range(len(time)):
        # Loading the depth data for all the particles
        depth = dataset.variables['z'][:, t]
        # Saving the binned concentrations
        concentrations, bin_edges = np.histogram(depth, bins=depth_bins)
        output_dir[t] = concentrations
        # Check if particles have been removed, if yes raise an assertion error
        removed_frac = np.sum(depth.mask) / depth.shape[0] * 100.
        str_format = diffusion_type, boundary, w_10, w_rise, alpha
        assert removed_frac == 0.0, "Particle number was not conserved for {} with {}, (w_10={}, w_rise={}, a={})".format(*str_format)
    # Saving the bin edges, and specifying which time slice was the last for easy loading later
    output_dir['bin_edges'] = bin_edges
    output_dir['last_time_slice'] = t

    # Pickling the concentration array
    utils.save_obj(filename=utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary, alpha=alpha),
                   item=output_dir)
    # We don't need the parcels file for any subsequent analysis, so I'm removing it to save storage on my computer
    if remove_file:
        utils.remove_file(conduct=True, file_name=parcels_file)
