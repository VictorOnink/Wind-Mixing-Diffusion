# A number of analysis functions for determining the e-folding scale of the vertical distribution of particles
import settings
import utils
from netCDF4 import Dataset
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt


def depth_concentration(w_10, w_rise, diffusion_type, boundary):
    dataset = Dataset(utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary=boundary, mld=settings.MLD))
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
    output_dir['bin_edges'] = bin_edges
    output_dir['last_time_slice'] = t
    utils.save_obj(filename=utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary),
                   object=output_dir)


def correlation_depth_concentration(exclude=None, get_r=[]):
    sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka']
    if exclude is not None:
        for source in exclude:
            sources.remove(source)

    # Loading in the data
    concentration, norm_depth, depth = [], [], []
    for source in sources:
        data_dict = utils.load_obj(utils.get_data_output_name(source))
        concentration += list(data_dict['concentration'])
        norm_depth += list(data_dict['depth_norm'])
        depth += list(data_dict['depth'])
        if source in get_r:
            determine_correlation(data_dict['concentration'], data_dict['depth'], data_dict['depth_norm'],
                                  subset=source)

    # Determine the data using all
    determine_correlation(np.array(concentration), np.array(depth), np.array(norm_depth))


def determine_correlation(concentration, depth, norm_depth, subset='all data', zero_correc=1.0):
    """
    Determining the correlation between concentrations and the natural log of depths. The zero_correc term is to assure
    that the logarithm doesn't blow up when there is zero value in depth
    """
    # Converting to log, and removing all points with nan values
    log_depth = np.log(depth[~np.isnan(norm_depth)] + zero_correc)
    log_norm = np.log(norm_depth[~np.isnan(norm_depth)] + zero_correc)
    concentration = concentration[~np.isnan(norm_depth)]

    # Calculating the correlation coefficient
    r, p = stats.pearsonr(concentration, log_depth)
    r_norm, p_norm = stats.pearsonr(concentration, log_norm)
    print('For concentrations and log of depth for {}, r = {}, p = {}, '.format(subset, r, p) + r'r^2 = ' +
          '{}'.format(r ** 2))
    print('For concentrations and log of norm. depth for {}, r = {}, p = {}, '.format(subset, r_norm, p_norm)
          + r'r^2 = ' + '{}\n'.format(r_norm ** 2))
