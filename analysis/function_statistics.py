import settings
import utils
from netCDF4 import Dataset
import numpy as np
import scipy.stats as stats
import analysis


def correlation_depth_concentration(exclude=None, get_r=[]):
    sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka', 'Egger']
    sources = utils.exclude_field_data(exclude, sources)

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
    analysis.determine_correlation(np.array(concentration), np.array(depth), np.array(norm_depth))


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
    print('For concentrations and log of depth for {}, r = {:.2f}, p = {:.2f}, '.format(subset, r, p) + r'r^2 = ' +
          '{}'.format(r ** 2))
    print('For concentrations and log of norm. depth for {}, r = {:.2f}, p = {:.2f}, '.format(subset, r_norm, p_norm)
          + r'r^2 = ' + '{}\n'.format(r_norm ** 2))


def range_MLD_values(exclude=None):
    sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka', 'Egger']
    sources = utils.exclude_field_data(exclude, sources)

    # Loading all the determined MLD depths, and then get the mean, std, max and min value over all data
    MLD = []
    for source in sources:
        data_dict = utils.load_obj(utils.utils_filenames.get_data_output_name(source))
        MLD += list(data_dict['MLD'])
    # Remove all nan values
    MLD = np.array(MLD)[~np.isnan(MLD)]
    Max, Min, STD, Mean = np.max(MLD), np.min(MLD), np.std(MLD), np.mean(MLD)
    print('The average MLD over all field data is {:.2f} ± {:.2f}m, min = {:.2f}m, max = {:.2f}m'.format(Mean, STD, Min,
                                                                                                         Max))
