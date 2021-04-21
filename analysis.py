# A number of analysis functions for determining the e-folding scale of the vertical distribution of particles
import settings
import utils
from netCDF4 import Dataset
import numpy as np
import scipy.stats as stats


def depth_concentration(w_10, w_rise, diffusion_type, boundary, alpha, print_removed=True):
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


def determine_RMSE(w_10, w_rise, diffusion_type, boundary, alpha, exclude=None, output=False):
    # Loading the concentration profile and the depths
    conc_dict = utils.load_obj(filename=utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary,
                               alpha=alpha))
    concentration = conc_dict[conc_dict['last_time_slice']]
    concentration = concentration / concentration.max()
    concentration_depth = conc_dict['bin_edges'][:-1]

    # Loading the field measurements
    sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka', 'Egger']
    sources = utils.exclude_field_data(exclude, sources)
    field_data, depth, wind = np.array([]), np.array([]), np.array([])
    for source in sources:
        data_dict = utils.load_obj(utils.utils_filenames.get_data_output_name(source))
        field_data = np.append(field_data, data_dict['concentration'])
        wind = np.append(wind, data_dict['wind_speed'])
        depth = np.append(depth, data_dict['depth'])

    # Select only the measurements measured under the selected wind conditions
    w_10_select = {0.85: utils.beaufort_limits()[1], 2.4: utils.beaufort_limits()[2], 4.35: utils.beaufort_limits()[3],
                   6.65: utils.beaufort_limits()[4], 9.3: utils.beaufort_limits()[5]}
    w_min, w_max = w_10_select[w_10]
    wind_select = (w_min < wind) & (wind <= w_max)
    field_data, depth = field_data[wind_select], depth[wind_select]

    # For each field data point, calculate the index of concentration_depth that is closest
    nearest_point = np.zeros(depth.shape, dtype=np.int32)
    for ind, Z in enumerate(depth):
        nearest_point[ind] = utils.utils_files.find_nearest_index(concentration_depth, Z)

    # Now, calculate the RMSE
    RMSE = np.sqrt(np.sum(np.square(field_data - concentration[nearest_point])) / field_data.size)

    # And then printing or returning the RMSE value the output
    if not output:
        str_format = diffusion_type, boundary, w_rise, w_10, alpha, RMSE, field_data.size
        print('For the {} profile with {}, w_r = {}, w_10 = {}, alpha = {}, RMSE = {} over {} data points'.format(
            *str_format))
    elif output:
        return RMSE


def determine_RMSE_eulerian(w_10, w_rise, diffusion_type, exclude=None):
    eul_dict = utils.load_obj(utils.get_eulerian_output_name(w_10=w_10, w_rise=w_rise, diffusion_type=diffusion_type))
    concentration = eul_dict['C']
    concentration_depth = np.abs(eul_dict['Z'])

    # Loading the field measurements
    sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka', 'Egger']
    sources = utils.exclude_field_data(exclude, sources)
    field_data, depth, wind = np.array([]), np.array([]), np.array([])
    for source in sources:
        data_dict = utils.load_obj(utils.utils_filenames.get_data_output_name(source))
        field_data = np.append(field_data, data_dict['concentration'])
        wind = np.append(wind, data_dict['wind_speed'])
        depth = np.append(depth, data_dict['depth'])

    # Select only the measurements measured under the selected wind conditions
    w_10_select = {0.85: utils.beaufort_limits()[1], 2.4: utils.beaufort_limits()[2], 4.35: utils.beaufort_limits()[3],
                   6.65: utils.beaufort_limits()[4], 9.3: utils.beaufort_limits()[5]}
    w_min, w_max = w_10_select[w_10]
    wind_select = (w_min < wind) & (wind <= w_max)
    field_data, depth = field_data[wind_select], depth[wind_select]

    # For each field data point, calculate the index of concentration_depth that is closest
    nearest_point = np.zeros(depth.shape, dtype=np.int32)
    for ind, Z in enumerate(depth):
        nearest_point[ind] = utils.utils_files.find_nearest_index(concentration_depth, Z)

    # Now, calculate the RMSE
    RMSE = np.sqrt(np.sum(np.square(field_data - concentration[nearest_point])) / field_data.size)

    return RMSE