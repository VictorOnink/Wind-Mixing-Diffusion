import numpy as np
import scipy.stats as stats
import analysis, utils


def correlation_depth_concentration(conduct=True, exclude=None, get_r=[]):
    """
    Calculating the correlation between the particle concentration and the particle depth
    :param exclude: Specifying any field data we don't want to include in the analysis
    :param get_r: list containing field data sources where we want the correlation for just that specific data source
    :return:
    """
    if conduct:
        sources = ['Kooi', 'Pieper', 'Zettler', 'Kukulka', 'Egger']
        sources = utils.exclude_field_data(exclude, sources)

        # Loading in the data
        concentration, norm_depth, depth = [], [], []
        for source in sources:
            data_dict = utils.load_obj(utils.get_data_output_name(source))
            concentration += list(data_dict['concentration'])
            norm_depth += list(data_dict['depth_norm'])
            depth += list(data_dict['depth'])
            # If source is in the list of get_r, determine the correlation for just that specific data source
            if source in get_r:
                determine_correlation(data_dict['concentration'], data_dict['depth'], data_dict['depth_norm'],
                                      subset=source)

        # Determine the data using all field data
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


def range_MLD_values(conduct = True, exclude=None):
    """
    Determining some basic statistics on the MLD
    :param exclude: specifying any field data that we don't want to include in the analysis
    :return:
    """
    if conduct:
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


def correlation_field_model_data(w_10, w_rise, diffusion_type, boundary, alpha, conduct=False):
    """
    Computing the correlation between the average concentrations and the modelled distributions
    :param w_10:
    :param w_rise:
    :param diffusion_type:
    :param boundary:
    :param alpha:
    :return:
    """
    if conduct:
        # First, load the concentration array for the given parameters
        conc_dict = utils.load_obj(filename=utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary,
                                                                                alpha=alpha))
        concentration = conc_dict[conc_dict['last_time_slice']]
        concentration = concentration / concentration.sum()
        concentration_depth = (conc_dict['bin_edges'][:-1] + conc_dict['bin_edges'][1:]) / 2

        # Next, load the averaged field data and the standard deviations, which we convert to variance
        data_dict = utils.load_obj(utils.utils_filenames.get_data_output_name('average'))
        mean_field = data_dict['average'][w_10]
        data_depth = data_dict['depth']

        # For each field data point, find the nearest model point
        nearest_point = np.zeros(data_depth.shape, dtype=np.int32)
        for ind, Z in enumerate(data_depth):
            nearest_point[ind] = utils.utils_files.find_nearest_index(concentration_depth, Z)

        # Calculate the correlation between the field data and the model concentrations
        r, p = stats.pearsonr(concentration[nearest_point], mean_field)

        str_format = diffusion_type, boundary, w_rise, w_10, alpha, r, p
        print('For the {} profile with {}, w_r = {}, w_10 = {}, alpha = {}, r = {:.3f} (p = {:.3f})'.format(*str_format))

