import numpy as np
import scipy.stats as stats
import analysis
import utils
import settings


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


def correlation_field_model_data(w_10, w_rise, diffusion_type, boundary, alpha, theta, conduct=True, to_print=True,
                                 wave_roughness=False):
    """
    Computing the correlation between the average concentrations and the modelled distributions
    :param w_10:
    :param w_rise:
    :param diffusion_type:
    :param boundary:
    :param alpha:
    :param theta: Langmuir circulation amplification term
    :param conduct: if True, run the code
    :param to_print: if True, print the r and p values, otherwise return the r and p values
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :return:
    """
    if conduct:
        # First, load the concentration array for the given parameters
        conc_dict = utils.load_obj(filename=utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary,
                                                                                alpha=alpha, theta=theta,
                                                                                wave_roughness=wave_roughness))
        concentration = conc_dict['mean_profile']
        concentration = concentration / concentration.sum()
        concentration_depth = conc_dict['bin_edges'][:-1]

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

        if to_print:
            str_format = diffusion_type, boundary, w_rise, w_10, alpha, r, p
            print('For the {} profile with {}, w_r = {}, w_10 = {}, alpha = {}, r = {:.3f} (p = {:.3f})'.format(*str_format))
        else:
            return r, p


def rise_velocity_turbulence_ratio(wind, rise, diffusion, theta=1, conduct=False, wave_roughness=False):
    """
    This calculates the max w' within each depth profile (calculated using equation 4 in the manuscript), and divides
    it by the rise velocity
    :param wind:
    :param rise:
    :param diffusion:
    :param theta: Langmuir circulation amplification term
    :param conduct:
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :return:
    """
    if conduct:
        rise = np.abs(rise)
        # Computing the profiles for the diffusion and diffusion gradient
        depth = np.linspace(0, settings.max_depth, num=settings.depth_levels)
        Kz_profile = utils.get_vertical_diffusion_profile(wind, depth, diffusion, theta=theta,
                                                          wave_roughness=wave_roughness)
        dKz_profile = utils.get_vertical_diffusion_gradient_profile(wind, depth, diffusion, theta=theta,
                                                                    wave_roughness=wave_roughness)

        # Calculating the w' profile, and the maximum value
        dt = settings.dt_int.total_seconds()
        w_prime = dKz_profile + 1 / dt * np.sqrt(2 * Kz_profile) * np.sqrt(dt * 3)
        w_prime_max = np.nanmax(w_prime)

        # Print statement!
        str_format = wind, rise, diffusion, theta, w_prime_max, rise / w_prime_max
        print('For w_10 = {}, rise = {} and diffusion {} with theta {}, we have max w prime {:.5f} and ratio {:.5f}'.format(*str_format))


def mean_percentage_deviation(wind, rise, diffusion, theta=1, conduct=False, wave_roughness=False, boundary='Ceiling',
                              alpha=0.0, to_print=True):
    """
    The RMSE can be skewed by large deviations at the surface where concentrations are significantly higher than at
    depth. Therefore, we will calculate the average percentage difference in
    :param wind:
    :param rise:
    :param diffusion:
    :param theta:
    :param conduct:
    :param wave_roughness:
    :param boundary
    :param alpha
    :param to_print
    :return:
    """
    if conduct:
        # First, load the concentration array for the given parameters
        conc_dict = utils.load_obj(filename=utils.get_concentration_output_name(wind, rise, diffusion,
                                                                                boundary=boundary,
                                                                                alpha=alpha, theta=theta,
                                                                                wave_roughness=wave_roughness))
        concentration = conc_dict['mean_profile']
        concentration = concentration / concentration.sum()
        concentration_depth = conc_dict['bin_edges'][1:]

        # Next, load the averaged field data
        data_dict = utils.load_obj(utils.utils_filenames.get_data_output_name('average'))
        mean_field = data_dict['average'][wind]
        data_depth = data_dict['depth']

        # Just keep the field data above the MLD
        depth_selection = data_depth > -1 * settings.MLD
        mean_field = mean_field[depth_selection]
        data_depth = data_depth[depth_selection]

        # For each field data point, find the nearest model point
        nearest_point = np.array([], dtype=np.int32)
        field_data = np.array([], dtype=np.int32)
        for ind, Z in enumerate(data_depth):
            if mean_field[ind] > 0:
                field_data = np.append(field_data, mean_field[ind])
                nearest_point = np.append(nearest_point, utils.utils_files.find_nearest_index(concentration_depth, Z))

        # Looping through the data and finding the percentage difference relative to the field data
        percent_diff = np.zeros(data_depth.shape, dtype=np.int32)
        abs_percent_diff = np.zeros(data_depth.shape, dtype=np.int32)
        for index, field_point in enumerate(field_data):
            percent_diff[index] = (field_point - concentration[nearest_point[index]]) / field_point * 100.0
            abs_percent_diff[index] = np.abs(percent_diff[index])

        # Calculating the mean percentage difference and mean absolute percentage difference
        mean_percent_diff = np.nanmean(percent_diff)
        std_percent_diff = np.nanstd(percent_diff)
        mean_abs_percent_diff = np.nanmean(abs_percent_diff)
        std_abs_percent_diff = np.nanstd(abs_percent_diff)

        if to_print:
            str_format = diffusion, boundary, rise, wind, alpha, mean_percent_diff, std_percent_diff
            print('For the {} profile with {}, w_r = {}, w_10 = {}, alpha = {}, mean diff = {:.3f}±{:.3f}'.format(*str_format))
            str_format = diffusion, boundary, rise, wind, alpha, mean_abs_percent_diff, std_abs_percent_diff
            print('For the {} profile with {}, w_r = {}, w_10 = {}, alpha = {}, |mean diff| = {:.3f}±{:.3f}'.format(*str_format))
        else:
            return mean_percent_diff, mean_abs_percent_diff


