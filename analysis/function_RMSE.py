import settings
import utils
from netCDF4 import Dataset
import numpy as np
import scipy.stats as stats
import pandas as pd
from copy import deepcopy


def determine_RMSE(w_10, w_rise, diffusion_type, boundary, alpha, exclude=None, output=False):
    # Loading the concentration profile and the depths
    conc_dict = utils.load_obj(filename=utils.get_concentration_output_name(w_10, w_rise, diffusion_type, boundary,
                               alpha=alpha))
    concentration = conc_dict[conc_dict['last_time_slice']]
    concentration = concentration / concentration.sum()
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


def timestep_dependent_RMSE(conduct: bool, diffusion_type: str):
    diffusion_dict = {'KPP': 'KPP', "Kukulka": 'SWB'}
    output_name = settings.data_dir + '{}_M0_RMSE_timestep.xlsx'.format(diffusion_dict[diffusion_type])
    if conduct is True and not utils.check_file_exist(output_name):
        w_10_list = [0.85, 2.4, 4.35, 6.65, 9.3]
        w_rise_list = [-0.03, -0.003, -0.0003]
        dt_list = [30, 15, 1]

        # Creating a pandas dataframe
        core_dataframe = pd.DataFrame(index=w_10_list, columns=w_rise_list)
        RMSE_dict = {1: deepcopy(core_dataframe), 15: deepcopy(core_dataframe), 30: deepcopy(core_dataframe)}
        for dt in dt_list:
            for w_10 in w_10_list:
                for w_rise in w_rise_list:
                    RMSE_dict[dt][w_rise][w_10] = reference_RMSE_difference(w_rise, w_10, dt, diffusion_type)

        # Saving everything to the output excel file
        writer = pd.ExcelWriter(output_name)
        for sheet in RMSE_dict.keys():
            RMSE_dict[sheet].to_excel(writer, sheet_name='dt={}'.format(sheet))
        writer.save()


def reference_RMSE_difference(w_rise, w_10, dt, diffusion_type, boundary='Reflect'):
    # Getting the reference concentration
    ref_file = utils.get_concentration_output_name(w_rise=w_rise, w_10=w_10, boundary=boundary, dt=1,
                                                   diffusion_type=diffusion_type)
    reference_concentration = utils.load_obj(ref_file)
    reference_concentration = reference_concentration[reference_concentration['last_time_slice']]
    reference_concentration = reference_concentration / reference_concentration.sum()

    # Loading the concentration that we are checking
    candidate_file = utils.get_concentration_output_name(w_rise=w_rise, w_10=w_10, boundary=boundary,
                                                         dt=dt, diffusion_type=diffusion_type)
    candidate_concentration = utils.load_obj(candidate_file)
    candidate_concentration = candidate_concentration[candidate_concentration['last_time_slice']]
    candidate_concentration = candidate_concentration / candidate_concentration.sum()

    # Get the depth, so we can just select the concentrations within the MLD
    depth = utils.load_obj(ref_file)['bin_edges'][:-1]
    within_ML = depth < settings.MLD

    # Computing the RMSE between the two
    RMSE = RMSE_calculation(reference_concentration[within_ML],
                            candidate_concentration[within_ML])
    return RMSE


def RMSE_calculation(reference_array, candidate_array):
    """
    This works assuming the reference and candidate array are both of the same size and are aligned with regards
    to the depth levels.
    """
    return np.sqrt(np.sum(np.square(reference_array - candidate_array)) / reference_array.size)